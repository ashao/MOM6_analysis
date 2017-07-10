import numpy as np
import pandas as pd

def interpolate_for_nondim_position(dRhoTop, Pneg, dRhoBot, Ppos):
  if (Ppos<=Pneg):
    interpolate_for_nondim_position = 0.
  elif ( dRhoBot - dRhoTop > 0. ):
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoTop / ( dRhoBot - dRhoTop ) ) )
  elif ( dRhoBot - dRhoTop == 0):
    if (dRhoTop>0.):
      interpolate_for_nondim_position = 0.
    elif (dRhoTop<0.):
      interpolate_for_nondim_position = 1.
    else:
      interpolate_for_nondim_position = 0.5
  else:
    interpolate_for_nondim_position = 0.5
  return interpolate_for_nondim_position

def absolute_position(Pint, Karr, NParr, k_surface):
    k = Karr[k_surface]
    return Pint[k] + NParr[k_surface] * ( Pint[k+1] - Pint[k] )

def interface_nondim_position(k):
    if np.floor(0.5*k) == np.floor(0.5*(k+1)):
      return 0.
    else:
      return 1.

def set_neutral_surface_position( Pres_l, Tint_lt, Tint_lb, Sint_lt, Sint_lb, dRdT_lt, dRdT_lb, dRdS_lt, dRdS_lb, Pres_r, Tint_rt, Tint_rb, Sint_rt, Sint_rb, dRdT_rt, dRdT_rb, dRdS_rt, dRdS_rb ):

  nk = Tint_lt.size
  Pl = np.zeros(2*nk) ; Pr = np.zeros(2*nk)
  Tl = np.zeros(2*nk) ; Tr = np.zeros(2*nk)
  Sl = np.zeros(2*nk) ; Sr = np.zeros(2*nk)
  dRdT_l = np.zeros(2*nk) ; dRdT_r = np.zeros(2*nk) ;
  dRdS_l = np.zeros(2*nk) ; dRdS_r = np.zeros(2*nk) ;

  PoL = np.zeros(4*nk) ; PoR = np.zeros(4*nk)
  PoL_abs = np.zeros(4*nk) ; PoR_abs = np.zeros(4*nk)
  KoL = np.zeros(4*nk, dtype=np.int8) ; KoR = np.zeros(4*nk, dtype=np.int8)
  hL = np.zeros(4*nk) ; hR = np.zeros(4*nk)
  hEff = np.zeros(4*nk-1)

  Pl[::2]  = Pres_l[0:-1] ; Pl[1::2] = Pres_l[1:]
  Sl[::2]  = Sint_lt ; Sl[1::2] = Sint_lb
  Tl[::2]  = Tint_lt ; Tl[1::2] = Tint_lb
  dRdT_l[::2]  = dRdT_lt ; dRdT_l[1::2] = dRdT_lb
  dRdS_l[::2]  = dRdS_lt ; dRdS_l[1::2] = dRdS_lb

  Pr[::2]  = Pres_r[0:-1] ; Pr[1::2] = Pres_r[1:]
  Sr[::2]  = Sint_rt ; Sr[1::2] = Sint_rb
  Tr[::2]  = Tint_rt ; Tr[1::2] = Tint_rb
  dRdT_r[::2]  = dRdT_rt ; dRdT_r[1::2] = dRdT_rb
  dRdS_r[::2]  = dRdS_rt ; dRdS_r[1::2] = dRdS_rb

  kr = 0 ; lastK_right = 0 ; lastP_right = -1.
  kl = 0 ; lastK_left = 0 ; lastP_left = -1.
  reached_bottom = False
  left_column = pd.DataFrame({"Pl" : Pl, "Tl": Tl, "Sl" : Sl, "dRdT_l" : dRdT_l, "dRdS_l" : dRdS_l})
  right_column = pd.DataFrame({"Pr" : Pr, "Tr": Tr, "Sr" : Sr, "dRdT_r" : dRdT_r, "dRdS_r" : dRdS_r})
  print(left_column)
  print(right_column)
  for k_surface in np.arange(0,4*nk):
    klm1 = max(kl-1,0)
    krm1 = max(kr-1,0)
    dRho = 0.5 * ( ( dRdT_r[kr] + dRdT_l[kl] ) * ( Tr[kr] - Tl[kl] )
                  + (dRdS_r[kr] + dRdS_r[kr] ) * ( Sr[kr] - Sl[kl] ) )
    print( "\nWorking on k_surface %d: dRho: %f Tl[%d]: %f Tr[%d]: %f" % (k_surface, dRho, kl, Tl[kl], kr, Tr[kr]))

    if not reached_bottom:
      if dRho < 0.:
        searching_left_column = True
        searching_right_column = False
      elif dRho > 0.:
        searching_right_column = True
        searching_left_column = False
      else:
        if (kl+kr==0):
          searching_left_column = True
          searching_right_column = False
        else:
          searching_left_column = not searching_left_column
          searching_right_column = not searching_right_column

    if searching_left_column:
      same_k = np.floor(0.5*(klm1)) == np.floor(0.5*kl)
      dRhoTop = 0.5 * ( ( dRdT_l[klm1] + dRdT_r[kr] ) * ( Tl[klm1] - Tr[kr] )
                      + ( dRdS_l[klm1] + dRdS_r[kr] ) * ( Sl[klm1] - Sr[kr] ) )
      dRhoBot = 0.5 * ( ( dRdT_l[klm1+1] + dRdT_r[kr] ) * ( Tl[klm1+1] - Tr[kr] )
                      + ( dRdS_l[klm1+1] + dRdS_r[kr] ) * ( Sl[klm1+1] - Sr[kr] ) )
      print("Searching left: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
      print("klm1: %d kl: %d kr: %d" % (klm1,kl,kr))
      # Search left
      if kr + kl == 0.:
        PoL[k_surface] = 0.
        lastP_left = 0.
      elif same_k: # Searching within a layer
        if dRhoTop > 0.: # Right surface always lighter, point to top interface of this layer
          PoL[k_surface] = 0.
          print("Point to top interface")
        elif dRhoTop >= dRhoBot: # Left layer unstratified, point to top unless it was pointed to last
          if lastP_left == 0. or np.floor(0.5*klm1)==(nk-1):
            PoL[k_surface] = 1.
          else:
            PoL[k_surface] = 0.
          print("dRhoTop >= dRhoBot")
        else:
          # dRhoTop is negative and dRhoTop is postive, so interpolate to find neutral surface
          PoL[k_surface] = interpolate_for_nondim_position( dRhoTop, Pl[klm1], dRhoBot, Pl[klm1+1] )
          print("Interpolate for position")
        lastP_left = PoL[k_surface]
      else: # Searching across a discontinuity, similar logic as an existing layer except for interpolation case where we're going to just use a nearest neighbor approach
        PoL[k_surface] = 0.
        klm1 += 1
#        PoL[k_surface] = 1.

        lastP_left = -1.
        print("At discontinuity")

      KoL[k_surface] = np.floor(0.5*klm1)
      KoR[k_surface] = np.floor(0.5*kr)
      if kr <= (2*nk-1):
        PoR[k_surface] = interface_nondim_position(kr)
      else:
        PoR[k_surface] = 1.
        KoR[k_surface] = nk-1
      if kr < 2*nk-1:
          kr = kr + 1
      else:
        reached_bottom = True
        searching_right_column = True
        searching_left_column = False


    elif (searching_right_column):
      same_k = np.floor( krm1*0.5 ) == np.floor( (krm1+1)*0.5 )
      same_p = Pr[krm1] == Pr[krm1+1]
      dRhoTop = 0.5 * ( ( dRdT_r[krm1] + dRdT_l[kl] ) * ( Tr[krm1] - Tl[kl] )
                      + ( dRdS_r[krm1] + dRdS_l[kl] ) * ( Sr[krm1] - Sl[kl] ) )
      dRhoBot = 0.5 * ( ( dRdT_r[krm1+1] + dRdT_l[kl] ) * ( Tr[krm1+1] - Tl[kl] )
                      + ( dRdS_r[krm1+1] + dRdS_l[kl] ) * ( Sr[krm1+1] - Sl[kl] ) )
      print("Searching right: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
      print("krm1: %d kr: %d kl: %d" % (krm1,kr,kl))
      if kr + kl == 0.:
        PoR[k_surface] = 0.
        lastP_right = 0.
        print("At surface")
      elif same_k: # Searching within a layer
        if dRhoTop > 0.:
          PoR[k_surface] = 0.
          print("dRhoTop>0")
        elif dRhoTop >= dRhoBot:
          if lastP_right == 0. or np.floor(0.5*krm1)==(nk-1):
            PoR[k_surface] = 1.
          else:
            PoR[k_surface] = 0.
          print("dRhoTop>=dRhoBot")
        else:
          PoR[k_surface] = interpolate_for_nondim_position( dRhoTop, Pr[krm1], dRhoBot, Pl[krm1+1] )
          print("Interpolate for position")
        lastP_right = PoR[k_surface]
      else:        # Searching across a discontinuity
        PoR[k_surface] = 1.
        lastP_right = -1
        print("At discontinuity")

      KoR[k_surface] = np.floor(0.5*krm1)
      KoL[k_surface] = np.floor(0.5*kl)
      if kl <= (2*nk-1):
        PoL[k_surface] = interface_nondim_position(kl)
      else:
        PoL[k_surface] = 1.
        KoL[k_surface] = nk-1
      if kl < 2*nk-1:
        kl = kl + 1
      else:
        reached_bottom = True
        searching_right_column = False
        searching_left_column = True
    else:
      print("ERROR")

    print("Position on left : %f" % PoL[k_surface])
    print("Position on right: %f" % PoR[k_surface])
    if k_surface>0:
      PoL_abs[k_surface] = absolute_position(Pres_l, KoL, PoL, k_surface)
      PoR_abs[k_surface] = absolute_position(Pres_r, KoR, PoR, k_surface)
      hL[k_surface] = absolute_position(Pres_l, KoL, PoL, k_surface) - absolute_position(Pres_l, KoL, PoL, k_surface-1)
      hR[k_surface] = absolute_position(Pres_r, KoR, PoR, k_surface) - absolute_position(Pres_r, KoR, PoR, k_surface-1)
      if (hL[k_surface] + hR[k_surface] > 0.):
          hEff[k_surface-1] = 2. * hL[k_surface] * hR[k_surface] / ( hL[k_surface]+hR[k_surface] )
      else:
          hEff[k_surface-1] = 0.
      print("hL: %f hR: %f hEff: %f" % (hL[k_surface],hR[k_surface],hEff[k_surface-1]))
  print(PoL,PoR)
  return PoL, PoR, PoL_abs, PoR_abs, KoL, KoR, hEff, hL, hR


def set_neutral_surface_position2( Pres_l, Tint_lt, Tint_lb, Sint_lt, Sint_lb, dRdT_lt, dRdT_lb, dRdS_lt, dRdS_lb, Pres_r, Tint_rt, Tint_rb, Sint_rt, Sint_rb, dRdT_rt, dRdT_rb, dRdS_rt, dRdS_rb ):
    nk = Tint_lt.size
    Pl = np.zeros((nk,2)) ; Pr = np.zeros((nk,2))
    Tl = np.zeros((nk,2)) ; Tr = np.zeros((nk,2))
    Sl = np.zeros((nk,2)) ; Sr = np.zeros((nk,2))
    dRdT_l = np.zeros((nk,2)) ; dRdT_r = np.zeros((nk,2)) ;
    dRdS_l = np.zeros((nk,2)) ; dRdS_r = np.zeros((nk,2)) ;
    search_dir = None
    PoL = np.zeros(4*nk) ; PoR = np.zeros(4*nk)
    PoL_abs = np.zeros(4*nk) ; PoR_abs = np.zeros(4*nk)
    KoL = np.zeros(4*nk, dtype=np.int8) ; KoR = np.zeros(4*nk, dtype=np.int8)
    hL = np.zeros(4*nk) ; hR = np.zeros(4*nk)
    hEff = np.zeros(4*nk-1)

    Pl[:,0] = Pres_l[0:-1]
    Pl[:,1] = Pres_l[1:]
    Sl[:,0] = Sint_lt
    Sl[:,1] = Sint_lb
    Tl[:,0] = Tint_lt
    Tl[:,1] = Tint_lb
    dRdT_l[:,0] = dRdT_lt
    dRdT_l[:,1] = dRdT_lb
    dRdS_l[:,0] = dRdS_lt
    dRdS_l[:,1] = dRdS_lb

    Pr[:,0] = Pres_r[0:-1]
    Pr[:,1] = Pres_r[1:]
    Sr[:,0] = Sint_rt
    Sr[:,1] = Sint_rb
    Tr[:,0] = Tint_rt
    Tr[:,1] = Tint_rb
    dRdT_r[:,0] = dRdT_rt
    dRdT_r[:,1] = dRdT_rb
    dRdS_r[:,0] = dRdS_rt
    dRdS_r[:,1] = dRdS_rb

    k_search_left = 0
    k_search_right = 0
    kl_r = 0 ; ki_r = 0 ; lastK_right = -1 ; lastP_right = -1.
    kl_l = 0 ; ki_l = 0 ; lastK_left = -1 ; lastP_left = -1.
    reached_bottom = False

    for k_surface in np.arange(0,4*nk):

        print( "\nWorking on k_surface %d: Tl[%d,%d]: %f Tr[%d,%d]: %f" % (k_surface, kl_l, ki_l, Tl[kl_l,ki_l], kl_r, ki_r, Tr[kl_r,ki_r]))
        dRho = 0.5 * ( ( dRdT_r[kl_r,ki_r] + dRdT_l[kl_l,ki_l] ) * ( Tr[kl_r,ki_r] - Tl[kl_l,ki_l] )
                      + (dRdS_r[kl_r,ki_r] + dRdS_l[kl_l,ki_l] ) * ( Sr[kl_r,ki_r] - Sl[kl_l,ki_l] ) )
        if not reached_bottom:
          if dRho < 0.:
            searching_left_column = True
            searching_right_column = False
          elif dRho > 0.:
            searching_right_column = True
            searching_left_column = False
          else:
            if (kl_l+kl_r ==0 and ki_l+ki_r==0):
              searching_left_column = True
              searching_right_column = False
            else:
              searching_left_column = not searching_left_column
              searching_right_column = not searching_right_column

        if searching_left_column:
          search_dir = np.append(search_dir,0)
          dRhoTop = 0.5 * ( ( dRdT_l[k_search_left,0] + dRdT_r[kl_r,ki_r] ) * ( Tl[k_search_left,0] - Tr[kl_r,ki_r] )
                          + ( dRdS_l[k_search_left,0] + dRdS_r[kl_r,ki_r] ) * ( Sl[k_search_left,0] - Sr[kl_r,ki_r] ) )
          dRhoBot = 0.5 * ( ( dRdT_l[k_search_left,1] + dRdT_r[kl_r,ki_r] ) * ( Tl[k_search_left,1] - Tr[kl_r,ki_r] )
                          + ( dRdS_l[k_search_left,1] + dRdS_r[kl_r,ki_r] ) * ( Sl[k_search_left,1] - Sr[kl_r,ki_r] ) )
          print("Searching left layer %d: dRhoTop: %f dRhoBot: %f" % (k_search_left, dRhoTop, dRhoBot))
          # Set positions on left
          if kl_r + kl_l + ki_r + ki_l == 0.:
            PoL[k_surface] = 0.
          elif dRhoTop > 0.:
            PoL[k_surface] = 0.
          elif dRhoTop >= dRhoBot: # Perfectly unstratified
            if lastP_left == 0.: # We've already pointed to the top interface so now point to the bottom one
              PoL[k_surface] = 1.
              print("Unstratified, point bottom")
            elif lastP_left == 1.:
              PoL[k_surface] = 1.
              print("Unstratified, point bottom")
            else:             # Point to the top
              PoL[k_surface] = 0.
              print("Unstratified, point to top")
          elif dRhoTop > dRhoBot: # This layer is unstable (bottom of layer is lighter than top). Point to bottom always
            PoL[k_surface] = 1.
          elif dRhoBot > 0. and dRhoTop >0.: # Neutral surface lgihter than anything in layer, point to bottom
            PoL[k_surface] = 0.
          else:
            PoL[k_surface] = interpolate_for_nondim_position( dRhoTop, Pl[k_search_left,0], dRhoBot, Pl[k_search_left,1] )
          KoL[k_surface] = k_search_left
          KoR[k_surface] = kl_r

          lastP_left = PoL[k_surface]
          if PoL[k_surface] == 1. and kl_l < nk-1:
            k_search_left = k_search_left + 1
#            k_search_right = k_search_right + 1
            lastP_left = -1

          # Set position on right
          if ki_r == 0:
            PoR[k_surface] = 0.
          else:
            PoR[k_surface] = 1.

          # Determine which interface on the right will need to be connected
          if ki_r == 0:
            ki_r = 1
          elif kl_r < (nk-1): # Set to top interface of next layer unless at the bottom of the column
            ki_r = 0
            kl_r = kl_r + 1
          elif kl_r == (nk-1) and ki_r == 1: # Reached the bottom
            reached_bottom = True
            searching_right_column = True
            searching_left_column = False
          else:
            print("AHHHH!")

        elif (searching_right_column):
          search_dir = np.append(search_dir,1)
          dRhoTop = 0.5 * ( ( dRdT_r[k_search_right,0] + dRdT_l[kl_l,ki_l] ) * ( Tr[k_search_right,0] - Tl[kl_l,ki_l] )
                          + ( dRdS_r[k_search_right,0] + dRdS_l[kl_l,ki_l] ) * ( Sr[k_search_right,0] - Sl[kl_l,ki_l] ) )
          dRhoBot = 0.5 * ( ( dRdT_r[k_search_right,1] + dRdT_l[kl_l,ki_l] ) * ( Tr[k_search_right,1] - Tl[kl_l,ki_l] )
                          + ( dRdS_r[k_search_right,1] + dRdS_l[kl_l,ki_l] ) * ( Sr[k_search_right,1] - Sl[kl_l,ki_l] ) )
          print("Searching right layer %d: dRhoTop: %f dRhoBot: %f" % (k_search_right, dRhoTop,dRhoBot))
          # Set positions on right
          if k_search_right + ki_r == 0.:
            PoR[k_surface] = 0.
          elif dRhoTop > 0.:
            PoR[k_surface] == 0.
          elif dRhoTop >= dRhoBot: # Perfectly unstratified
            if lastP_right == 0.: # We've already pointed to the top interface so now point to the bottom one
              PoR[k_surface] = 1.
              print("Unstratified, point bottom")
            elif lastP_right == 1.:
              PoR[k_surface] = 1.
              print("Unstratified, point bottom")
            else:             # Point to the top
              PoR[k_surface] = 0.
              print("Unstratified, point top")
          elif dRhoTop > dRhoBot: # This layer is unstable (bottom of layer is lighter than top). Point to bottom always
            PoR[k_surface] = 1.
          elif dRhoBot > 0. and dRhoTop >0.: # Neutral surface lgihter than anything in layer, point to bottom
            PoR[k_surface] = 0.
          else:
            PoR[k_surface] = interpolate_for_nondim_position( dRhoTop, Pr[k_search_right,0], dRhoBot, Pr[k_search_right,1] )

          lastP_right = PoR[k_surface]
          KoL[k_surface] = kl_l
          KoR[k_surface] = k_search_right

          if PoR[k_surface] == 1.:
            k_search_right = k_search_right + 1
            lastP_right = -1

          # Set position on left
          if ki_l == 0:
            PoL[k_surface] = 0.
          else:
            PoL[k_surface] = 1.

          # Determine which interface on the left will need to be connected
          if ki_l == 0: # Move to bottom interface of same layer
            ki_l = 1
          elif kl_l < (nk-1): # Set to top interface of next layer unless at the bottom of the column
            ki_l = 0
            kl_l = kl_l + 1
          elif kl_l == (nk-1) and ki_l == 1: # Reached the bottom
              reached_bottom = True
              searching_right_column = False
              searching_left_column = True
          else:
            print("AHHHH!")
        else:
            print("ERROR")

        print("Position on left : %f" % PoL[k_surface])
        print("Position on right: %f" % PoR[k_surface])
        if k_surface>0:
          print(Pres_l,KoL)
          PoL_abs[k_surface] = absolute_position(Pres_l, KoL, PoL, k_surface)
          PoR_abs[k_surface] = absolute_position(Pres_r, KoR, PoR, k_surface)
          hL[k_surface] = absolute_position(Pres_l, KoL, PoL, k_surface) - absolute_position(Pres_l, KoL, PoL, k_surface-1)
          hR[k_surface] = absolute_position(Pres_r, KoR, PoR, k_surface) - absolute_position(Pres_r, KoR, PoR, k_surface-1)
          if (hL[k_surface] + hR[k_surface] > 0.):
            hEff[k_surface-1] = 2. * hL[k_surface] * hR[k_surface] / ( hL[k_surface]+hR[k_surface] )
          else:
            hEff[k_surface-1] = 0.
          print("hL: %f hR: %f hEff: %f" % (hL[k_surface],hR[k_surface],hEff[k_surface-1]))
    return PoL, PoR, PoL_abs, PoR_abs, KoL, KoR, hEff, hL, hR, search_dir

def find_neutral_surface_positions_continuous(Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, dRdTr, dRdSr):
  nk = Pl.size - 1
  kr = 0 ; lastK_right = 0 ; lastP_right = 0.
  kl = 0 ; lastK_left = 0 ; lastP_left = 0.
  reached_bottom = False

  PoL = np.zeros(2*nk+2) ; PoR = np.zeros(2*nk+2)
  PoL_abs = np.zeros(2*nk+2) ; PoR_abs = np.zeros(2*nk+2)
  KoL = np.zeros(2*nk+2, dtype = np.int8) ; KoR = np.zeros(2*nk+2, dtype = np.int8)
  hEff = np.zeros(2*nk+1)
  hL = np.zeros(2*nk+2)
  hR = np.zeros(2*nk+2)


  # Loop over each neutral surface, working from top to bottom
  for k_surface in np.arange(0, 2*nk+2):
    klm1 = max(kl-1, 0)
    krm1 = max(kr-1, 0)

    print( "\nWorking on k_surface %d: Tl[%d]: %f Tr[%d]: %f" % (k_surface, kl, Tl[kl], kr, Tr[kr]))
    # Potential density difference, rho[kr] - rho[kl]
    dRho = 0.5 * ( ( dRdTr[kr] + dRdTl[kl] ) * ( Tr[kr] - Tl[kl] )
                 + ( dRdSr[kr] + dRdSl[kl] ) * ( Sr[kr] - Sl[kl] ) )
    # Which column has the lighter surface for the current indexes, krandkl
    if (not reached_bottom):
      if (dRho < 0.):
        searching_left_column = True
        searching_right_column = False
      elif (dRho > 0.):
        searching_right_column = True
        searching_left_column = False
      else: # dRho == 0.
        if (kl + kr == 0): # Still at surface
          searching_left_column = True
          searching_right_column = False
        else: #notthe surface so we simply change direction
          searching_left_column = not  searching_left_column
          searching_right_column = not  searching_right_column

    if (searching_left_column):
      # Interpolate for the neutral surface position within the left column, layer klm1
      # Potential density difference, rho(kl-1) - rho[kr] (should be negative)
      dRhoTop = 0.5 * ( ( dRdTl[klm1] + dRdTr[kr] ) * ( Tl[klm1] - Tr[kr] )
                      + ( dRdSl[klm1] + dRdSr[kr] ) * ( Sl[klm1] - Sr[kr] ) )
      # Potential density difference, rho[kl] - rho[kr] (will be positive)
      dRhoBot = 0.5 * ( ( dRdTl[klm1+1] + dRdTr[kr] ) * ( Tl[klm1+1] - Tr[kr] )
                      + ( dRdSl[klm1+1] + dRdSr[kr] ) * ( Sl[klm1+1] - Sr[kr] ) )
      print("Searching left interfaces (%d, %d): dRhoTop: %f dRhoBot: %f" % (klm1, klm1+1, dRhoTop, dRhoBot))

      # Because we are looking left, the right surface, kr, is lighter than klm1+1andshould be denser than klm1
      # unless we are still at the top of the left column (kl=1)
      if (dRhoTop > 0. or kr+kl==0):
        PoL[k_surface] = 0. # The right surface is lighter than anything in layer klm1
        print("At surface or left surface is lighter than layer krm1")
      elif (dRhoTop >= dRhoBot): # Left layer is unstratified
        PoL[k_surface] = 1.
        print("dRhoTop>=dRhoBot")
      else:
        # Linearly interpolate for the position between Pl(kl-1)andPl[kl] where the density difference
        # between rightandleft is zero.
        PoL[k_surface] = interpolate_for_nondim_position( dRhoTop, Pl[klm1], dRhoBot, Pl[klm1+1] )
        print("Interpolating for position")

      if (PoL[k_surface]>=1. and klm1<nk-1): # >= is really ==, when PoL==1 we point to the bottom of the cell
        klm1 = klm1 + 1
        PoL[k_surface] = PoL[k_surface] - 1.
        print("Point to bottom of cell")

      if ((klm1-lastK_left)+(PoL[k_surface]-lastP_left)<0.):
        PoL[k_surface] = lastP_left
        klm1 = lastK_left
        print("Rewind for some reason?")

      KoL[k_surface] = klm1
      if (kr <= nk - 1):
        PoR[k_surface] = 0.
        KoR[k_surface] = kr
      else:
        PoR[k_surface] = 1.
        KoR[k_surface] = nk-1

      if (kr <= nk-1):
        kr = kr + 1
      else:
        reached_bottom = True
        searching_right_column = True
        searching_left_column = False

    elif (searching_right_column):
      # Interpolate for the neutral surface position within the right column, layer krm1
      # Potential density difference, rho(kr-1) - rho[kl] (should be negative)
      dRhoTop = 0.5 * ( ( dRdTr[krm1] + dRdTl[kl] ) * ( Tr[krm1] - Tl[kl] )
                    + ( dRdSr[krm1] + dRdSl[kl] ) * ( Sr[krm1] - Sl[kl] ) )
      # Potential density difference, rho[kr] - rho[kl] (will be positive)
      dRhoBot = 0.5 * ( ( dRdTr[krm1+1] + dRdTl[kl] ) * ( Tr[krm1+1] - Tl[kl] )
                    + ( dRdSr[krm1+1] + dRdSl[kl] ) * ( Sr[krm1+1] - Sl[kl] ) )
      print("Searching right interfaces (%d, %d): dRhoTop: %f dRhoBot: %f" % (krm1, krm1+1, dRhoTop, dRhoBot))

      # Because we are looking right, the left surface, kl, is lighter than krm1+1andshould be denser than krm1
      # unless we are still at the top of the right column (kr=1)
      if (dRhoTop >= 0. or kr+kl==0):
        PoR[k_surface] = 0. # The left surface is lighter than anything in layer krm1
        print("At surface or left surface is lighter than layer krm1")
      elif (dRhoTop >= dRhoBot): # Right layer is unstratified
        PoR[k_surface] = 1.
        print("dRhoTop>=dRhoBot")
      else:
        # Linearly interpolate for the position between Pr(kr-1)andPr[kr] where the density difference
        # between rightandleft is zero.
        PoR[k_surface] = interpolate_for_nondim_position( dRhoTop, Pr[krm1], dRhoBot, Pr[krm1+1] )
        print("Interpolating for position")

      if (PoR[k_surface]>=1. and krm1<nk-1): # >= is really ==, when PoR==1 we point to the bottom of the cell
        krm1 = krm1 + 1
        PoR[k_surface] = PoR[k_surface] - 1.
        print("Point to bottom of cell")

      if ((krm1-lastK_right)+(PoR[k_surface]-lastP_right)<0.):
        PoR[k_surface] = lastP_right
        krm1 = lastK_right
        print("Rewind for some reason?")

      KoR[k_surface] = krm1
      if (kl <= nk -1):
        PoL[k_surface] = 0.
        KoL[k_surface] = kl
      else:
        PoL[k_surface] = 1.
        KoL[k_surface] = nk-1

      if (kl <= nk-1):
        kl = kl + 1
      else:
        reached_bottom = True
        searching_right_column = False
        searching_left_column = True

    else:
      print('Else what?')

    lastK_left = KoL[k_surface] ; lastP_left = PoL[k_surface]
    lastK_right = KoR[k_surface] ; lastP_right = PoR[k_surface]

    # Effective thickness
    #not: This would be better expressed in terms of the layers thicknesses rather
    # than as differences of position - AJA
    if (k_surface>0):
      PoL_abs[k_surface] = absolute_position(Pl, KoL, PoL, k_surface)
      PoR_abs[k_surface] = absolute_position(Pr, KoR, PoR, k_surface)
      hL[k_surface] = absolute_position(Pl,KoL,PoL,k_surface) - absolute_position(Pl,KoL,PoL,k_surface-1)
      hR[k_surface] = absolute_position(Pr,KoR,PoR,k_surface) - absolute_position(Pr,KoR,PoR,k_surface-1)
      if ( hL[k_surface] + hR[k_surface] > 0.):
        hEff[k_surface-1] = 2. * hL[k_surface] * hR[k_surface] / ( hL[k_surface] + hR[k_surface] ) # Analogous of effective resistance for two resistors
      else:
        hEff[k_surface-1] = 0.

  return PoL, PoR, PoL_abs, PoR_abs, KoL, KoR, hEff, hL, hR
