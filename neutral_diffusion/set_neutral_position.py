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

def absolute_position_discontinuous(Pint, Karr, NParr, k_surface):
    k = Karr[k_surface]
    return Pint[k] + NParr[k_surface] * ( Pint[k+1] - Pint[k] )

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

  kr = 0 ; lastK_right = 0 ; lastP_right = 0.
  kl = 0 ; lastK_left = 0 ; lastP_left = 0.
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
      same_k = np.floor( klm1*0.5 ) == np.floor( (klm1+1)*0.5 )
      same_p = Pl[klm1] == Pl[klm1+1]
      dRhoTop = 0.5 * ( ( dRdT_l[klm1] + dRdT_r[kr] ) * ( Tl[klm1] - Tr[kr] )
                      + ( dRdS_l[klm1] + dRdS_r[kr] ) * ( Sl[klm1] - Sr[kr] ) )
      dRhoBot = 0.5 * ( ( dRdT_l[klm1+1] + dRdT_r[kr] ) * ( Tl[klm1+1] - Tr[kr] )
                      + ( dRdS_l[klm1+1] + dRdS_r[kr] ) * ( Sl[klm1+1] - Sr[kr] ) )
      print("Searching left: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
      print("klm1: %d kl: %d kr: %d" % (klm1,kl,kr))
      # Search left
      if kr+kl==0:                # At surface
        PoL[k_surface] = 0.
      elif dRhoTop >= dRhoBot:
        PoL[k_surface] = 1.
      elif ( dRhoTop<0. and dRhoBot>0.) and ( not same_k ):
        PoL[k_surface] = 1.
      else:
        PoL[k_surface] = interpolate_for_nondim_position( dRhoTop, Pl[klm1], dRhoBot, Pl[klm1+1] )

      KoL[k_surface] = np.floor(0.5*klm1)
      KoR[k_surface] = np.floor(0.5*kr)
      if kr <= (2*nk-1):
        if np.floor(0.5*kr) == np.floor(0.5*(kr+1)):
          PoR[k_surface] = 0.
        else:
          PoR[k_surface] = 1.
      else:
        PoR[k_surface] = 1.
        KoR[k_surface] = nk-1
      if kr < 2*nk-1:
          kr = kr + 1
      else:
        reached_bottom = True
        searching_right_column = True
        searching_left_column = False
      lastK_left = klm1
      lastK_right = kr
    elif (searching_right_column):
      same_k = np.floor( krm1*0.5 ) == np.floor( (krm1+1)*0.5 )
      same_p = Pr[krm1] == Pr[krm1+1]
      dRhoTop = 0.5 * ( ( dRdT_r[krm1] + dRdT_l[kl] ) * ( Tr[krm1] - Tl[kl] )
                      + ( dRdS_r[krm1] + dRdS_l[kl] ) * ( Sr[krm1] - Sl[kl] ) )
      dRhoBot = 0.5 * ( ( dRdT_r[krm1+1] + dRdT_l[kl] ) * ( Tr[krm1+1] - Tl[kl] )
                      + ( dRdS_r[krm1+1] + dRdS_l[kl] ) * ( Sr[krm1+1] - Sl[kl] ) )
      print("Searching right: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
      print("krm1: %d kr: %d kl: %d" % (krm1,kr,kl))
      if kr+kl==0:                # At surface
        PoR[k_surface] = 0.       # Point to top of top layer
      elif dRhoTop >= dRhoBot:    # Top is unstratified or unstable
        PoR[k_surface] = 0.
      elif ( dRhoTop<0. and dRhoBot>0.) and ( not same_k ):
        PoR[k_surface] = 1.      # Point to bottom of layer of krm1
      else:
        PoR[k_surface] = interpolate_for_nondim_position( dRhoTop, Pr[krm1], dRhoBot, Pr[krm1+1] )

      KoR[k_surface] = np.floor(0.5*krm1)
      KoL[k_surface] = np.floor(0.5*kl)
      if kl <= (2*nk-1):
        if np.floor(0.5*kl) == np.floor(0.5*(kl+1)):
          PoL[k_surface] = 0.
        else:
          PoL[k_surface] = 1.
      else:
        PoL[k_surface] = 1.
        KoL[k_surface] = nk-1
      if kl < 2*nk-1:
        kl = kl + 1
      else:
        reached_bottom = True
        searching_right_column = False
        searching_left_column = True
      lastK_left = kl
      lastK_right = krm1
    else:
      print("ERROR")

    lastP_left = PoL[k_surface]
    lastP_right = PoR[k_surface]

    print("Position on left : %f" % PoL[k_surface])
    print("Position on right: %f" % PoR[k_surface])
    if k_surface>0:
      PoL_abs[k_surface] = absolute_position_discontinuous(Pres_l, KoL, PoL, k_surface)
      PoR_abs[k_surface] = absolute_position_discontinuous(Pres_r, KoR, PoR, k_surface)
      hL[k_surface] = absolute_position_discontinuous(Pres_l, KoL, PoL, k_surface) - absolute_position_discontinuous(Pres_l, KoL, PoL, k_surface-1)
      hR[k_surface] = absolute_position_discontinuous(Pres_r, KoR, PoR, k_surface) - absolute_position_discontinuous(Pres_r, KoR, PoR, k_surface-1)
      if (hL[k_surface] + hR[k_surface] > 0.):
          hEff[k_surface-1] = 2. * hL[k_surface] * hR[k_surface] / ( hL[k_surface]+hR[k_surface] )
      else:
          hEff[k_surface-1] = 0.
      print("hL: %f hR: %f hEff: %f" % (hL[k_surface],hR[k_surface],hEff[k_surface-1]))
  print(PoL,PoR)
  return PoL, PoR, PoL_abs, PoR_abs, KoL, KoR, hEff, hL, hR


def set_neutral_surface_position2( Pres_l, Tint_lt, Tint_lb, Sint_lt, Sint_lb, dRdT_lt, dRdT_lb, dRdS_lt, dRdS_lb, Pres_r, Tint_rt, Tint_rb, Sint_rt, Sint_rb, dRdT_rt, dRdT_rb, dRdS_rt, dRdS_rb ):
    nk = Tint_lt.size
    Pl = np.zeros(2*nk) ; Pr = np.zeros(2*nk)
    Tl = np.zeros(2*nk) ; Tr = np.zeros(2*nk)
    Sl = np.zeros(2*nk) ; Sr = np.zeros(2*nk)
    dRdT_l = np.zeros(2*nk) ; dRdT_r = np.zeros(2*nk) ;
    dRdS_l = np.zeros(2*nk) ; dRdS_r = np.zeros(2*nk) ;

    PoL = np.zeros(4*nk) ; PoR = np.zeros(4*nk)
    KoL = np.zeros(4*nk) ; KoR = np.zeros(4*nk)
    hL = np.zeros(4*nk) ; hR = np.zeros(4*nk)
    hEff = np.zeros(4*nk-1)

    for kl in np.arange(0,nk):
        Pl[2*kl-1] = Pres_l[kl]
        Pl[2*kl] = Pres_l[kl+1]
        Sl[2*kl-1] = Sint_lt[kl]
        Sl[2*kl] = Sint_lb[kl]
        Tl[2*kl-1] = Tint_lt[kl]
        Tl[2*kl] = Tint_lb[kl]
        dRdT_l[2*kl-1] = dRdT_lt[kl]
        dRdT_l[2*kl] = dRdT_lb[kl]
        dRdS_l[2*kl-1] = dRdS_lt[kl]
        dRdS_l[2*kl] = dRdS_lb[kl]
    for kr in np.arange(0,nk):
        Pr[2*kr-1] = Pres_r[kr]
        Pr[2*kr] = Pres_r[kr+1]
        Sr[2*kr-1] = Sint_rt[kr]
        Sr[2*kr] = Sint_rb[kr]
        Tr[2*kr-1] = Tint_rt[kr]
        Tr[2*kr] = Tint_rb[kr]
        dRdT_r[2*kr-1] = dRdT_rt[kr]
        dRdT_r[2*kr] = dRdT_rb[kr]
        dRdS_r[2*kr-1] = dRdS_rt[kr]
        dRdS_r[2*kr] = dRdS_rb[kr]
    kr = 0 ; lastK_right = 0 ; lastP_right = 0.
    kl = 0 ; lastK_left = 0 ; lastP_left = 0.
    reached_bottom = False

    for k_surface in np.arange(0,4*nk):
        print( "\nWorking on k_surface %d: Tl[%d]: %f Tr[%d]: %f" % (k_surface, kl, Tl[kl], kr, Tr[kr]))
        klm1 = max(kl-1,0)
        krm1 = max(kr-1,0)
        print("klm1: %d kl: %d krm1: %d kr: %d" % (klm1, kl, krm1, kr))
        dRho = 0.5 * ( ( dRdT_r[kr] + dRdT_l[kl] ) * ( Tr[kr] - Tl[kl] )
                      + (dRdS_r[kr] + dRdS_r[kr] ) * ( Sr[kr] - Sl[kl] ) )
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
            dRhoTop = 0.5 * ( ( dRdT_l[klm1] + dRdT_r[kr] ) * ( Tl[klm1] - Tr[kr] )
                            + ( dRdS_l[klm1] + dRdS_r[kr] ) * ( Sl[klm1] - Sr[kr] ) )
            dRhoBot = 0.5 * ( ( dRdT_l[klm1+1] + dRdT_r[kr] ) * ( Tl[klm1+1] - Tr[kr] )
                            + ( dRdS_l[klm1+1] + dRdS_r[kr] ) * ( Sl[klm1+1] - Sr[kr] ) )
            print("Searching left: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
            print("klm1: %d kr: %d" % (klm1,kr))
            # Handle how to handle where the position on the left is

            if dRhoTop > 0. or kr+kl==0:
              PoL[k_surface] = 0.
            elif dRhoTop >= dRhoBot:
              PoL[k_surface] = 1.
            elif ( dRhoTop<0. and dRhoBot>0.) and (Pl[klm1] == Pl[klm1+1]):
              PoL[k_surface] = 1.
            else:
              PoL[k_surface] = interpolate_for_nondim_position( dRhoTop, Pl[klm1], dRhoBot, Pl[klm1+1] )

            print("Position on left: %f" % PoL[k_surface])
            # When the position is 1, then we're at the bottom of a cell
            if PoL[k_surface]>=1 and klm1<(2*nk-1):
                klm1 += 1
                PoL[k_surface] = PoL[k_surface] - 1
                print("Point to bottom of cell")
            if ((np.floor(0.5*klm1) - np.floor(0.5*lastK_left)) + (PoL[k_surface]-lastP_left)) < 0.:
                PoL[k_surface] = lastP_left
                klm1 = lastK_left
                print("Point to bottom of previous cell")
            KoL[k_surface] = np.floor(0.5*klm1)
            if kr <= (2*nk-1):
                PoR[k_surface] = 0.
                KoR[k_surface] = np.floor( 0.5*kr )
            else:
                PoR[k_surface] = 1.
                KoR[k_surface] = nk-1
            if kr < 2*nk-1:
                kr = kr + 1
            else:
                reached_bottom = True
                searching_right_column = True
                searching_left_column = False
            lastK_left = klm1
            lastK_right = kr
        elif (searching_right_column):
            dRhoTop = 0.5 * ( ( dRdT_r[krm1] + dRdT_l[kl] ) * ( Tr[krm1] - Tl[kl] )
                            + ( dRdS_r[krm1] + dRdS_l[kl] ) * ( Sr[krm1] - Sl[kl] ) )
            dRhoBot = 0.5 * ( ( dRdT_r[krm1+1] + dRdT_l[kl] ) * ( Tr[krm1+1] - Tl[kl] )
                            + ( dRdS_r[krm1+1] + dRdS_l[kl] ) * ( Sr[krm1+1] - Sl[kl] ) )
            print("Searching right: dRhoTop: %f dRhoBot: %f" % (dRhoTop, dRhoBot))
            print("krm1: %d kl: %d" % (krm1,kl))
            if dRhoTop > 0. or kr+kl==0:
                PoR[k_surface] = 0.
            elif dRhoTop >= dRhoBot:
                PoR[k_surface] = 1.
            elif ( dRhoTop<0. and dRhoBot>0.) and ( Pr[krm1]==Pr[kr] ):
                PoR[k_surface] = 1.
            else:
                PoR[k_surface] = interpolate_for_nondim_position( dRhoTop, Pr[krm1], dRhoBot, Pr[krm1+1] )
            print("Position on right: %f" % PoR[k_surface])
            if PoR[k_surface]>=1 and krm1<(2*nk-1):
                krm1 += 1
#                PoR[k_surface] = PoR[k_surface] - 1
            if ((np.floor(0.5*krm1) - np.floor(0.5*lastK_right)) + (PoR[k_surface]-lastP_right)) < 0.:
                PoR[k_surface] = lastP_right
                krm1 = lastK_right
            KoR[k_surface] = np.floor(0.5*krm1)
            if kl <= (2*nk-1):
                PoL[k_surface] = 0.
                KoL[k_surface] = np.floor( 0.5*kl )
            else:
                PoL[k_surface] = 1.
                KoL[k_surface] = nk-1
            if kl < 2*nk-1:
                kl = kl + 1
            else:
                reached_bottom = True
                searching_right_column = False
                searching_left_column = True
            lastK_left = kl
            lastK_right = krm1
        else:
            print("ERROR")

        lastP_left = PoL[k_surface]
        lastP_right = PoR[k_surface]

        print("Position on left : %f" % PoL[k_surface])
        print("Position on right: %f" % PoR[k_surface])
        if k_surface>0:
            hL[k_surface] = absolute_position_discontinuous(Pl, KoL, PoL, k_surface) - absolute_position_discontinuous(Pl, KoL, PoL, k_surface-1)
            hR[k_surface] = absolute_position_discontinuous(Pr, KoR, PoR, k_surface) - absolute_position_discontinuous(Pr, KoR, PoR, k_surface-1)
            if (hL[k_surface] + hR[k_surface] > 0.):
                hEff[k_surface-1] = 2. * hL[k_surface] * hR[k_surface] / ( hL[k_surface]+hR[k_surface] )
            else:
                hEff[k_surface-1] = 0.
            print("hL: %f hR: %f hEff: %f" % (hL[k_surface],hR[k_surface],hEff[k_surface-1]))
    print(PoL,PoR)
    return PoL, PoR, KoL, KoR, hEff, hL, hR

