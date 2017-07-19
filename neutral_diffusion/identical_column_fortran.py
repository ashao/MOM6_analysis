#!/local/home/miniconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import ndiff
import pandas as pd
from imp import reload

ppm_routines = ndiff.ppm_routines
# Some simple functions to set up the and perform the neutral diffusion flux calculations
def ppm_reconstruction(h, Slayer):
    Sinterface = np.zeros(np.size(Slayer)+1, dtype=np.float32)
    ppm_routines.interface_scalar(h, Slayer, Sinterface, 2)
    Sl = np.zeros(np.size(Slayer), dtype=np.float32)
    Sr = np.zeros(np.size(Slayer), dtype=np.float32)
    ppm_routines.ppm_left_right_edge_values(Slayer.astype(np.float32), Sinterface, Sl, Sr)
    return Sl, Sr

def construct_column(h, T):
    zi = np.array(0) ; zi = np.append(zi,hl.cumsum())
    z_t = h.cumsum() - h[0]
    z_b = h.cumsum()
    z_c = h.cumsum() - 0.5*h
    T_t, T_b = ppm_reconstruction(h,T)

    return zi, z_t, z_b, z_c, T_t, T_b

def interp_reconstruction(T_t, T_c, T_b, P, k):
    if P == 0.:
        return T_t[k]
    if P == 1.:
        return T_b[k]
    elif P>0. and P<1.:
        T_int = np.interp( P, np.array( (0., 0.5, 1.) ), np.array( (T_t[k], T_c[k], T_b[k]) ) )
        return T_int
    else:
        print(k, P)

def plot_neutral_surfaces(zr_t, zr_c, zr_b, Tr_t, Tr, Tr_b, PoR_abs, KoR,
                          zl_t, zl_c, zl_b, Tl_t, Tl, Tl_b, PoL_abs, KoL, search_dir):
    # Plot PPM reconstructon of left column
    for k in np.arange(0,hl.size):
        z = np.array( [zl_t[k], zl_c[k], zl_b[k]])
        T = np.array( [Tl_t[k], Tl[k], Tl_b[k]] )
        plt.plot(T,z, color='green')
        plt.scatter(Tl_t, zl_t, marker='^')
        plt.scatter(Tl_b, zl_b, marker='v')
    # PPM reconstruction of right column
    for k in np.arange(0,hl.size):
        z = np.array( [zr_t[k], zr_c[k], zr_b[k]] )
        T = np.array( [Tr_t[k], Tr[k], Tr_b[k]] )
        plt.plot(T,z, color = 'cyan')
        plt.scatter(Tr_t, zr_t, marker='^')
        plt.scatter(Tr_b, zr_b, marker='v')

    # Plot neutral surfaces
    for k in np.arange(0,PoR_abs.size):
        kr = KoR[k]
        Tr_int = interp_reconstruction(Tr_t, Tr, Tr_b, PoR[k], kr)
        kl = KoL[k]
        Tl_int = interp_reconstruction(Tl_t, Tl, Tl_b, PoL[k], kl)
        if search_dir[k] == 0:
            #plt.plot(np.array( (Tl_int, Tr_int) ), np.array( (PoL_abs[k],PoR_abs[k]) ) , ':', color="green" )
            plt.arrow( Tl_int, PoL_abs[k], Tr_int - Tl_int , PoR_abs[k] - PoL_abs[k], head_width=0.5, color="green" )
        else:
            plt.arrow( Tr_int, PoR_abs[k], Tl_int - Tr_int , PoL_abs[k] - PoR_abs[k], head_width=0.5, color="cyan" )
            #plt.plot(np.array( (Tl_int, Tr_int) ), np.array( (PoL_abs[k],PoR_abs[k]) ) , ':', color="cyan" )
        print("Surface %d:" % k, Tl_int, Tr_int, PoL_abs[k], PoR_abs[k])
    plt.xlim( (5,22.5) )
    plt.grid(ls='dotted')
    plt.xlabel('Temperature')
    plt.gca().invert_yaxis()

# Set up some examples based on Alistair's schematics of
nk = 4
Sl_t = np.zeros(nk) ; Sl_b = np.zeros(nk)
Sr_t = np.zeros(nk) ; Sr_b = np.zeros(nk)
drdt_lt = -1*np.ones(nk) ; drdt_lb = -1*np.ones(nk)
drds_lt = np.zeros(nk) ; drds_lb = np.zeros(nk)
drdt_rt = -1*np.ones(nk) ; drdt_rb = -1*np.ones(nk)
drds_rt = np.zeros(nk) ; drds_rb = np.zeros(nk)

hl = np.array([10.,10.,10.,10.])
hr = np.array([10.,10.,10.,10.])
Tr = np.array([20.,16.,12.,10.])+2 ; Sr = np.zeros_like(Tr)
Tl = np.array([20.,16.,12.,10.]) ; Sl = np.zeros_like(Tl)
zil, zl_t, zl_b, zl_c, Tl_t, Tl_b = construct_column(hl, Tl)
zir, zr_t, zr_b, zr_c, Tr_t, Tr_b = construct_column(hr, Tr)
print(Tl_t,Tl_b)
print(Tr_t, Tr_b)


dRdTl = np.ones_like(Tl)*-1. ; dRdTr = np.ones_like(Tr)*-1.
dRdSl = np.zeros_like(Tl) ; dRdSr = np.zeros_like(Tr)
KoL = np.zeros(4*hl.size,dtype=np.int32) ; PoL = np.zeros(4*hl.size,dtype=np.float32)
KoR = np.zeros(4*hl.size,dtype=np.int32) ; PoR = np.zeros(4*hl.size,dtype=np.float32)
hEff = np.zeros(4*hl.size-1,dtype=np.float32)

Til = np.zeros((nk,2)) ; Tir = np.zeros((nk,2))
Sil = np.zeros((nk,2)) ; Sir = np.zeros((nk,2))
drdt_l = np.ones((nk,2))*-1 ; drdt_r = np.ones((nk,2))*-1
drds_l = np.zeros((nk,2)) ; drds_r = np.zeros((nk,2))

Til[:,0] = Tl_t ; Til[:,1] = Tl_b
Sil[:,0] = Sl_t ; Sil[:,1] = Tl_b
Tir[:,0] = Tr_t ; Tir[:,1] = Tr_b
Sir[:,0] = Sr_t ; Sir[:,1] = Tr_b


ndiff.ppm_routines.find_neutral_surface_positions_discontinuous(
    zil, Til, Sil, drdt_l, drds_l,
    zir, Tir, Sir, drdt_r, drds_r, PoL, PoR, KoL, KoR, hEff)
hEff = np.append(hEff,0.)
df_identical_fortran = pd.DataFrame({'PoL': PoL, 'PoR': PoR, 'KoL': (KoL), 'KoR': (KoR), 'hEff' : hEff})
interfaces = pd.DataFrame({'T TopLeft' : Tl_t, 'T BottomLeft' : Tl_b, 'T TopRight' : Tr_t, 'T BottomRight' : Tr_b})

print(df_identical_fortran)
print(interfaces)
