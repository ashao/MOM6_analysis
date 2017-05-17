!> A column-wise toolbox for implementing neutral diffusion
module PPM_routines

implicit none 

logical, parameter :: debug_this_module = .false.

contains

!> Returns interface scalar, Si, for a column of layer values, S.
subroutine interface_scalar(nk, h, S, Si, i_method)
  integer,               intent(in)    :: nk       !< Number of levels
  real, dimension(nk),   intent(in)    :: h        !< Layer thickness (H units)
  real, dimension(nk),   intent(in)    :: S        !< Layer scalar (conc, e.g. ppt)
  real, dimension(nk+1), intent(inout) :: Si       !< Interface scalar (conc, e.g. ppt)
  integer,               intent(in)    :: i_method !< =1 use average of PLM edges
                                                   !! =2 use continuous PPM edge interpolation
  ! Local variables
  integer :: k, km2, kp1
  real, dimension(nk) :: diff
  real :: Sb, Sa

  call PLM_diff(nk, h, S, 2, 1, diff)
  Si(1) = S(1) - 0.5 * diff(1)
  if (i_method==1) then
    do k = 2, nk
      ! Average of the two edge values (will be bounded and,
      ! when slopes are unlimited, notionally second-order accurate)
      Sa = S(k-1) + 0.5 * diff(k-1) ! Lower edge value of a PLM reconstruction for layer above
      Sb = S(k) - 0.5 * diff(k) ! Upper edge value of a PLM reconstruction for layer below
      Si(k) = 0.5 * ( Sa + Sb )
    enddo
  elseif (i_method==2) then
    do k = 2, nk
      ! PPM quasi-fourth order interpolation for edge values following
      ! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
      km2 = max(1, k-2)
      kp1 = min(nk, k+1)
      Si(k) = ppm_edge(h(km2), h(k-1), h(k), h(kp1),  S(k-1), S(k), diff(k-1), diff(k))
    enddo
  endif
  Si(nk+1) = S(nk) + 0.5 * diff(nk)

end subroutine interface_scalar

!> Returns the PPM quasi-fourth order edge value at k+1/2 following
!! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
real function ppm_edge(hkm1, hk, hkp1, hkp2,  Ak, Akp1, Pk, Pkp1)
  real, intent(in) :: hkm1 !< Width of cell k-1
  real, intent(in) :: hk   !< Width of cell k
  real, intent(in) :: hkp1 !< Width of cell k+1
  real, intent(in) :: hkp2 !< Width of cell k+2
  real, intent(in) :: Ak   !< Average scalar value of cell k
  real, intent(in) :: Akp1 !< Average scalar value of cell k+1
  real, intent(in) :: Pk   !< PLM slope for cell k
  real, intent(in) :: Pkp1 !< PLM slope for cell k+1

  ! Local variables
  real :: R_hk_hkp1, R_2hk_hkp1, R_hk_2hkp1, f1, f2, f3, f4
  real, parameter :: h_neglect = 1.e-30

  R_hk_hkp1 = hk + hkp1
  if (R_hk_hkp1 <= 0.) then
    ppm_edge = 0.5 * ( Ak + Akp1 )
    return
  endif
  R_hk_hkp1 = 1. / R_hk_hkp1
  if (hk<hkp1) then
    ppm_edge = Ak + ( hk * R_hk_hkp1 ) * ( Akp1 - Ak )
  else
    ppm_edge = Akp1 + ( hkp1 * R_hk_hkp1 ) * ( Ak - Akp1 )
  endif

  R_2hk_hkp1 = 1. / ( ( 2. * hk + hkp1 ) + h_neglect )
  R_hk_2hkp1 = 1. / ( ( hk + 2. * hkp1 ) + h_neglect )
  f1 = 1./ ( ( hk + hkp1) + ( hkm1 + hkp2 ) )
  f2 = 2. * ( hkp1 * hk ) * R_hk_hkp1 * &
            ( ( hkm1 + hk ) * R_2hk_hkp1  - ( hkp2 + hkp1 ) * R_hk_2hkp1 )
  f3 = hk * ( hkm1 + hk ) * R_2hk_hkp1
  f4 = hkp1 * ( hkp1 + hkp2 ) * R_hk_2hkp1

  ppm_edge = ppm_edge + f1 * ( f2 * ( Akp1 - Ak ) - ( f3 * Pkp1 - f4 * Pk ) )

end function ppm_edge

!> Returns the average of a PPM reconstruction between two
!! fractional positions.
real function ppm_ave(xL, xR, aL, aR, aMean)
  real, intent(in) :: xL    !< Fraction position of left bound (0,1)
  real, intent(in) :: xR    !< Fraction position of right bound (0,1)
  real, intent(in) :: aL    !< Left edge scalar value, at x=0
  real, intent(in) :: aR    !< Right edge scalar value, at x=1
  real, intent(in) :: aMean !< Average scalar value of cell

  ! Local variables
  real :: dx, xave, a6, a6o3

  dx = xR - xL
  xave = 0.5 * ( xR + xL )
  a6o3 = 2. * aMean - ( aL + aR ) ! a6 / 3.
  a6 = 3. * a6o3

  if (dx<0.) then
    stop 'ppm_ave: dx<0 should not happend!'
  elseif (dx>1.) then
    stop 'ppm_ave: dx>1 should not happend!'
  elseif (dx==0.) then
    ppm_ave = aL + ( aR - aL ) * xR + a6 * xR * ( 1. - xR )
  else
    ppm_ave = ( aL + xave * ( ( aR - aL ) + a6 ) )  - a6o3 * ( xR**2 + xR * xL + xL**2 )
  endif

end function ppm_ave

!> A true signum function that returns either -abs(a), when x<0; or abs(a) when x>0; or 0 when x=0.
real function signum(a,x)
  real, intent(in) :: a !< The magnitude argument
  real, intent(in) :: x !< The sign (or zero) argument

  signum = sign(a,x)
  if (x==0.) signum = 0.

end function signum

!> Returns PLM slopes for a column where the slopes are the difference in value across each cell.
!! The limiting follows equation 1.8 in Colella & Woodward, 1984: JCP 54, 174-201.
subroutine PLM_diff(nk, h, S, c_method, b_method, diff)
  integer,             intent(in)    :: nk       !< Number of levels
  real, dimension(nk), intent(in)    :: h        !< Layer thickness (H units)
  real, dimension(nk), intent(in)    :: S        !< Layer salinity (conc, e.g. ppt)
  integer,             intent(in)    :: c_method !< Method to use for the centered difference
  integer,             intent(in)    :: b_method !< =1, use PCM in first/last cell, =2 uses linear extrapolation
  real, dimension(nk), intent(inout) :: diff     !< Scalar difference across layer (conc, e.g. ppt)
                                                 !! determined by the following values for c_method:
                                                 !!   1. Second order finite difference (not recommended)
                                                 !!   2. Second order finite volume (used in original PPM)
                                                 !!   3. Finite-volume weighted least squares linear fit
                                                 !! \todo  The use of c_method to choose a scheme is inefficient
                                                 !! and should eventually be moved up the call tree.

  ! Local variables
  integer :: k
  real :: hkm1, hk, hkp1, Skm1, Sk, Skp1, diff_l, diff_r, diff_c

  do k = 2, nk-1
    hkm1 = h(k-1)
    hk = h(k)
    hkp1 = h(k+1)

    if ( ( hkp1 + hk ) * ( hkm1 + hk ) > 0.) then
      Skm1 = S(k-1)
      Sk = S(k)
      Skp1 = S(k+1)
      if (c_method==1) then
        ! Simple centered diff (from White)
        if ( hk + 0.5 * (hkm1 + hkp1) /= 0. ) then
          diff_c = ( Skp1 - Skm1 ) * ( hk / ( hk + 0.5 * (hkm1 + hkp1) ) )
        else
          diff_c = 0.
        endif
      elseif (c_method==2) then
        ! Second order accurate centered FV slope (from Colella and Woodward, JCP 1984)
        diff_c = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      elseif (c_method==3) then
        ! Second order accurate finite-volume least squares slope
        diff_c = hk * fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      endif
      ! Limit centered slope by twice the side differenced slopes
      diff_l = 2. * ( Sk - Skm1 )
      diff_r = 2. * ( Skp1 - Sk )
      if ( signum(1., diff_l) * signum(1., diff_r) <= 0. ) then
        diff(k) = 0. ! PCM for local extrema
      else
        diff(k) = sign( min( abs(diff_l), abs(diff_c), abs(diff_r) ), diff_c )
      endif
    else
      diff(k) = 0. ! PCM next to vanished layers
    endif
  enddo
  if (b_method==1) then ! PCM for top and bottom layer
    diff(1) = 0.
    diff(nk) = 0.
  elseif (b_method==2) then ! Linear extrapolation for top and bottom interfaces
    diff(1) = ( S(2) - S(1) ) * 2. * ( h(1) / ( h(1) + h(2) ) )
    diff(nk) = S(nk) - S(nk-1) * 2. * ( h(nk) / ( h(nk-1) + h(nk) ) )
  endif

end subroutine PLM_diff

!> Returns the cell-centered second-order finite volume (unlimited PLM) slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a difference across the central cell (i.e. units of scalar S).
!! Discretization follows equation 1.7 in Colella & Woodward, 1984: JCP 54, 174-201.
real function fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: h_sum, hp, hm

  h_sum = ( hkm1 + hkp1 ) + hk
  if (h_sum /= 0.) h_sum = 1./ h_sum
  hm =  hkm1 + hk
  if (hm /= 0.) hm = 1./ hm
  hp =  hkp1 + hk
  if (hp /= 0.) hp = 1./ hp
  fv_diff = ( hk * h_sum ) * &
            (   ( 2. * hkm1 + hk ) * hp * ( Skp1 - Sk ) &
              + ( 2. * hkp1 + hk ) * hm * ( Sk - Skm1 ) )
end function fv_diff


!> Returns the cell-centered second-order weighted least squares slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a gradient (i.e. units of scalar S over width units).
real function fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: xkm1, xkp1
  real :: h_sum, hx_sum, hxsq_sum, hxy_sum, hy_sum, det

  xkm1 = -0.5 * ( hk + hkm1 )
  xkp1 = 0.5 * ( hk + hkp1 )
  h_sum = ( hkm1 + hkp1 ) + hk
  hx_sum = hkm1*xkm1 + hkp1*xkp1
  hxsq_sum = hkm1*(xkm1**2) + hkp1*(xkp1**2)
  hxy_sum = hkm1*xkm1*Skm1 + hkp1*xkp1*Skp1
  hy_sum = ( hkm1*Skm1 + hkp1*Skp1 ) + hk*Sk
  det = h_sum * hxsq_sum - hx_sum**2
  if (det /= 0.) then
    !a = ( hxsq_sum * hy_sum - hx_sum*hxy_sum ) / det ! a would be mean of straight line fit
    fvlsq_slope = ( h_sum * hxy_sum - hx_sum*hy_sum ) / det ! Gradient of straight line fit
  else
    fvlsq_slope = 0. ! Adcroft's reciprocal rule
  endif
end function fvlsq_slope


!> Returns positions within left/right columns of combined interfaces
subroutine find_neutral_surface_positions(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Fractional position of neutral surface within layer KoR of right column
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces (Pa)

  ! Local variables
  integer :: k_surface              ! Index of neutral surface
  integer :: kl                     ! Index of left interface
  integer :: kr                     ! Index of right interface
  real    :: dRdT, dRdS             ! dRho/dT and dRho/dS for the neutral surface
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  integer :: krm1, klm1
  real    :: dRho, dRhoTop, dRhoBot, hL, hR
  integer :: lastK_left, lastK_right
  real    :: lastP_left, lastP_right

  ! Initialize variables for the search
  kr = 1 ; lastK_right = 1 ; lastP_right = 0.
  kl = 1 ; lastK_left = 1 ; lastP_left = 0.
  reached_bottom = .false.

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, 2*nk+2
    klm1 = max(kl-1, 1)
    if (klm1>nk) stop 'find_neutral_surface_positions(): klm1 went out of bounds!'
    krm1 = max(kr-1, 1)
    if (krm1>nk) stop 'find_neutral_surface_positions(): krm1 went out of bounds!'

    ! Potential density difference, rho(kr) - rho(kl)
    dRho = 0.5 * ( ( dRdTr(kr) + dRdTl(kl) ) * ( Tr(kr) - Tl(kl) ) &
                 + ( dRdSr(kr) + dRdSl(kl) ) * ( Sr(kr) - Sl(kl) ) )
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (.not. reached_bottom) then
      if (dRho < 0.) then
        searching_left_column = .true.
        searching_right_column = .false.
      elseif (dRho > 0.) then
        searching_right_column = .true.
        searching_left_column = .false.
      else ! dRho == 0.
        if (kl + kr == 2) then ! Still at surface
          searching_left_column = .true.
          searching_right_column = .false.
        else ! Not the surface so we simply change direction
          searching_left_column = .not.  searching_left_column
          searching_right_column = .not.  searching_right_column
        endif
      endif
    endif

    if (searching_left_column) then
      ! Interpolate for the neutral surface position within the left column, layer klm1
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTl(klm1) + dRdTr(kr) ) * ( Tl(klm1) - Tr(kr) ) &
                     + ( dRdSl(klm1) + dRdSr(kr) ) * ( Sl(klm1) - Sr(kr) ) )
      ! Potential density difference, rho(kl) - rho(kr) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTl(klm1+1) + dRdTr(kr) ) * ( Tl(klm1+1) - Tr(kr) ) &
                   + ( dRdSl(klm1+1) + dRdSr(kr) ) * ( Sl(klm1+1) - Sr(kr) ) )

      ! Because we are looking left, the right surface, kr, is lighter than klm1+1 and should be denser than klm1
      ! unless we are still at the top of the left column (kl=1)
      if (dRhoTop > 0. .or. kr+kl==2) then
        PoL(k_surface) = 0. ! The right surface is lighter than anything in layer klm1
      elseif (dRhoTop >= dRhoBot) then ! Left layer is unstratified
        PoL(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pl(kl-1) and Pl(kl) where the density difference
        ! between right and left is zero.
        PoL(k_surface) = interpolate_for_nondim_position( dRhoTop, Pl(klm1), dRhoBot, Pl(klm1+1) )
      endif
      if (PoL(k_surface)>=1. .and. klm1<nk) then ! >= is really ==, when PoL==1 we point to the bottom of the cell
        klm1 = klm1 + 1
        PoL(k_surface) = PoL(k_surface) - 1.
      endif
      if (real(klm1-lastK_left)+(PoL(k_surface)-lastP_left)<0.) then
        PoL(k_surface) = lastP_left
        klm1 = lastK_left
      endif
      KoL(k_surface) = klm1
      if (kr <= nk) then
        PoR(k_surface) = 0.
        KoR(k_surface) = kr
      else
        PoR(k_surface) = 1.
        KoR(k_surface) = nk
      endif
      if (kr <= nk) then
        kr = kr + 1
      else
        reached_bottom = .true.
        searching_right_column = .true.
        searching_left_column = .false.
      endif
    elseif (searching_right_column) then
      ! Interpolate for the neutral surface position within the right column, layer krm1
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTr(krm1) + dRdTl(kl) ) * ( Tr(krm1) - Tl(kl) ) &
                     + ( dRdSr(krm1) + dRdSl(kl) ) * ( Sr(krm1) - Sl(kl) ) )
      ! Potential density difference, rho(kr) - rho(kl) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTr(krm1+1) + dRdTl(kl) ) * ( Tr(krm1+1) - Tl(kl) ) &
                   + ( dRdSr(krm1+1) + dRdSl(kl) ) * ( Sr(krm1+1) - Sl(kl) ) )

      ! Because we are looking right, the left surface, kl, is lighter than krm1+1 and should be denser than krm1
      ! unless we are still at the top of the right column (kr=1)
      if (dRhoTop >= 0. .or. kr+kl==2) then
        PoR(k_surface) = 0. ! The left surface is lighter than anything in layer krm1
      elseif (dRhoTop >= dRhoBot) then ! Right layer is unstratified
        PoR(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference
        ! between right and left is zero.
        PoR(k_surface) = interpolate_for_nondim_position( dRhoTop, Pr(krm1), dRhoBot, Pr(krm1+1) )
      endif
      if (PoR(k_surface)>=1. .and. krm1<nk) then ! >= is really ==, when PoR==1 we point to the bottom of the cell
        krm1 = krm1 + 1
        PoR(k_surface) = PoR(k_surface) - 1.
      endif
      if (real(krm1-lastK_right)+(PoR(k_surface)-lastP_right)<0.) then
        PoR(k_surface) = lastP_right
        krm1 = lastK_right
      endif
      KoR(k_surface) = krm1
      if (kl <= nk) then
        PoL(k_surface) = 0.
        KoL(k_surface) = kl
      else
        PoL(k_surface) = 1.
        KoL(k_surface) = nk
      endif
      if (kl <= nk) then
        kl = kl + 1
      else
        reached_bottom = .true.
        searching_right_column = .false.
        searching_left_column = .true.
      endif
    else
      stop 'Else what?'
    endif

    lastK_left = KoL(k_surface) ; lastP_left = PoL(k_surface)
    lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)

    ! Effective thickness
    ! NOTE: This would be better expressed in terms of the layers thicknesses rather
    ! than as differences of position - AJA
    if (k_surface>1) then
      hL = absolute_position(nk,Pl,KoL,PoL,k_surface) - absolute_position(nk,Pl,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,Pr,KoR,PoR,k_surface) - absolute_position(nk,Pr,KoR,PoR,k_surface-1)
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR )
      else
        hEff(k_surface-1) = 0.
      endif
    endif

  enddo neutral_surfaces

end subroutine find_neutral_surface_positions

!> Converts non-dimensional position within a layer to absolute position (for debugging)
real function absolute_position(n,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces (Pa)
  integer, intent(in) :: Karr(2*n+2)  !< Index of interface above position
  real,    intent(in) :: NParr(2*n+2) !< Non-dimensional position within layer Karr(:)

  ! Local variables
  integer :: k_surface, k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Converts non-dimensional positions within layers to absolute positions (for debugging)
function absolute_positions(n,Pint,Karr,NParr)
  integer, intent(in) :: n            !< Number of levels
  real,    intent(in) :: Pint(n+1)    !< Position of interface (Pa)
  integer, intent(in) :: Karr(2*n+2)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(2*n+2) !< Non-dimensional positions within layers Karr(:)

  real,  dimension(2*n+2) :: absolute_positions ! Absolute positions (Pa)

  ! Local variables
  integer :: k_surface, k

  do k_surface = 1, 2*n+2
    absolute_positions(k_surface) = absolute_position(n,Pint,Karr,NParr,k_surface)
  enddo

end function absolute_positions

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) stop 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg'
  if (dRhoNeg>dRhoPos) write(0,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
  if (dRhoNeg>dRhoPos) stop 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos'
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRhoPos - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

subroutine ppm_left_right_edge_values(nk, Tl, Ti, aL, aR)
  integer,                    intent(in)    :: nk !< Number of levels
  real, dimension(nk),        intent(in)    :: Tl !< Layer tracer (conc, e.g. degC)
  real, dimension(nk+1),      intent(in)    :: Ti !< Interface tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aL !< Left edge value of tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aR !< Right edge value of tracer (conc, e.g. degC)

  integer :: k 
  ! Setup reconstruction edge values
  do k = 1, nk
    aL(k) = Ti(k)
    aR(k) = Ti(k+1)
    if ( signum(1., aR(k) - Tl(k))*signum(1., Tl(k) - aL(k)) <= 0.0 ) then
      aL(k) = Tl(k)
      aR(k) = Tl(k)
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) > abs(aR(k) - aL(k)) ) then
      aL(k) = Tl(k) + 2.0 * ( Tl(k) - aR(k) )
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) < -abs(aR(k) - aL(k)) ) then
      aR(k) = Tl(k) + 2.0 * ( Tl(k) - aL(k) )
    endif
  enddo
end subroutine ppm_left_right_edge_values

!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, hEff, Flx)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk),        intent(in)    :: hl    !< Left-column layer thickness (Pa)
  real, dimension(nk),        intent(in)    :: hr    !< Right-column layer thickness (Pa)
  real, dimension(nk),        intent(in)    :: Tl    !< Left-column layer tracer (conc, e.g. degC)
  real, dimension(nk),        intent(in)    :: Tr    !< Right-column layer tracer (conc, e.g. degC)
  real, dimension(2*nk+2),    intent(in)    :: PiL   !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(2*nk+2),    intent(in)    :: PiR   !< Fractional position of neutral surface within layer KoR of right column
  integer, dimension(2*nk+2), intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(in)    :: hEff  !< Effective thickness between two neutral surfaces (Pa)
  real, dimension(2*nk+1),    intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)

  ! Local variables
  integer :: k_sublayer, klb, klt, krb, krt, k
  real :: T_right_top, T_right_bottom, T_right_layer
  real :: T_left_top, T_left_bottom, T_left_layer
  real :: dT_top, dT_bottom, dT_layer, dT_ave
  real, dimension(nk+1) :: Til !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk+1) :: Tir !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_l !< Left-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_l !< Left-column right edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_r !< Right-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_r !< Right-column right edge value of tracer (conc, e.g. degC)

  call interface_scalar(nk, hl, Tl, Til, 2)
  call interface_scalar(nk, hr, Tr, Tir, 2)

  ! Setup reconstruction edge values
  do k = 1, nk
    aL_l(k) = Til(k)
    aR_l(k) = Til(k+1)
    if ( signum(1., aR_l(k) - Tl(k))*signum(1., Tl(k) - aL_l(k)) <= 0.0 ) then
      aL_l(k) = Tl(k)
      aR_l(k) = Tl(k)
    elseif ( sign(3., aR_l(k) - aL_l(k)) * ( (Tl(k) - aL_l(k)) + (Tl(k) - aR_l(k))) > abs(aR_l(k) - aL_l(k)) ) then
      aL_l(k) = Tl(k) + 2.0 * ( Tl(k) - aR_l(k) )
    elseif ( sign(3., aR_l(k) - aL_l(k)) * ( (Tl(k) - aL_l(k)) + (Tl(k) - aR_l(k))) < -abs(aR_l(k) - aL_l(k)) ) then
      aR_l(k) = Tl(k) + 2.0 * ( Tl(k) - aL_l(k) )
    endif
    aL_r(k) = Tir(k)
    aR_r(k) = Tir(k+1)
    if ( signum(1., aR_r(k) - Tr(k))*signum(1., Tr(k) - aL_r(k)) <= 0.0 ) then
      aL_r(k) = Tr(k)
      aR_r(k) = Tr(k)
    elseif ( sign(3., aR_r(k) - aL_r(k)) * ( (Tr(k) - aL_r(k)) + (Tr(k) - aR_r(k))) > abs(aR_r(k) - aL_r(k)) ) then
      aL_r(k) = Tr(k) + 2.0 * ( Tr(k) - aR_r(k) )
    elseif ( sign(3., aR_r(k) - aL_r(k)) * ( (Tr(k) - aL_r(k)) + (Tr(k) - aR_r(k))) < -abs(aR_r(k) - aL_r(k)) ) then
      aR_r(k) = Tr(k) + 2.0 * ( Tr(k) - aL_r(k) )
    endif
  enddo

  do k_sublayer = 1, 2*nk+1
    if (hEff(k_sublayer) == 0.) then
      Flx(k_sublayer) = 0.
    else

      klb = KoL(k_sublayer+1)
      T_left_bottom = ( 1. - PiL(k_sublayer+1) ) * Til(klb) + PiL(k_sublayer+1) * Til(klb+1)

      klt = KoL(k_sublayer)
      T_left_top = ( 1. - PiL(k_sublayer) ) * Til(klt) + PiL(k_sublayer) * Til(klt+1)

      T_left_layer = ppm_ave(PiL(k_sublayer), PiL(k_sublayer+1) + real(klb-klt), &
                             aL_l(klt), aR_l(klt), Tl(klt))

      krb = KoR(k_sublayer+1)
      T_right_bottom = ( 1. - PiR(k_sublayer+1) ) * Tir(krb) + PiR(k_sublayer+1) * Tir(krb+1)

      krt = KoR(k_sublayer)
      T_right_top = ( 1. - PiR(k_sublayer) ) * Tir(krt) + PiR(k_sublayer) * Tir(krt+1)

      T_right_layer = ppm_ave(PiR(k_sublayer), PiR(k_sublayer+1) + real(krb-krt), &
                              aL_r(krt), aR_r(krt), Tr(krt))

      dT_top = T_right_top - T_left_top
      dT_bottom = T_right_bottom - T_left_bottom
      dT_ave = 0.5 * ( dT_top + dT_bottom )
      dT_layer = T_right_layer - T_left_layer
      if (signum(1.,dT_top) * signum(1.,dT_bottom) <= 0. .or. signum(1.,dT_ave) * signum(1.,dT_layer) <= 0. ) then
        dT_ave = 0.
      else
        dT_ave = dT_layer
      endif
      Flx(k_sublayer) = dT_ave * hEff(k_sublayer)
    endif
  enddo

end subroutine neutral_surface_flux

!> Compares a single row, k, of output from find_neutral_surface_positions()
logical function compare_nsp_row(KoL, KoR, pL, pR, KoL0, KoR0, pL0, pR0)
  integer,  intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer,  intent(in) :: KoR   !< Index of first right interface above neutral surface
  real,     intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real,     intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  integer,  intent(in) :: KoL0  !< Correct value for KoL
  integer,  intent(in) :: KoR0  !< Correct value for KoR
  real,     intent(in) :: pL0   !< Correct value for pL
  real,     intent(in) :: pR0   !< Correct value for pR

  compare_nsp_row = .false.
  if (KoL /= KoL0) compare_nsp_row = .true.
  if (KoR /= KoR0) compare_nsp_row = .true.
  if (pL /= pL0) compare_nsp_row = .true.
  if (pR /= pR0) compare_nsp_row = .true.
end function compare_nsp_row

end module PPM_routines
