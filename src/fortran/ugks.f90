!-------------------------------------------------------------------
! Unified Gas Kinetic Scheme
!-------------------------------------------------------------------

subroutine flux_ugks_1f1v(fluxw, fluxh, wL, hL, wR, hR, unum, uspace, weight, &
                          ink, gamma, muref, omega, prandtl, dt, lenL, lenR, shL, shR)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum) !// interface fluxes
    real(kind=8), intent(in) :: wL(3), wR(3) !// conservative variables
    real(kind=8), intent(in) :: hL(unum), shL(unum) !// distribution functions and their slopes in left cell
    real(kind=8), intent(in) :: hR(unum), shR(unum) !// distribution functions and their slopes in right cell
    real(kind=8), intent(in) :: lenL, lenR !// cell lengths
    real(kind=8), intent(in) :: uspace(unum), weight(unum) !// velocity quadrature points and weights
    real(kind=8), intent(in) :: ink, gamma !// internal degrees of freedom of gas and Poisson ratio
    real(kind=8), intent(in) :: muref, omega, prandtl !// reference viscosity, VHS model index, and Prandtl number
    real(kind=8), intent(in) :: dt !// time step

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    ! interface variable
    real(kind=8) :: h(unum)
    real(kind=8) :: H0(unum)
    real(kind=8) :: H_plus(unum)
    real(kind=8) :: sh(unum)
    real(kind=8) :: w(3), prim(3)
    real(kind=8) :: qf
    real(kind=8) :: sw(3)
    real(kind=8) :: aL(3), aR(3), aT(3)

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mxi(0:2)
    real(kind=8) :: Mau_0(3), Mau_L(3), Mau_R(3), Mau_T(3)
    real(kind=8) :: tau
    real(kind=8) :: Mt(5)

    ! heaviside step function
    real(kind=8) :: delta(unum)
    delta = (sign(1.d0, uspace) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! upwind reconstruction
    !--------------------------------------------------
    h = (hL + 0.5 * lenL * shL) * delta + &
        (hR - 0.5 * lenR * shR) * (1.d0 - delta)

    sh = shL * delta + shR * (1.d0 - delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    ! conservative variables W_0 
    w(1) = sum(weight * h)
    w(2) = sum(weight * uspace * h)
    w(3) = 0.5d0 * sum(weight * uspace**2 * h)

    ! convert to primary variables
    prim = conserve_prim_1d(w, gamma)

    ! heat flux
    qf = heat_flux_1f1v(h, prim, uspace, weight)

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw = (w - wL) / (0.5 * lenL)
    aL = micro_slope_1d(prim, sw, ink)

    sw = (wR - w) / (0.5 * lenR)
    aR = micro_slope_1d(prim, sw, ink)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    ! <u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment_1d(prim, Mu, Mxi, Mu_L, Mu_R, ink) 

    Mau_L = moment_au_1d(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
    Mau_R = moment_au_1d(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

    sw = -prim(1) * (Mau_L + Mau_R) !time slope of W
    aT = micro_slope_1d(prim, sw, ink) !calculate A

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_1d(prim, muref, omega)

    Mt(4) = tau * (1.d0 - exp(-dt / tau))
    Mt(5) = -tau * dt * exp(-dt / tau) + tau * Mt(4)
    Mt(1) = dt - Mt(4)
    Mt(2) = -tau * Mt(1) + Mt(5) 
    Mt(3) = dt**2 / 2.d0 - tau * Mt(1)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Mau_0 = moment_uv_1d(Mu, Mxi, 1, 0) !<u*\psi>
    Mau_L = moment_au_1d(aL, Mu_L, Mxi, 2) !<aL*u^2*\psi>_{>0}
    Mau_R = moment_au_1d(aR, Mu_R, Mxi, 2) !<aR*u^2*\psi>_{<0}
    Mau_T = moment_au_1d(aT, Mu, Mxi, 1) !<A*u*\psi>

    fluxw = Mt(1) * prim(1) * Mau_0 + Mt(2) * prim(1) * (Mau_L + Mau_R) + Mt(3) * prim(1) * Mau_T

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_1f1v(H0, prim, uspace)

    ! Shakhov part H+ and B+
    call shakhov_1f1v(H0, qf, prim, H_plus, uspace, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * uspace * H_plus) + &
               Mt(4) * sum(weight * uspace * h) - Mt(5) * sum(weight * uspace**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * uspace**2 * H_plus) + &
               Mt(4) * sum(weight * uspace**2 * h) - Mt(5) * sum(weight * uspace**3 * sh)
    fluxw(3) = fluxw(3) + &
               Mt(1) * 0.5d0 * sum(weight * uspace * uspace**2 * H_plus) + &
               Mt(4) * 0.5d0 * sum(weight * uspace * uspace**2 * h) - &
               Mt(5) * 0.5d0 * sum(weight * uspace**2 * uspace**2 * sh)

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * uspace * (H0 + H_plus) + &
            Mt(2) * uspace**2*(aL(1) * H0 + aL(2) * uspace * H0 + 0.5d0 * aL(3) * uspace**2 * H0) * delta + &
            Mt(2) * uspace**2*(aR(1) * H0 + aR(2) * uspace * H0 + 0.5d0 * aR(3) * uspace**2 * H0) * (1.d0 - delta) + &
            Mt(3) * uspace *(aT(1) * H0 + aT(2) * uspace * H0 + 0.5d0 * aT(3) * uspace**2 * H0) + &
            Mt(4) * uspace * h - Mt(5) * uspace**2 * sh

end

!-------------------------------------------------------------------

subroutine flux_ugks_2f1v(fluxw, fluxh, fluxb, wL, hL, bL, wR, hR, bR, unum, uspace, weight, &
                          ink, gamma, muref, omega, prandtl, dt, lenL, lenR, shL, sbL, shR, sbR)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum), fluxb(unum) !// interface fluxes
    real(kind=8), intent(in) :: wL(3), wR(3) !// conservative variables
    real(kind=8), intent(in) :: hL(unum), bL(unum), shL(unum), sbL(unum) !// distribution functions and their slopes in left cell
    real(kind=8), intent(in) :: hR(unum), bR(unum), shR(unum), sbR(unum) !// distribution functions and their slopes in right cell
    real(kind=8), intent(in) :: lenL, lenR !// cell lengths
    real(kind=8), intent(in) :: uspace(unum), weight(unum) !// velocity quadrature points and weights
    real(kind=8), intent(in) :: ink, gamma !// internal degrees of freedom of gas and Poisson ratio
    real(kind=8), intent(in) :: muref, omega, prandtl !// reference viscosity, VHS model index, and Prandtl number
    real(kind=8), intent(in) :: dt !// time step

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    ! interface variable
    real(kind=8) :: h(unum), b(unum)
    real(kind=8) :: H0(unum), B0(unum)
    real(kind=8) :: H_plus(unum), B_plus(unum)
    real(kind=8) :: sh(unum), sb(unum)
    real(kind=8) :: w(3), prim(3)
    real(kind=8) :: qf
    real(kind=8) :: sw(3)
    real(kind=8) :: aL(3), aR(3), aT(3)

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mxi(0:2)
    real(kind=8) :: Mau_0(3), Mau_L(3), Mau_R(3), Mau_T(3)
    real(kind=8) :: tau
    real(kind=8) :: Mt(5)

    ! heaviside step function
    real(kind=8) :: delta(unum)
    delta = (sign(1.d0, uspace) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! upwind reconstruction
    !--------------------------------------------------
    h = hL * delta + hR * (1.d0 - delta)
    b = bL * delta + bR * (1.d0 - delta)

    sh = shL * delta + shR * (1.d0 - delta)
    sb = sbL * delta + sbR * (1.d0 - delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    ! conservative variables W_0 
    w(1) = sum(weight * h)
    w(2) = sum(weight * uspace * h)
    w(3) = 0.5d0 * (sum(weight * uspace**2 * h) + sum(weight * b))

    ! convert to primary variables
    prim = conserve_prim_1d(w, gamma)

    ! heat flux
    qf = heat_flux_2f1v(h, b, prim, uspace, weight)

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw = (w - wL) / (0.5 * lenL)
    aL = micro_slope_1d(prim, sw, ink)

    sw = (wR - w) / (0.5 * lenR)
    aR = micro_slope_1d(prim, sw, ink)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    ! <u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment_1d(prim, Mu, Mxi, Mu_L, Mu_R, ink) 

    Mau_L = moment_au_1d(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
    Mau_R = moment_au_1d(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

    sw = -prim(1) * (Mau_L + Mau_R) !time slope of W
    aT = micro_slope_1d(prim, sw, ink) !calculate A

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_1d(prim, muref, omega)

    Mt(4) = tau * (1.d0 - exp(-dt / tau))
    Mt(5) = -tau * dt * exp(-dt / tau) + tau * Mt(4)
    Mt(1) = dt - Mt(4)
    Mt(2) = -tau * Mt(1) + Mt(5) 
    Mt(3) = dt**2 / 2.d0 - tau * Mt(1)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Mau_0 = moment_uv_1d(Mu, Mxi, 1, 0) !<u*\psi>
    Mau_L = moment_au_1d(aL, Mu_L, Mxi, 2) !<aL*u^2*\psi>_{>0}
    Mau_R = moment_au_1d(aR, Mu_R, Mxi, 2) !<aR*u^2*\psi>_{<0}
    Mau_T = moment_au_1d(aT, Mu, Mxi, 1) !<A*u*\psi>

    fluxw = Mt(1) * prim(1) * Mau_0 + Mt(2) * prim(1) * (Mau_L + Mau_R) + Mt(3) * prim(1) * Mau_T

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_2f1v(H0, B0, prim, uspace, ink)

    ! Shakhov part H+ and B+
    call shakhov_2f1v(H0, B0, qf, prim, H_plus, B_plus, uspace, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * uspace * H_plus) + &
               Mt(4) * sum(weight * uspace * h) - Mt(5) * sum(weight * uspace**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * uspace**2 * H_plus) + &
               Mt(4) * sum(weight * uspace**2 * h) - Mt(5) * sum(weight * uspace**3 * sh)
    fluxw(3) = fluxw(3) + &
               Mt(1) * 0.5d0 * (sum(weight * uspace * uspace**2 * H_plus) + sum(weight * uspace * B_plus)) + &
               Mt(4) * 0.5d0 * (sum(weight * uspace * uspace**2 * h) + sum(weight * uspace * b)) - &
               Mt(5) * 0.5d0 * (sum(weight * uspace**2 * uspace**2 * sh) + sum(weight * uspace**2 * sb))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * uspace * (H0 + H_plus) + &
            Mt(2) * uspace**2*(aL(1) * H0 + aL(2) * uspace * H0 + 0.5d0 * aL(3) * (uspace**2 * H0 + B0)) * delta + &
            Mt(2) * uspace**2*(aR(1) * H0 + aR(2) * uspace * H0 + 0.5d0 * aR(3) * (uspace**2 * H0 + B0)) * (1.d0 - delta) + &
            Mt(3) * uspace *(aT(1) * H0 + aT(2) * uspace * H0 + 0.5d0 * aT(3) * (uspace**2 * H0 + B0)) + &
            Mt(4) * uspace * h - Mt(5) * uspace**2 * sh

    fluxb = Mt(1) * uspace * (B0 + B_plus) + &
            Mt(2) * uspace**2*(aL(1) * B0 + aL(2) * uspace * B0 + 0.5 * aL(3) * (uspace**2 * B0 + Mxi(2) * H0)) * delta + &
            Mt(2) * uspace**2*(aR(1) * B0 + aR(2) * uspace * B0 + 0.5 * aR(3) * (uspace**2 * B0 + Mxi(2) * H0)) * (1.d0 - delta) + &
            Mt(3) * uspace*(aT(1) * B0 + aT(2) * uspace * B0 + 0.5 * aT(3) * (uspace**2 * B0 + Mxi(2) * H0)) + &
            Mt(4) * uspace * b - Mt(5) * uspace**2 * sb

end

!-------------------------------------------------------------------

subroutine flux_ugks_2f2v(fluxw, fluxh, fluxb, &
                          wL, hL, bL, &
                          wR, hR, bR, &
                          unum, vnum, vn, vt, weight, &
                          ink, gamma, muref, omega, prandtl, &
                          dt, lenL, lenR, lenFace, shL, sbL, shR, sbR)

    integer, intent(in) :: unum, vnum
    real(kind=8), intent(inout) :: fluxw(4), fluxh(unum, vnum), fluxb(unum, vnum)
    real(kind=8), intent(in) :: wL(4), wR(4)
    real(kind=8), intent(in) :: hL(unum, vnum), bL(unum, vnum), shL(unum, vnum), sbL(unum, vnum)
    real(kind=8), intent(in) :: hR(unum, vnum), bR(unum, vnum), shR(unum, vnum), sbR(unum, vnum)
    real(kind=8), intent(in) :: lenL, lenR
    real(kind=8), intent(in) :: vn(unum, vnum), vt(unum, vnum), weight(unum, vnum)
    real(kind=8), intent(in) :: ink, gamma
    real(kind=8), intent(in) :: muref, omega, prandtl
    real(kind=8), intent(in) :: dt, lenFace

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    real(kind=8) :: delta(unum, vnum)

    ! interface variable
    real(kind=8) :: h(unum, vnum), b(unum, vnum)
    real(kind=8) :: H0(unum, vnum), B0(unum, vnum)
    real(kind=8) :: H_plus(unum, vnum), B_plus(unum, vnum)
    real(kind=8) :: sh(unum, vnum), sb(unum, vnum)
    real(kind=8) :: w(4), prim(4)
    real(kind=8) :: qf(2)
    real(kind=8) :: sw(4)
    real(kind=8) :: aL(4), aR(4), aT(4)

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    real(kind=8) :: Muv(4)
    real(kind=8) :: Mau_L(4), Mau_R(4)
    real(kind=8) :: Mau_T(4)
    real(kind=8) :: tau
    real(kind=8) :: Mt(5)

    delta = (sign(1.d0, vn) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! reconstruct initial distribution
    !--------------------------------------------------
    h = hL * delta + hR * (1.d0 - delta)
    b = bL * delta + bR * (1.d0 - delta)

    sh = shL * delta + shR * (1.d0 - delta)
    sb = sbL * delta + sbR * (1.d0 - delta)

    !--------------------------------------------------
    ! obtain macroscopic variables at interface
    !--------------------------------------------------
    ! conservative variables W_0 
    w(1) = sum(weight * h)
    w(2) = sum(weight * vn * h)
    w(3) = sum(weight * vt * h)
    w(4) = 0.5d0 * (sum(weight * (vn**2 + vt**2) * h) + sum(weight * b))

    ! convert to primary variables
    prim = conserve_prim_2d(w, gamma)

    ! heat flux
    qf = heat_flux_2f2v(h, b, prim, vn, vt, weight) 

    !--------------------------------------------------
    ! calculate a^L,a^R
    !--------------------------------------------------
    sw = (w - wL) / (0.5 * lenL)
    aL = micro_slope_2d(prim, sw, ink)

    sw = (wR - w) / (0.5 * lenR)
    aR = micro_slope_2d(prim, sw, ink)

    !--------------------------------------------------
    ! calculate time slope of W and A
    !--------------------------------------------------
    ! <u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
    call calc_moment_2d(prim, Mu, Mv, Mxi, Mu_L, Mu_R, ink) 

    Mau_L = moment_au_2d(aL, Mu_L, Mv, Mxi, 1, 0)
    Mau_R = moment_au_2d(aR, Mu_R, Mv, Mxi, 1, 0)

    sw = -prim(1) * (Mau_L + Mau_R)
    aT = micro_slope_2d(prim, sw, ink)

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_2d(prim, muref, omega)

    Mt(4) = tau*(1.d0 - exp(-dt / tau))
    Mt(5) = -tau * dt * exp(-dt / tau) + tau * Mt(4)
    Mt(1) = dt - Mt(4)
    Mt(2) = -tau * Mt(1) + Mt(5)
    Mt(3) = dt**2 / 2.d0 - tau * Mt(1)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Muv = moment_uv_2d(Mu, Mv, Mxi, 1, 0, 0)
    Mau_L = moment_au_2d(aL, Mu_L, Mv, Mxi, 2, 0)
    Mau_R = moment_au_2d(aR, Mu_R, Mv, Mxi, 2, 0)
    Mau_T = moment_au_2d(aT, Mu, Mv, Mxi, 1, 0)

    fluxw = Mt(1) * prim(1) * Muv + Mt(2) * prim(1) * (Mau_L + Mau_R) + Mt(3) * prim(1) * Mau_T

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_2f2v(H0, B0, prim, vn, vt, ink)

    ! Shakhov part H+ and B+
    call shakhov_2f2v(H0, B0, qf, prim, H_plus, B_plus, vn, vt, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * vn * H_plus) + &
               Mt(4) * sum(weight * vn * h) - Mt(5) * sum(weight * vn**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * vn * vn * H_plus) + &
               Mt(4) * sum(weight * vn * vn * h) - Mt(5) * sum(weight * vn * vn**2 * sh)
    fluxw(3) = fluxw(3) + Mt(1) * sum(weight * vt * vn * H_plus) + &
               Mt(4) * sum(weight * vt * vn * h) - Mt(5) * sum(weight * vt * vn**2 * sh)
    fluxw(4) = fluxw(4) + &
               Mt(1) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * H_plus) + sum(weight * vn * B_plus)) + &
               Mt(4) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * h) + sum(weight * vn * b)) - &
               Mt(5) * 0.5d0 * (sum(weight * vn**2 * (vn**2 + vt**2) * sh) + sum(weight * vn**2 * sb))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * vn * (H0 + H_plus) + &
            Mt(2) * vn**2 * (aL(1) * H0 + aL(2) * vn * H0 + aL(3) * vt * H0 + &
            0.5d0 * aL(4) * ((vn**2 + vt**2) * H0 + B0)) * delta + & 
            Mt(2) * vn**2 * (aR(1) * H0 + aR(2) * vn * H0 + aR(3) * vt * H0 + &
            0.5d0 * aR(4) * ((vn**2 + vt**2) * H0 + B0)) * (1.d0 - delta) + &
            Mt(3) * vn * (aT(1) * H0 + aT(2) * vn * H0 + aT(3) * vt * H0 + &
            0.5d0 * aT(4) * ((vn**2 + vt**2) * H0 + B0)) + &
            Mt(4) * vn * h - Mt(5) * vn**2 * sh

    fluxb = Mt(1) * vn * (B0 + B_plus) + &
            Mt(2) * vn**2 * (aL(1) * B0 + aL(2) * vn * B0 + aL(3) * vt * B0 + &
            0.5d0 * aL(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) * delta + &
            Mt(2) * vn**2 * (aR(1) * B0 + aR(2) * vn * B0 + aR(3) * vt * B0 + &
            0.5d0 * aR(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) * (1.d0 - delta) + &
            Mt(3) * vn * (aT(1) * B0 + aT(2) * vn * B0 + aT(3) * vt * B0 + &
            0.5d0 * aT(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) + &
            Mt(4) * vn * b - Mt(5) * vn**2 * sb

    !--------------------------------------------------
    ! final flux
    !--------------------------------------------------
    fluxw = lenFace * fluxw 
    fluxh = lenFace * fluxh
    fluxb = lenFace * fluxb

end
