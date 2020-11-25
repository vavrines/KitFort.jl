!-------------------------------------------------------------------
! Unified Gas Kinetic Scheme
!-------------------------------------------------------------------

subroutine flux_kcu_1f1v(fluxw, fluxh, hL, shL, lenL, hR, shR, lenR, &
                         unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum) !// interface fluxes
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

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mxi(0:2)
    real(kind=8) :: Mau_0(3)
    real(kind=8) :: tau
    real(kind=8) :: Mt(3)

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

    ! moments
    call calc_moment_1d(prim, Mu, Mxi, Mu_L, Mu_R, ink) 

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_1d(prim, muref, omega)

    Mt(2) = tau * (1.d0 - exp(-dt / tau))
    Mt(3) = -tau * dt * exp(-dt / tau) + tau * Mt(2)
    Mt(1) = dt - Mt(2)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Mau_0 = moment_uv_1d(Mu, Mxi, 1, 0) !<u*\psi>
    fluxw = Mt(1) * prim(1) * Mau_0

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_1f1v(H0, prim, uspace)

    ! Shakhov part H+ and B+
    call shakhov_1f1v(H0, qf, prim, H_plus, uspace, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * uspace * H_plus) + &
               Mt(2) * sum(weight * uspace * h) - Mt(3) * sum(weight * uspace**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * uspace**2 * H_plus) + &
               Mt(2) * sum(weight * uspace**2 * h) - Mt(3) * sum(weight * uspace**3 * sh)
    fluxw(3) = fluxw(3) + &
               Mt(1) * 0.5d0 * sum(weight * uspace * uspace**2 * H_plus) + &
               Mt(2) * 0.5d0 * sum(weight * uspace * uspace**2 * h) - &
               Mt(3) * 0.5d0 * sum(weight * uspace**2 * uspace**2 * sh)

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * uspace * (H0 + H_plus) + &
            Mt(2) * uspace * h - Mt(3) * uspace**2 * sh

end

!-------------------------------------------------------------------

subroutine flux_kcu_2f1v(fluxw, fluxh, fluxb, hL, bL, shL, sbL, lenL, hR, bR, shR, sbR, lenR, &
                         unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum), fluxb(unum) !// interface fluxes
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

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mxi(0:2)
    real(kind=8) :: Mau_0(3)
    real(kind=8) :: tau
    real(kind=8) :: Mt(3)

    ! heaviside step function
    real(kind=8) :: delta(unum)
    delta = (sign(1.d0, uspace) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! upwind reconstruction
    !--------------------------------------------------
    h = (hL + 0.5 * lenL * shL) * delta + &
        (hR - 0.5 * lenR * shR) * (1.d0 - delta)
    b = (bL + 0.5 * lenL * sbL) * delta + &
        (bR - 0.5 * lenR * sbR) * (1.d0 - delta)

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

    ! moments
    call calc_moment_1d(prim, Mu, Mxi, Mu_L, Mu_R, ink) 

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_1d(prim, muref, omega)

    Mt(2) = tau * (1.d0 - exp(-dt / tau))
    Mt(3) = -tau * dt * exp(-dt / tau) + tau * Mt(2)
    Mt(1) = dt - Mt(2)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Mau_0 = moment_uv_1d(Mu, Mxi, 1, 0) !<u*\psi>
    fluxw = Mt(1) * prim(1) * Mau_0

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_2f1v(H0, B0, prim, uspace, ink)

    ! Shakhov part H+ and B+
    call shakhov_2f1v(H0, B0, qf, prim, H_plus, B_plus, uspace, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * uspace * H_plus) + &
               Mt(2) * sum(weight * uspace * h) - Mt(3) * sum(weight * uspace**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * uspace**2 * H_plus) + &
               Mt(2) * sum(weight * uspace**2 * h) - Mt(3) * sum(weight * uspace**3 * sh)
    fluxw(3) = fluxw(3) + &
               Mt(1) * 0.5d0 * (sum(weight * uspace * uspace**2 * H_plus) + sum(weight * uspace * B_plus)) + &
               Mt(2) * 0.5d0 * (sum(weight * uspace * uspace**2 * h) + sum(weight * uspace * b)) - &
               Mt(3) * 0.5d0 * (sum(weight * uspace**2 * uspace**2 * sh) + sum(weight * uspace**2 * sb))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * uspace * (H0 + H_plus) + &
            Mt(2) * uspace * h - Mt(3) * uspace**2 * sh
    fluxb = Mt(1) * uspace * (B0 + B_plus) + &
            Mt(2) * uspace * b - Mt(3) * uspace**2 * sb

end

!-------------------------------------------------------------------

subroutine flux_kcu_2f2v(fluxw, fluxh, fluxb, &
                         hL, bL, shL, sbL, lenL, &
                         hR, bR, shR, sbR, lenR, &
                         unum, vnum, uspace, vspace, weight, &
                         ink, gamma, muref, omega, prandtl, &
                         dt, lenFace, cosa, sina)

    integer, intent(in) :: unum, vnum
    real(kind=8), intent(inout) :: fluxw(4), fluxh(unum, vnum), fluxb(unum, vnum)
    real(kind=8), intent(in) :: hL(unum, vnum), bL(unum, vnum), shL(unum, vnum), sbL(unum, vnum)
    real(kind=8), intent(in) :: hR(unum, vnum), bR(unum, vnum), shR(unum, vnum), sbR(unum, vnum)
    real(kind=8), intent(in) :: lenL, lenR
    real(kind=8), intent(in) :: uspace(unum, vnum), vspace(unum, vnum), weight(unum, vnum)
    real(kind=8), intent(in) :: ink, gamma
    real(kind=8), intent(in) :: muref, omega, prandtl
    real(kind=8), intent(in) :: dt, lenFace
    real(kind=8), intent(in) :: cosa, sina

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    real(kind=8) :: delta(unum, vnum)
    real(kind=8) :: vn(unum, vnum), vt(unum, vnum)

    ! interface variable
    real(kind=8) :: h(unum, vnum), b(unum, vnum)
    real(kind=8) :: H0(unum, vnum), B0(unum, vnum)
    real(kind=8) :: H_plus(unum, vnum), B_plus(unum, vnum)
    real(kind=8) :: sh(unum, vnum), sb(unum, vnum)
    real(kind=8) :: w(4), prim(4)
    real(kind=8) :: qf(2)

    ! moments variable
    real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    real(kind=8) :: Muv(4)
    real(kind=8) :: tau
    real(kind=8) :: Mt(3)

    vn = uspace * cosa + vspace * sina
    vt = vspace * cosa - uspace * sina

    delta = (sign(1.d0, vn) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! reconstruct initial distribution
    !--------------------------------------------------
    h = (hL + 0.5 * lenL * shL) * delta + &
        (hR - 0.5 * lenR * shR) * (1.d0 - delta)
    b = (bL + 0.5 * lenL * sbL) * delta + &
        (bR - 0.5 * lenR * sbR) * (1.d0 - delta)

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

    ! moments
    call calc_moment_2d(prim, Mu, Mv, Mxi, Mu_L, Mu_R, ink) 

    !--------------------------------------------------
    ! calculate collision time and some time integration terms
    !--------------------------------------------------
    tau = collision_time_2d(prim, muref, omega)

    Mt(2) = tau * (1.d0 - exp(-dt / tau))
    Mt(3) = -tau * dt * exp(-dt / tau) + tau * Mt(2)
    Mt(1) = dt - Mt(2)

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g0
    !--------------------------------------------------
    Muv = moment_uv_2d(Mu, Mv, Mxi, 1, 0, 0)
    fluxw = Mt(1) * prim(1) * Muv

    !--------------------------------------------------
    ! calculate the flux of conservative variables related to g+ and f0
    !--------------------------------------------------
    ! Maxwellian distribution H0 and B0
    call maxwell_2f2v(H0, B0, prim, vn, vt, ink)

    ! Shakhov part H+ and B+
    call shakhov_2f2v(H0, B0, qf, prim, H_plus, B_plus, vn, vt, ink, prandtl)

    ! macro flux related to g+ and f0
    fluxw(1) = fluxw(1) + Mt(1) * sum(weight * vn * H_plus) + &
               Mt(2) * sum(weight * vn * h) - Mt(3) * sum(weight * vn**2 * sh)
    fluxw(2) = fluxw(2) + Mt(1) * sum(weight * vn * vn * H_plus) + &
               Mt(2) * sum(weight * vn * vn * h) - Mt(3) * sum(weight * vn * vn**2 * sh)
    fluxw(3) = fluxw(3) + Mt(1) * sum(weight * vt * vn * H_plus) + &
               Mt(2) * sum(weight * vt * vn * h) - Mt(3) * sum(weight * vt * vn**2 * sh)
    fluxw(4) = fluxw(4) + &
               Mt(1) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * H_plus) + sum(weight * vn * B_plus)) + &
               Mt(2) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * h) + sum(weight * vn * b)) - &
               Mt(3) * 0.5d0 * (sum(weight * vn**2 * (vn**2 + vt**2) * sh) + sum(weight * vn**2 * sb))

    !--------------------------------------------------
    ! calculate flux of distribution function
    !--------------------------------------------------
    fluxh = Mt(1) * vn * (H0 + H_plus) + &
            Mt(2) * vn * h - Mt(3) * vn**2 * sh

    fluxb = Mt(1) * vn * (B0 + B_plus) + &
            Mt(2) * vn * b - Mt(3) * vn**2 * sb

    !--------------------------------------------------
    ! final flux
    !--------------------------------------------------
    ! convert to global frame
    fluxw = global_frame(fluxw, cosa, sina)

    ! total flux
    fluxw = lenFace * fluxw 
    fluxh = lenFace * fluxh
    fluxb = lenFace * fluxb

end
