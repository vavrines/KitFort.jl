!-------------------------------------------------------------------
! Kinetic Flux Vector Splitting Method
!-------------------------------------------------------------------

subroutine flux_kfvs_1f1v(fluxw, fluxh, hL, shL, lenL, hR, shR, lenR, &
                          unum, uspace, weight, dt)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum) !// interface fluxes
    real(kind=8), intent(in) :: hL(unum), shL(unum) !// distribution functions and their slopes in left cell
    real(kind=8), intent(in) :: hR(unum), shR(unum) !// distribution functions and their slopes in right cell
    real(kind=8), intent(in) :: lenL, lenR !// cell lengths
    real(kind=8), intent(in) :: uspace(unum), weight(unum) !// velocity quadrature points and weights
    real(kind=8), intent(in) :: dt !// time step

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    ! interface variable
    real(kind=8) :: h(unum)
    real(kind=8) :: sh(unum)

    ! time integral
    real(kind=8) :: Mt(2)

    ! heaviside step function
    real(kind=8) :: delta(unum)
    delta = (sign(1.d0, uspace) + 1.d0) / 2.d0

    !--------------------------------------------------
    ! upwind reconstruction
    !--------------------------------------------------
    h = (hL + 0.5 * lenL * shL) * delta + &
        (hR - 0.5 * lenR * shR) * (1.d0 - delta)
    sh = shL * delta + shR * (1.d0 - delta)

    Mt(1) = dt
    Mt(2) = 0.5d0 * dt**2

    !--------------------------------------------------
    ! fluxes
    !--------------------------------------------------
    ! macro flux
    fluxw(1) = Mt(1) * sum(weight * uspace * h) - Mt(2) * sum(weight * uspace**2 * sh)
    fluxw(2) = Mt(1) * sum(weight * uspace**2 * h) - Mt(2) * sum(weight * uspace**3 * sh)
    fluxw(3) = Mt(1) * 0.5d0 * sum(weight * uspace * uspace**2 * h) - &
               Mt(2) * 0.5d0 * sum(weight * uspace**2 * uspace**2 * sh)

    ! micro flux
    fluxh = Mt(1) * uspace * h - Mt(2) * uspace**2 * sh

end

!-------------------------------------------------------------------

subroutine flux_kfvs_2f1v(fluxw, fluxh, fluxb, hL, bL, shL, sbL, lenL, hR, bR, shR, sbR, lenR, &
                          unum, uspace, weight, dt)

    integer, intent(in) :: unum !// number of velocity grids
    real(kind=8), intent(inout) :: fluxw(3), fluxh(unum), fluxb(unum) !// interface fluxes
    real(kind=8), intent(in) :: hL(unum), bL(unum), shL(unum), sbL(unum) !// distribution functions and their slopes in left cell
    real(kind=8), intent(in) :: hR(unum), bR(unum), shR(unum), sbR(unum) !// distribution functions and their slopes in right cell
    real(kind=8), intent(in) :: lenL, lenR !// cell lengths
    real(kind=8), intent(in) :: uspace(unum), weight(unum) !// velocity quadrature points and weights
    real(kind=8), intent(in) :: dt !// time step

    !--------------------------------------------------
    ! initialize
    !--------------------------------------------------
    ! interface variable
    real(kind=8) :: h(unum)
    real(kind=8) :: b(unum)
    real(kind=8) :: sh(unum)
    real(kind=8) :: sb(unum)

    ! time integral
    real(kind=8) :: Mt(2)

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

    Mt(1) = dt
    Mt(2) = 0.5d0 * dt**2

    !--------------------------------------------------
    ! fluxes
    !--------------------------------------------------
    ! macro flux
    fluxw(1) = Mt(1) * sum(weight * uspace * h) - Mt(2) * sum(weight * uspace**2 * sh)
    fluxw(2) = Mt(1) * sum(weight * uspace**2 * h) - Mt(2) * sum(weight * uspace**3 * sh)
    fluxw(3) = Mt(1) * 0.5d0 * (sum(weight * uspace * uspace**2 * h) + sum(weight * uspace * b)) - &
               Mt(2) * 0.5d0 * (sum(weight * uspace**2 * uspace**2 * sh) + sum(weight * uspace**2 * sb))

    ! micro flux
    fluxh = Mt(1) * uspace * h - Mt(2) * uspace**2 * sh
    fluxb = Mt(1) * uspace * b - Mt(2) * uspace**2 * sb

end