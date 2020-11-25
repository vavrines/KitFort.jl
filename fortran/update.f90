!-------------------------------------------------------------------
! Update Algorithm
!-------------------------------------------------------------------

subroutine update_1d2f1v(fwL, fhL, fbL, w, prim, h, b, fwR, fhR, fbR, unum, u, weights, &
                         inK, gamma, muref, omega, pr, dx, dt, res, avg)

    integer, intent(in):: unum
    real(kind=8), intent(in) :: fwL(3), fhL(unum), fbL(unum)
    real(kind=8), intent(inout) :: w(3), h(unum), b(unum)
    real(kind=8), intent(in) :: fwR(3), fhR(unum), fbR(unum)
    real(kind=8), intent(in) :: u(unum), weights(unum)
    real(kind=8), intent(in) :: inK, gamma
    real(kind=8), intent(in) :: muref, omega, pr
    real(kind=8), intent(in) :: dx, dt
    real(kind=8), intent(inout) :: res(3), avg(3)

    ! t=t^n
    real(kind=8) :: MH_old(unum), MB_old(unum)
    real(kind=8) :: w_old(3), prim_old(3)
    real(kind=8) :: tau_old

    ! t=t^n+1
    real(kind=8) :: MH(unum), MB(unum)
    real(kind=8) :: MH_plus(unum), MB_plus(unum)

    real(kind=8) :: prim(3)
    real(kind=8) :: qf
    real(kind=8) :: tau

    !--------------------------------------------------
    ! store W^n and calculate H^n,B^n,\tau^n
    !--------------------------------------------------
    w_old = w
    prim_old = conserve_prim_1d(w_old, gamma)

    call maxwell_2f1v(MH_old, MB_old, prim_old, u, inK)
    tau_old = collision_time_1d(prim_old, muref, omega)

    !--------------------------------------------------
    ! update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
    !--------------------------------------------------
    w = w + (fwL - fwR) / dx
    prim = conserve_prim_1d(w, gamma)

    call maxwell_2f1v(MH, MB, prim, u, inK)
    tau = collision_time_1d(prim, muref, omega)

    !--------------------------------------------------
    ! record residual
    !--------------------------------------------------
    res = res + (w_old - w)**2
    avg = avg + abs(w)

    !--------------------------------------------------
    ! Shakhov part
    !--------------------------------------------------
    ! heat flux at t=t^n
    qf = heat_flux_2f1v(h, b, prim_old, u, weights)

    ! h^+ = H+H^+ at t=t^n
    call shakhov_part(MH_old, MB_old, qf, prim_old, MH_plus, MB_plus, u, weights, inK, pr)
    MH_old = MH_old + MH_plus !h^+
    MB_old = MB_old + MB_plus !b^+

    ! h^+ = H+H^+ at t=t^{n+1}
    call shakhov_part(MH, MB, qf, prim, MH_plus, MB_plus, u, weights, inK, pr)
    MH = MH + MH_plus
    MB = MB + MB_plus

    !--------------------------------------------------
    ! update distribution function
    !--------------------------------------------------
    h = (h + (fhL - fhR) / dx + &
        0.5 * dt * (MH / tau + (MH_old - h) / tau_old))/(1.0 + 0.5 * dt / tau)
    b = (b + (fbL - fbR) / dx + &
        0.5 * dt * (MB / tau + (MB_old - b) / tau_old))/(1.0 + 0.5 * dt / tau)

end