!-------------------------------------------------------------------
! 1D Methods
!-------------------------------------------------------------------

function conserve_prim_1d(w, gamma)

    real(kind=8), intent(in) :: w(3), gamma
    real(kind=8) :: conserve_prim_1d(3)

    conserve_prim_1d(1) = w(1)
    conserve_prim_1d(2) = w(2) / w(1)
    conserve_prim_1d(3) = 0.5d0 * w(1) / (gamma - 1.d0) / (w(3) - 0.5d0 * w(2)**2 / w(1))

end

!-------------------------------------------------------------------

subroutine maxwell_1f1v(H, prim, uspace)

    real(kind=8), dimension(:), intent(out) :: H
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace

    H = prim(1) * (prim(3) / PI)**(1.d0/2.d0) * exp(-prim(3) * (uspace - prim(2))**2)

end

!-------------------------------------------------------------------

subroutine maxwell_2f1v(H, B, prim, uspace, ink)

    real(kind=8), dimension(:), intent(out) :: H, B
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace
    real(kind=8), intent(in) :: ink

    H = prim(1) * (prim(3) / PI)**(1.d0/2.d0) * exp(-prim(3) * (uspace - prim(2))**2)
    B = h * ink / (2.d0 * prim(3))

end

!-------------------------------------------------------------------

subroutine shakhov_1f1v(H, qf, prim, H_plus, uspace, ink, prandtl)

    real(kind=8), dimension(:), intent(in) :: H
    real(kind=8), intent(in) :: qf
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(out) :: H_plus
    real(kind=8), dimension(:), intent(in) :: uspace
    real(kind=8), intent(in) :: ink
    real(kind=8), intent(in) :: prandtl
    
    H_plus = 0.8d0 * (1.d0 - prandtl) * prim(3)**2 / prim(1) * &
             (uspace - prim(2)) * qf * (2.d0 * prim(3) * (uspace - prim(2))**2 + ink - 5.d0) * H

end

!-------------------------------------------------------------------

subroutine shakhov_2f1v(H, B, qf, prim, H_plus, B_plus, uspace, ink, prandtl)

    real(kind=8), dimension(:), intent(in) :: H, B
    real(kind=8), intent(in) :: qf
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(out) :: H_plus, B_plus
    real(kind=8), dimension(:), intent(in) :: uspace
    real(kind=8), intent(in) :: ink
    real(kind=8), intent(in) :: prandtl
    
    H_plus = 0.8d0 * (1.d0 - prandtl) * prim(3)**2 / prim(1) * &
             (uspace - prim(2)) * qf * (2.d0 * prim(3) * (uspace - prim(2))**2 + ink - 5.d0) * H
    B_plus = 0.8d0 * (1.d0 - prandtl) * prim(3)**2 / prim(1) * &
             (uspace - prim(2)) * qf * (2.d0 * prim(3) * (uspace - prim(2))**2 + ink - 3.d0) * B

end

!-------------------------------------------------------------------

function collision_time_1d(prim, muref, omega)

    real(kind=8), intent(in) :: prim(3), muref, omega
    real(kind=8) :: collision_time_1d

    collision_time_1d = muref * 2.d0 * prim(3)**(1.d0 - omega) / prim(1)

end

!-------------------------------------------------------------------

function heat_flux_1f1v(h, prim, uspace, weight)
    
    real(kind=8), dimension(:), intent(in) :: h
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace, weight
    real(kind=8) :: heat_flux_1f1v

    heat_flux_1f1v = 0.5d0 * sum(weight * (uspace - prim(2)) * (uspace - prim(2))**2 * h)

end

!-------------------------------------------------------------------

function heat_flux_2f1v(h, b, prim, uspace, weight)
    
    real(kind=8), dimension(:), intent(in) :: h, b
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace, weight
    real(kind=8) :: heat_flux_2f1v

    heat_flux_2f1v = 0.5d0 * (sum(weight * (uspace - prim(2)) * (uspace - prim(2))**2 * h) + &
                     sum(weight * (uspace - prim(2)) * b))

end

!-------------------------------------------------------------------

function micro_slope_1d(prim, sw, ink)

    real(kind=8), intent(in) :: prim(3), sw(3), ink
    real(kind=8) :: micro_slope_1d(3)

    micro_slope_1d(3) = 4.d0 * prim(3)**2 / (ink + 1.d0) / prim(1) * &
                        (2.d0 * sw(3) - 2.d0 * prim(2) * sw(2) + sw(1) * (prim(2)**2 - 0.5d0 * (ink + 1.d0) / prim(3)))
    micro_slope_1d(2) = 2.d0 * prim(3) / prim(1) * (sw(2) - prim(2) * sw(1)) - prim(2) * micro_slope_1d(3)
    micro_slope_1d(1) = sw(1) / prim(1) - prim(2) * micro_slope_1d(2) - &
                        0.5d0 * (prim(2)**2 + 0.5d0 * (ink + 1.d0) / prim(3)) * micro_slope_1d(3)

end

!-------------------------------------------------------------------

subroutine calc_moment_1d(prim, Mu, Mxi, Mu_L, Mu_R, inK)

    real(kind=8), intent(in) :: prim(3)
    real(kind=8), intent(out) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM)
    real(kind=8), intent(out) :: Mxi(0:2)
    real(kind=8), intent(in) :: inK
    integer :: i

    !moments of normal velocity
    Mu_L(0) = 0.5d0 * erfc(-sqrt(prim(3)) * prim(2))
    Mu_L(1) = prim(2) * Mu_L(0) + 0.5d0 * exp(-prim(3) * prim(2)**2) / sqrt(PI * prim(3))
    Mu_R(0) = 0.5d0 * erfc(sqrt(prim(3)) * prim(2))
    Mu_R(1) = prim(2) * Mu_R(0) - 0.5d0 * exp(-prim(3) * prim(2)**2) / sqrt(PI * prim(3))

    do i=2,MNUM
        Mu_L(i) = prim(2) * Mu_L(i-1) + 0.5d0 * (i-1) * Mu_L(i-2) / prim(3)
        Mu_R(i) = prim(2) * Mu_R(i-1) + 0.5d0 * (i-1) * Mu_R(i-2) / prim(3)
    end do

    Mu = Mu_L + Mu_R

    !moments of \xi
    Mxi(0) = 1.0 !<\xi^0>
    Mxi(1) = 0.5 * inK / prim(3) !<\xi^2>
    Mxi(2) = (inK**2 + 2.d0 * inK) / (4.d0 * prim(3)**2) !<\xi^4>

end

!-------------------------------------------------------------------

function moment_uv_1d(Mu, Mxi, alpha, delta)

    real(kind=8), intent(in) :: Mu(0:MNUM), Mxi(0:2)
    integer, intent(in) :: alpha, delta
    real(kind=8) :: moment_uv_1d(3)

    moment_uv_1d(1) = Mu(alpha) * Mxi(delta/2)
    moment_uv_1d(2) = Mu(alpha+1) * Mxi(delta/2)
    moment_uv_1d(3) = 0.5d0 * (Mu(alpha+2) * Mxi(delta/2) + Mu(alpha) * Mxi((delta+2)/2))

end

!-------------------------------------------------------------------

function moment_au_1d(a, Mu, Mxi, alpha)

    real(kind=8), intent(in) :: a(3)
    real(kind=8), intent(in) :: Mu(0:MNUM), Mxi(0:2)
    integer, intent(in) :: alpha
    real(kind=8) :: moment_au_1d(3)

    moment_au_1d = a(1) * moment_uv_1d(Mu, Mxi, alpha+0, 0) + &
                   a(2) * moment_uv_1d(Mu, Mxi, alpha+1, 0) + &
                   0.5d0 * a(3) * moment_uv_1d(Mu, Mxi, alpha+2, 0) + &
                   0.5d0 * a(3) * moment_uv_1d(Mu, Mxi, alpha+0, 2)

end