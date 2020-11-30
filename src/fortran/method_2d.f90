!-------------------------------------------------------------------
! 2D Methods
!-------------------------------------------------------------------

function global_frame(w, cosx, cosy)

    real(kind=8), intent(in) :: w(4)
    real(kind=8), intent(in) :: cosx, cosy
    real(kind=8) :: global_frame(4)

    global_frame(1) = w(1)
    global_frame(2) = w(2) * cosx - w(3) * cosy
    global_frame(3) = w(2) * cosy + w(3) * cosx
    global_frame(4) = w(4)
    
end function global_frame

!-------------------------------------------------------------------

function local_frame(w, cosx, cosy)

    real(kind=8), intent(in) :: w(4)
    real(kind=8), intent(in) :: cosx, cosy
    real(kind=8) :: local_frame(4)

    local_frame(1) = w(1)
    local_frame(2) = w(2) * cosx + w(3) * cosy
    local_frame(3) = w(3) * cosx - w(2) * cosy
    local_frame(4) = w(4)
    
end function local_frame

!-------------------------------------------------------------------

function conserve_prim_2d(w, gamma)

    real(kind=8), intent(in) :: w(4), gamma
    real(kind=8) :: conserve_prim_2d(4)

    conserve_prim_2d(1) = w(1)
    conserve_prim_2d(2) = w(2) / w(1)
    conserve_prim_2d(3) = w(3) / w(1)
    conserve_prim_2d(4) = 0.5d0 * w(1) / (gamma - 1.d0) / (w(4) - 0.5d0 * (w(2)**2 + w(3)**2) / w(1))

end

!-------------------------------------------------------------------

subroutine maxwell_2f2v(h, b, prim, vn, vt, ink)

    real(kind=8), dimension(:,:), intent(out) :: h, b
    real(kind=8), dimension(:,:), intent(in) :: vn, vt
    real(kind=8), intent(in) :: prim(4), ink

    h = prim(1) * (prim(4) / PI) * exp(-prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2))
    b = h * ink / (2.d0 * prim(4))
    
end

!-------------------------------------------------------------------

subroutine shakhov_2f2v(H, B, qf, prim, H_plus, B_plus, vn, vt, ink, pr)
            
    real(kind=8), dimension(:,:), intent(in) :: H, B
    real(kind=8), intent(in) :: qf(2)
    real(kind=8), intent(in) :: prim(4)
    real(kind=8), dimension(:,:), intent(out) :: H_plus, B_plus
    real(kind=8), dimension(:,:), intent(in) :: vn, vt
    real(kind=8), intent(in) :: ink, pr

    H_plus = 0.8d0 * (1.d0 - pr) * prim(4)**2 / prim(1) * &
             ((vn - prim(2)) * qf(1) + (vt - prim(3)) * qf(2)) * &
             (2.d0 * prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2) + ink - 5.d0) * H
    B_plus = 0.8d0 * (1.d0 - pr) * prim(4)**2 / prim(1) * &
             ((vn - prim(2)) * qf(1) + (vt - prim(3)) * qf(2)) * &
             (2.d0 * prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2) + ink - 3.d0) * B
            
end

!-------------------------------------------------------------------

function collision_time_2d(prim, muref, omega)

    real(kind=8), intent(in) :: prim(4), muref, omega
    real(kind=8) :: collision_time_2d

    collision_time_2d = muref * 2.d0 * prim(4)**(1.d0 - omega) / prim(1)

end

!-------------------------------------------------------------------

function heat_flux_2f2v(h, b, prim, vn, vt, weight)

    real(kind=8), dimension(:,:), intent(in) :: h, b
    real(kind=8), dimension(:,:), intent(in) :: vn, vt, weight
    real(kind=8), intent(in) :: prim(4)
    real(kind=8) :: heat_flux_2f2v(2)

    heat_flux_2f2v(1) = 0.5d0 * (sum(weight * (vn - prim(2)) * ((vn - prim(2))**2 + (vt - prim(3))**2) * h) + &
                        sum(weight * (vn - prim(2)) * b))
    heat_flux_2f2v(2) = 0.5d0 * (sum(weight * (vt - prim(3)) * ((vn - prim(2))**2 + (vt - prim(3))**2) * h) + &
                        sum(weight * (vt - prim(3)) * b))
    
end

!-------------------------------------------------------------------

function micro_slope_2d(prim, sw, ink)
            
    real(kind=8), intent(in) :: prim(4), sw(4), ink
    real(kind=8) :: micro_slope_2d(4)

    micro_slope_2d(4) = 4.d0 * prim(4)**2 / (ink + 2.d0) / prim(1) * &
                        (2.d0 * sw(4) - 2.d0 * prim(2) * sw(2) - 2.d0 * prim(3) * sw(3) + &
                        sw(1) * (prim(2)**2 + prim(3)**2 - 0.5d0 * (ink + 2.d0) / prim(4)))
    micro_slope_2d(3) = 2.d0 * prim(4) / prim(1) * (sw(3) - prim(3) * sw(1)) - prim(3) * micro_slope_2d(4)
    micro_slope_2d(2) = 2.d0 * prim(4) / prim(1) * (sw(2) - prim(2) * sw(1)) - prim(2) * micro_slope_2d(4)
    micro_slope_2d(1) = sw(1) / prim(1) - prim(2) * micro_slope_2d(2) - prim(3) * micro_slope_2d(3) - &
                        0.5d0 * (prim(2)**2 + prim(3)**2 + 0.5d0 * (ink + 2.d0) / prim(4)) * micro_slope_2d(4)

end

!-------------------------------------------------------------------

subroutine calc_moment_2d(prim, Mu, Mv, Mxi, Mu_L, Mu_R, ink)

    real(kind=8), intent(in) :: prim(4)
    real(kind=8), intent(out) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM)
    real(kind=8), intent(out) :: Mv(0:MTUM)
    real(kind=8), intent(out) :: Mxi(0:2)
    real(kind=8), intent(in) :: ink
    integer :: i

    ! moments of normal velocity
    Mu_L(0) = 0.5d0 * erfc(-sqrt(prim(4)) * prim(2))
    Mu_L(1) = prim(2) * Mu_L(0) + 0.5d0 * exp(-prim(4) * prim(2)**2) / sqrt(PI * prim(4))
    Mu_R(0) = 0.5d0 * erfc(sqrt(prim(4)) * prim(2))
    Mu_R(1) = prim(2) * Mu_R(0) - 0.5d0 * exp(-prim(4) * prim(2)**2) / sqrt(PI * prim(4))

    do i=2,MNUM
        Mu_L(i) = prim(2) * Mu_L(i-1) + 0.5d0 * (i-1) * Mu_L(i-2) / prim(4)
        Mu_R(i) = prim(2) * Mu_R(i-1) + 0.5d0 * (i-1) * Mu_R(i-2) / prim(4)
    end do

    Mu = Mu_L + Mu_R

    ! moments of tangential velocity
    Mv(0) = 1.0
    Mv(1) = prim(3)

    do i=2,MTUM
        Mv(i) = prim(3) * Mv(i-1) + 0.5d0 * (i-1) * Mv(i-2) / prim(4)
    end do

    ! moments of \xi
    Mxi(0) = 1.d0 !<\xi^0>
    Mxi(1) = 0.5d0 * ink / prim(4) !<\xi^2>
    Mxi(2) = (ink**2 + 2.d0 * ink) / (4.d0 * prim(4)**2) !<\xi^4>
    
end

!-------------------------------------------------------------------

function moment_uv_2d(Mu, Mv, Mxi, alpha, beta, delta)

    real(kind=8), intent(in) :: Mu(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    integer, intent(in) :: alpha, beta, delta
    real(kind=8) :: moment_uv_2d(4)

    moment_uv_2d(1) = Mu(alpha) * Mv(beta) * Mxi(delta/2)
    moment_uv_2d(2) = Mu(alpha+1) * Mv(beta) * Mxi(delta/2)
    moment_uv_2d(3) = Mu(alpha) * Mv(beta+1) * Mxi(delta/2)
    moment_uv_2d(4) = 0.5d0 * (Mu(alpha+2) * Mv(beta) * Mxi(delta/2) + &
                      Mu(alpha) * Mv(beta+2) * Mxi(delta/2) + &
                      Mu(alpha) * Mv(beta) * Mxi((delta+2)/2))

end

!-------------------------------------------------------------------

function moment_au_2d(a, Mu, Mv, Mxi, alpha, beta)
    
    real(kind=8), intent(in) :: a(4)
    real(kind=8), intent(in) :: Mu(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    integer, intent(in) :: alpha, beta
    real(kind=8) :: moment_au_2d(4)

    moment_au_2d = a(1) * moment_uv_2d(Mu, Mv, Mxi, alpha+0, beta+0, 0) + &
                   a(2) * moment_uv_2d(Mu, Mv, Mxi, alpha+1, beta+0, 0) + &
                   a(3) * moment_uv_2d(Mu, Mv, Mxi, alpha+0, beta+1, 0) + &
                   0.5d0 * a(4) * moment_uv_2d(Mu, Mv, Mxi, alpha+2, beta+0, 0) + &
                   0.5d0 * a(4) * moment_uv_2d(Mu, Mv, Mxi, alpha+0, beta+2, 0) + &
                   0.5d0 * a(4) * moment_uv_2d(Mu, Mv, Mxi, alpha+0, beta+0, 2)
    
end