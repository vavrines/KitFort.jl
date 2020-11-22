!-------------------------------------------------------------------
! Fortran auxiliary module for Kinetic.jl
!-------------------------------------------------------------------

module KIT

implicit none

real(kind=8), parameter :: PI = 3.141592654d0

integer, parameter :: MNUM = 6 !number of normal moments
integer, parameter :: MTUM = 6 !number of tangential moments

contains

!-------------------------------------------------------------------

subroutine flux_ugks1d(fluxw, fluxh, fluxb, wL, hL, bL, shL, sbL, lenL, wR, hR, bR, shR, sbR, lenR, &
                       unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)

integer, intent(in):: unum !// number of velocity grids
real(kind=8), intent(inout) :: fluxw(3), fluxh(unum), fluxb(unum) !// interface fluxes
real(kind=8), intent(in) :: wL(3), wR(3) !// conservative variables
real(kind=8), intent(in) :: hL(unum), bL(unum), shL(unum), sbL(unum) !// distribution functions and their slopes in left cell
real(kind=8), intent(in) :: hR(unum), bR(unum), shR(unum), sbR(unum) !// distribution functions and their slopes in right cell
real(kind=8), intent(in) :: lenL, lenR !// cell lengths
real(kind=8), intent(in) :: uspace(unum), weight(unum) !// velocity quadrature points and weights
real(kind=8), intent(in) :: ink, gamma !// internal degrees of freedom of gas and Poisson ratio
real(kind=8), intent(in) :: muref, omega, prandtl !// reference viscosity, VHS model index, and Prandtl number
real(kind=8), intent(in) :: dt !// time step

!Heaviside step function
integer, allocatable, dimension(:) :: delta

!interface variable
real(kind=8) :: h(unum), b(unum)
real(kind=8) :: H0(unum), B0(unum)
real(kind=8) :: H_plus(unum), B_plus(unum)
real(kind=8) :: sh(unum), sb(unum)
real(kind=8) :: w(3), prim(3)
real(kind=8) :: qf
real(kind=8) :: sw(3)
real(kind=8) :: aL(3), aR(3), aT(3)

!moments variable
real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mv(0:MTUM), Mxi(0:2)
real(kind=8) :: Mau_0(3), Mau_L(3), Mau_R(3), Mau_T(3)
real(kind=8) :: tau
real(kind=8) :: Mt(5)

!--------------------------------------------------
! initialize
!--------------------------------------------------
!Heaviside step function
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
!conservative variables W_0 
w(1) = sum(weight * h)
w(2) = sum(weight * uspace * h)
w(3) = 0.5d0 * (sum(weight * uspace**2 * h) + sum(weight * b))

!convert to primary variables
prim = prim_variable(w, gamma)

!heat flux
qf = heat_flux(h, b, prim, uspace, weight) 

!--------------------------------------------------
! calculate a^L,a^R
!--------------------------------------------------
sw = (w - wL) / (0.5 * lenL)
aL = micro_slope(prim, sw, ink)

sw = (wR - w) / (0.5 * lenR)
aR = micro_slope(prim, sw, ink)

!--------------------------------------------------
! calculate time slope of W and A
!--------------------------------------------------
!<u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
call calc_moment_u(prim, Mu, Mxi, Mu_L, Mu_R, ink) 

Mau_L = moment_au(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
Mau_R = moment_au(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

sw = -prim(1) * (Mau_L + Mau_R) !time slope of W
aT = micro_slope(prim, sw, ink) !calculate A

!--------------------------------------------------
! calculate collision time and some time integration terms
!--------------------------------------------------
tau = collision_time(prim, muref, omega)

Mt(4) = tau * (1.d0 - exp(-dt / tau))
Mt(5) = -tau * dt * exp(-dt / tau) + tau * Mt(4)
Mt(1) = dt - Mt(4)
Mt(2) = -tau * Mt(1) + Mt(5) 
Mt(3) = dt**2 / 2.d0 - tau * Mt(1)

!--------------------------------------------------
! calculate the flux of conservative variables related to g0
!--------------------------------------------------
Mau_0 = moment_uv(Mu, Mxi, 1, 0) !<u*\psi>
Mau_L = moment_au(aL, Mu_L, Mxi, 2) !<aL*u^2*\psi>_{>0}
Mau_R = moment_au(aR, Mu_R, Mxi, 2) !<aR*u^2*\psi>_{<0}
Mau_T = moment_au(aT, Mu, Mxi, 1) !<A*u*\psi>

fluxw = Mt(1) * prim(1) * Mau_0 + Mt(2) * prim(1) * (Mau_L + Mau_R) + Mt(3) * prim(1) * Mau_T

!--------------------------------------------------
! calculate the flux of conservative variables related to g+ and f0
!--------------------------------------------------
!Maxwellian distribution H0 and B0
call maxwell(H0, B0, prim, uspace, ink)

!Shakhov part H+ and B+
call shakhov(H0, B0, qf, prim, H_plus, B_plus, uspace, ink, prandtl)

!macro flux related to g+ and f0
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

subroutine maxwell(H, B, prim, uspace, ink)

    real(kind=8), dimension(:), intent(out) :: H, B
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace
    real(kind=8), intent(in) :: ink

    H = prim(1) * (prim(3) / PI)**(1.d0/2.d0) * exp(-prim(3) * (uspace - prim(2))**2)
    B = h * ink / (2.d0 * prim(3))

end subroutine maxwell

!-------------------------------------------------------------------

subroutine shakhov(H, B, qf, prim, H_plus, B_plus, uspace, ink, prandtl)

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

end subroutine shakhov

!-------------------------------------------------------------------

function collision_time(prim, muref, omega)

    real(kind=8), intent(in) :: prim(3), muref, omega
    real(kind=8) :: collision_time

    collision_time = muref * 2.d0 * prim(3)**(1.d0 - omega) / prim(1)

end function collision_time

!-------------------------------------------------------------------

function prim_variable(w, gamma)

    real(kind=8), intent(in) :: w(3), gamma
    real(kind=8) :: prim_variable(3)

    prim_variable(1) = w(1)
    prim_variable(2) = w(2) / w(1)
    prim_variable(3) = 0.5d0 * w(1) / (gamma - 1.d0) / (w(3) - 0.5d0 * w(2)**2 / w(1))

end function prim_variable

!-------------------------------------------------------------------

function heat_flux(h, b, prim, uspace, weight)
    
    real(kind=8), dimension(:), intent(in) :: h, b
    real(kind=8), intent(in) :: prim(3)
    real(kind=8), dimension(:), intent(in) :: uspace, weight
    real(kind=8) :: heat_flux

    heat_flux = 0.5d0 * (sum(weight * (uspace - prim(2)) * (uspace - prim(2))**2 * h) + sum(weight * (uspace - prim(2)) * b))

end function heat_flux

!-------------------------------------------------------------------

function micro_slope(prim, sw, ink)

    real(kind=8), intent(in) :: prim(3), sw(3), ink
    real(kind=8) :: micro_slope(3)

    micro_slope(3) = 4.d0 * prim(3)**2 / (ink + 1.d0) / prim(1) * &
                     (2.d0 * sw(3) - 2.d0 * prim(2) * sw(2) + sw(1) * (prim(2)**2 - 0.5d0 * (ink + 1.d0) / prim(3)))
    micro_slope(2) = 2.d0 * prim(3) / prim(1) * (sw(2) - prim(2) * sw(1)) - prim(2) * micro_slope(3)
    micro_slope(1) = sw(1) / prim(1) - prim(2) * micro_slope(2) - &
                     0.5d0 * (prim(2)**2 + 0.5d0 * (ink + 1.d0) / prim(3)) * micro_slope(3)

end function micro_slope

!-------------------------------------------------------------------

subroutine calc_moment_u(prim, Mu, Mxi, Mu_L, Mu_R, inK)

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

    Mu = Mu_L+Mu_R

    !moments of \xi
    Mxi(0) = 1.0 !<\xi^0>
    Mxi(1) = 0.5 * inK / prim(3) !<\xi^2>
    Mxi(2) = (inK**2 + 2.d0 * inK) / (4.d0 * prim(3)**2) !<\xi^4>

end subroutine calc_moment_u

!-------------------------------------------------------------------

function moment_uv(Mu, Mxi, alpha, delta)

    real(kind=8), intent(in) :: Mu(0:MNUM), Mxi(0:2)
    integer, intent(in) :: alpha, delta
    real(kind=8) :: moment_uv(3)

    moment_uv(1) = Mu(alpha) * Mxi(delta/2)
    moment_uv(2) = Mu(alpha+1) * Mxi(delta/2)
    moment_uv(3) = 0.5d0 * (Mu(alpha+2) * Mxi(delta/2) + Mu(alpha) * Mxi((delta+2)/2))

end function moment_uv

!-------------------------------------------------------------------

function moment_au(a, Mu, Mxi, alpha)

    real(kind=8), intent(in) :: a(3)
    real(kind=8), intent(in) :: Mu(0:MNUM), Mxi(0:2)
    integer, intent(in) :: alpha
    real(kind=8) :: moment_au(3)

    moment_au = a(1) * moment_uv(Mu,Mxi,alpha+0,0) + &
                a(2) * moment_uv(Mu,Mxi,alpha+1,0) + &
                0.5d0 * a(3) * moment_uv(Mu,Mxi,alpha+2,0) + &
                0.5d0 * a(3) * moment_uv(Mu,Mxi,alpha+0,2)

end function moment_au

!-------------------------------------------------------------------

end