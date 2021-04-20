!-------------------------------------------------------------------
! Fast Spectral Method
!-------------------------------------------------------------------

subroutine boltzmann_fft(Q, f, kn, phi, psi, phipsi, unum, vnum, wnum, m)

    !include 'fftw3.f'
    include '/opt/intel/oneapi/mkl/latest/include/fftw/fftw3.f'
    
    integer, intent(in) :: unum, vnum, wnum, m
    real(kind=8), intent(inout) :: Q(unum, vnum, wnum)
    real(kind=8), intent(in) :: f(unum, vnum, wnum)
    real(kind=8), intent(in) :: kn
    real(kind=8), intent(in) :: phi(unum, vnum, wnum, m*(m-1)), psi(unum, vnum, wnum, m*(m-1))
    real(kind=8), intent(in) :: phipsi(unum, vnum, wnum)
    complex(kind=8) :: f_spec(unum, vnum, wnum), fc(unum, vnum, wnum), fc1(unum, vnum, wnum), fc2(unum, vnum, wnum), &
                       f_temp(unum, vnum, wnum), fc11(unum, vnum, wnum), fc22(unum, vnum, wnum)
    integer(kind=8) :: plan_backward, plan_backward2, plan_forward, plan_conv_f1, plan_conv_f2
    integer :: i
    
    ! initialize DFT
    ! forward(-1)/backward(+1) Fourier transformation
    call dfftw_plan_dft_3d(plan_backward, unum, vnum, wnum, f_spec, f_spec, FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_plan_dft_3d(plan_backward2, unum, vnum, wnum, f_temp, f_temp, FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_plan_dft_3d(plan_forward, unum, vnum, wnum, f_spec, fc, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_3d(plan_conv_f1, unum, vnum, wnum, fc1, fc11, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_3d(plan_conv_f2, unum, vnum, wnum, fc2, fc22, FFTW_FORWARD, FFTW_MEASURE)
    
    f_spec = cmplx(f, 0.0d0)
    call dfftw_execute_dft(plan_backward, f_spec, f_spec) ! backward
    f_spec = f_spec / (unum * vnum * wnum)
    
    ! shift zero-frequency component to center of spectrum
    call fftshift3d(f_spec)
    
    ! gain term
    f_temp = cmplx(0.0)
    do i=1,M*(M-1)
        fc1 = cmplx(f_spec * phi(:, :, :, i))
        fc2 = cmplx(f_spec * psi(:, :, :, i))
        call dfftw_execute_dft(plan_conv_f1, fc1, fc11)
        call dfftw_execute_dft(plan_conv_f2, fc2, fc22)
        !convolution
        f_temp = f_temp + fc11*fc22
    enddo
    
    ! loss term
    fc1 = cmplx(f_spec * phipsi)
    fc2 = f_spec
    call dfftw_execute_dft(plan_conv_f1, fc1, fc11)
    call dfftw_execute_dft(plan_conv_f2, fc2, fc22)
    
    ! gain-loss
    f_temp = f_temp - fc11 * fc22
    Q = 4.0d0 * PI**2 / kn / M**2 * real(f_temp)
    
    call dfftw_destroy_plan(plan_backward)
    call dfftw_destroy_plan(plan_backward2)
    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_conv_f1)
    call dfftw_destroy_plan(plan_conv_f2)
    
end