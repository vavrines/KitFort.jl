!-------------------------------------------------------------------
! Fortran auxiliary module for Kinetic.jl
!-------------------------------------------------------------------

module Kinetic

    implicit none

    real(kind=8), parameter :: PI = 3.141592654d0
    integer, parameter :: MNUM = 6 !number of normal moments
    integer, parameter :: MTUM = 6 !number of tangential moments

    contains

    include "method_1d.f90"
    include "method_2d.f90"
    include "kfvs.f90"
    include "ugks.f90"
    include "kcu.f90"
    !include "fsm.f90"

end