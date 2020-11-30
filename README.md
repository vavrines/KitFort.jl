# KitFort.jl

![CI](https://github.com/vavrines/KitFort.jl/workflows/CI/badge.svg)

The lightweight prototype with Fortran implementation of Kinetic.jl. It can be applied to cases with extreme pursuit of efficiency.

Compilation of dynamic Fortran library: `gfortran kitmod.f90 -o kitmod.so -shared -fPIC -O3 `
