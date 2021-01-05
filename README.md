# KitFort.jl

![CI](https://github.com/vavrines/KitFort.jl/workflows/CI/badge.svg)

This package serves as a lightweight module of Fortran methods in [Kinetic.jl](https://github.com/vavrines/Kinetic.jl) ecosystem. It can be applied to the cases with extreme pursuit of efficiency.

## Compilation of dynamic library
Execute `gfortran kitmod.f90 -o kitmod.so -shared -fPIC -O3 `
