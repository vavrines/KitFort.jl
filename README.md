# KitFort.jl

![CI](https://github.com/vavrines/KitFort.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/vavrines/KitFort.jl/branch/main/graph/badge.svg?token=67tfVc3AtW)](https://codecov.io/gh/vavrines/KitFort.jl)

This package serves as a lightweight module of Fortran methods in [Kinetic.jl](https://github.com/vavrines/Kinetic.jl) ecosystem. 
It's not included in the main module by default, and can be manually imported in the cases with extreme pursuit of efficiency.

## Dynamic library

The modern Fortran methods is provided by the shared library [kitmod.so](https://github.com/vavrines/KitFort.jl/blob/main/src/fortran/kitmod.so) and called from Julia with the help of `ccall` function.

## Recompilation

To generate the dynamic library file to be called from Julia, make sure the GNU Fortran compiler has been installed in the computer.
```bash
gfortran kitmod.f90 -o kitmod.so -shared -fPIC -O3
```
Alternatively, the Intel Fortran compiler `ifort` can be employed with the same command above.
Note that GNU and Intel compilers present slightly different behaviors on the function call.
For example, the low-level KFVS flux function takes:
- `:__kinetic_MOD_flux_kfvs_1f1v` for GNU
- `:kinetic_mp_flux_kfvs_1f1v_` for Intel
Please don't do the recompilation unless you're sure what's exactly going on.
