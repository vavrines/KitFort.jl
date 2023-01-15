# KitFort.jl

[![version](https://juliahub.com/docs/KitFort/version.svg)](https://juliahub.com/ui/Packages/KitFort/2mlJf)
![CI](https://img.shields.io/github/actions/workflow/status/vavrines/KitFort.jl/main.yml?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://xiaotianbai.com/Kinetic.jl/stable/)
[![codecov](https://img.shields.io/codecov/c/github/vavrines/KitFort.jl)](https://codecov.io/gh/vavrines/KitFort.jl)
[![deps](https://juliahub.com/docs/KitBase/deps.svg)](https://juliahub.com/ui/Packages/KitBase/YOFTS?t=2)

This lightweight module provides the Fortran backends in [Kinetic.jl](https://github.com/vavrines/Kinetic.jl) ecosystem.
It's not included in the main module by default, and can be manually imported in the extreme pursuit of efficiency.
[Check the documentation](https://xiaotianbai.com/Kinetic.jl/dev/) for information on the implementation and use of the package.

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
