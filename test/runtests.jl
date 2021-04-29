using Test
import KitFort

include("test_flux_1d.jl")

KitFort.gfortran("flux")
KitFort.ifort("flux")