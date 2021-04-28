module KitFort

export flux_kfvs!,
       flux_kcu!,
       flux_ugks!

include("kfvs.jl")
include("kcu.jl")
include("ugks.jl")

COMPILER = :gfortran

gfortran(file::AbstractString) = "__kinetic_MOD_" * file

ifort(file::AbstractString) = "kinetic_mp_" * file * "_"

end # module
