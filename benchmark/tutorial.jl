# Tutorial for calling 

# test function
x = 1.0
y = 2.0
println(ccall(
    (:__kit_MOD_add, "kitmod.so"),
    Float64,
    (Ptr{Float64}, Ptr{Float64}),
    Ref(x),
    Ref(y),
))
println(ccall(
    (:__kit_MOD_add, "kitmod.so"),
    Float64,
    (Ref{Float64}, Ref{Float64}),
    x,
    y,
))

# test subroutine
sf = zeros(3)
prim = Float64.([1.0, 0.0, 1.0])
sw = Float64.([1.0, 1.0, 1.0])
K = 2
ccall(
    (:__kit_MOD_micro_slope, "kitmod.so"),
    Nothing,
    (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    sf,
    prim,
    sw,
    K,
)
println(sf)

# for julia, it's not possible to use `dimension(:)`
# you have to specify concrete dimension args
ccall(
    (:__kit_MOD_flex_slope, "kitmod.so"),
    Cvoid,
    (Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    3,
    sf,
    prim,
    sw,
    K,
)
println(sf)

ccall(
    (:__kit_MOD_direct_slope, "kitmod.so"),
    #Float64,
    Array{Float64,1},
    (Ref{Float64}, Ref{Float64}, Ref{Float64}),
    prim,
    sw,
    K,
)