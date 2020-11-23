# ------------------------------------------
# Adjoint julia script
# ------------------------------------------
using Kinetic, BenchmarkTools

u = collect(-5.0:0.05:5.0)
nu = length(u)
weights = ones(nu) .* 0.5

fw = zeros(3)
fh = zeros(nu)
fb = zeros(nu)

inK = 2
γ = 5.0 / 3.0
primL = [1., 0., 1.]
wL = prim_conserve(primL, γ)
hL = maxwellian(u, primL) |> Array;
bL = hL .* 2 ./ (2.)
shL = zeros(nu)
sbL = zeros(nu)
lenL = 0.1

primR = [0.5, 0., 1.]
wR = prim_conserve(primR, γ)
hR = maxwellian(u, primR) |> Array;
bR = hR .* 2 ./ (2.)
shR = zeros(nu)
sbR = zeros(nu)
lenR = 0.1

muref = 0.001
omega = 0.72
prandtl = 1.0
dt = 1e-4

@btime ccall(
    (:__kit_MOD_flux_ugks1d, "./fortran/kitmod.so"),
    Nothing,
    (Ref{Float64}, Ref{Float64}, Ref{Float64}, 
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
    Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    fw,
    fh,
    fb,
    wL,
    hL,
    bL,
    shL,
    sbL,
    lenL,
    wR,
    hR,
    bR,
    shR,
    sbR,
    lenR,
    nu,
    u,
    weights,
    inK,
    γ,
    muref,
    omega,
    prandtl,
    dt,
)
println(fw)

@btime flux_ugks!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, weights, inK, γ, muref, omega, prandtl, dt, lenL, lenR, shL, sbL, shR, sbR)