begin
    u = collect(-5.0:0.1:5.0)
    nu = length(u)
    weights = ones(nu) .* 0.5

    fw = zeros(3)
    fh = zeros(nu)
    fb = zeros(nu)

    inK = 2
    γ = 5.0 / 3.0

    primL = [1., 0., 1.]
    wL = [
        primL[1],
        primL[1] * primL[2],
        0.5 * primL[1] / primL[3] / (γ - 1.0) + 0.5 * primL[1] * primL[2]^2,
    ]
    hL = @. primL[1] * sqrt(primL[3] / π) * exp(-primL[3] * (u - primL[2])^2)
    bL = hL .* 2 ./ (2.)
    shL = zeros(nu)
    sbL = zeros(nu)
    lenL = 0.1

    primR = [0.5, 0., 1.]
    wR = [
        primR[1],
        primR[1] * primR[2],
        0.5 * primR[1] / primR[3] / (γ - 1.0) + 0.5 * primR[1] * primR[2]^2,
    ]
    hR = @. primR[1] * sqrt(primR[3] / π) * exp(-primR[3] * (u - primR[2])^2)
    bR = hR .* 2 ./ (2.)
    shR = zeros(nu)
    sbR = zeros(nu)
    lenR = 0.1

    muref = 0.001
    omega = 0.72
    prandtl = 1.0
    dt = 1e-4
end
using KitFort
flux_kfvs!(
    fw,
    fh,
    hL,
    hR,
    u,
    weights,
    nu,
    dt,
    shL,
    shR,
)

flux_kfvs!(
    fw,
    fh,
    fb,
    hL,
    bL,
    hR,
    bR,
    u,
    weights,
    nu,
    dt,
    shL,
    sbL,
    shR,
    sbR,
)