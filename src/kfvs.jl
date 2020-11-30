#--- 1d1f1v ---#
function flux_kfvs!(
    fw::X,
    fh::Y,
    hL::Z,
    hR::Z,
    u::A,
    ω::A,
    nu,
    dt,
    shL = zeros(eltype(hL), axes(hL))::Z,
    shR = zeros(eltype(hR), axes(hR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
    A<:AbstractArray{<:AbstractFloat,1},
}
    cd(@__DIR__)
    ccall(
        (:__kinetic_MOD_flux_kfvs_1f1v, "fortran/kitmod.so"),
        Nothing,
        (Ref{Float64}, Ref{Float64},
        Ref{Float64}, Ref{Float64},
        Ref{Float64}, Ref{Float64},
        Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        fw,
        fh,
        hL,
        hR,
        u,
        ω,
        nu,
        dt,
        shL,
        shR,
    )
end

#--- 1d2f1v ---#
function flux_kfvs!(
    fw::X,
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    ω::A,
    nu,
    dt,
    shL = zeros(eltype(hL), axes(hL))::Z,
    sbL = zeros(eltype(bL), axes(bL))::Z,
    shR = zeros(eltype(hR), axes(hR))::Z,
    sbR = zeros(eltype(bR), axes(bR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
    A<:AbstractArray{<:AbstractFloat,1},
}
    cd(@__DIR__)
    ccall(
        (:__kinetic_MOD_flux_kfvs_2f1v, "fortran/kitmod.so"),
        Nothing,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, 
        Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        Ref{Float64}, Ref{Float64},
        Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        fw,
        fh,
        fb,
        hL,
        bL,
        hR,
        bR,
        u,
        ω,
        nu,
        dt,
        shL,
        sbL,
        shR,
        sbR,
    )
end