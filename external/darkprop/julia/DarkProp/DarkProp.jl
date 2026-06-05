module DarkProp

using HDF5
using StatsBase
using KernelDensity
using LinearAlgebra
using Interpolations
using MultiQuad: quad
using NumericalTools

using CxxWrap

if haskey(ENV, "DARKPROP_JULIA_PATH")
    @wrapmodule(() -> joinpath(ENV["DARKPROP_JULIA_PATH"], "libdarkpropjl"))
else
    @wrapmodule(() -> "libdarkpropjl")
end

function __init__()
    @initcxx
end

using HEPUnits; const units = HEPUnits

const rEarth = 6371units.km

export bl_to_xyz, hist_estimate, spectrumz, intensity_isoring

function bl_to_xyz(bl::Matrix)
    xyz = Array{eltype(bl)}(undef, size(bl)[1], 3)
    b = @views bl[:, 1]
    l = @views bl[:, 2]
    xyz[:, 1] = @. cos(b) * cos(l)
    xyz[:, 2] = @. cos(b) * sin(l)
    xyz[:, 3] = @. sin(b)
    return xyz
end


function hist_estimate(x; ws=nothing, nbins=201, method="linear")
    if method == "linear"
        bins = linspace(minimum(x), maximum(x), nbins)
        centers = bins[1:end-1] + diff(bins) / 2
    elseif method == "sqrt"
        bins = linspace(sqrt(minimum(x)), sqrt(maximum(x)), nbins).^2
        centers = bins[1:end-1] + diff(bins) / 2
    elseif method == "log"
        bins = geomspace(minimum(x), maximum(x), nbins)
        centers = sqrt.(bins[1:end-1] .* bins[2:end])
    end
    centers[begin] = bins[begin]
    centers[end] = bins[end]

    if isnothing(ws)
        ws = ones(size(x))
    end

    h = fit(Histogram, x, weights(ws), bins)
    normalize!(h, mode=:pdf)

    if method == "log"
        return loginterpolator(centers, h.weights)
    end
    return linear_interpolation(centers, h.weights, extrapolation_bc=0)
end

function spectrumz(filename, groupname; method="histogram", nbins=50, kwargs...)
    f = h5open(filename, "r")
    R0 = read(f, "R0")
    norm = read(f, "norm")

    group = f[groupname]
    Rd = read(group, "R")

    T::Vector{Float64} = read(group, "T")
    p::Matrix{Float64} = permutedims(read(group, "p"), (2, 1))
    r::Matrix{Float64} = permutedims(read(group, "r"), (2, 1))
    w::Vector{Float64} = read(group, "weight")
    Nsim = read(group, "Nsim")
    Nsample = read(group, "Nsample")

    valid = isfinite.(T)
    T = T[valid]
    p = p[valid, :]
    r = r[valid, :]
    w = w[valid]

    if Nsample != length(T)
        error("Nsample = $Nsample, Nsim = $Nsim")
    end

    ep = bl_to_xyz(p)
    er = bl_to_xyz(r)
    abscos = abs.(sum(ep .* er; dims=2))
    ws = w ./ abscos
    total_num = sum(ws)
    println("$filename, length(T) = $(length(T)), Nsim = $Nsim, total_num = $total_num")
    println("raw weights sum: $(sum(w))")
    if length(T) < 10
        return E -> zero(E)
    end

    if method == "kde"
        logkde = InterpKDE(kde(log.(T); weights=weights(ws), kwargs...))
        function kdewrapper(E)
            return pdf(logkde, log(E)) / E * total_num * R0^2 * norm / (4 * Rd^2 * Nsim)
        end
        return kdewrapper
    elseif method == "histogram"
        dens = hist_estimate(T; ws=ws, nbins=nbins, method="log", kwargs...)
        function histwrapper(E)
            return dens(E) * total_num * R0^2 * norm / (4 * Rd^2 * Nsim)
        end
        return histwrapper
    end
    error("Unknown method: $method")
end


function intensity_isoring(filename; nbins=100, check=false, kwargs...)
    fluxes = []
    println(filename)
    f = h5open(filename, "r")
    R0 = read(f, "R0")
    S0 = 4π * R0^2
    norm = read(f, "norm")
    depth = read(f, "depth")
    ring_bins = []
    for di in 1:length(depth)
        fluxes_d = []
        dgroup = f["depth_$(di-1)"]
        ring_bins = read(dgroup, "ring_bins")
        Rd = read(dgroup, "R")

        for ri in 1:length(ring_bins)-1
            rgroup = dgroup["ring_$(ri-1)"]

            T::Vector{Float64} = read(rgroup, "T")
            p::Matrix{Float64} = permutedims(read(rgroup, "p"), (2, 1))
            r::Matrix{Float64} = permutedims(read(rgroup, "r"), (2, 1))
            w::Vector{Float64} = read(rgroup, "weight")
            Nsim = read(rgroup, "Nsim")
            Nsample = read(rgroup, "Nsample")

            valid = isfinite.(T)
            if !all(valid)
                T = T[valid]
                p = p[valid, :]
                r = r[valid, :]
                w = w[valid]
            end

            if Nsample != length(T)
                error("Nsample = $Nsample, Nsim = $Nsim")
            end

            delta_S = 2π * Rd^2 * (cos(ring_bins[ri]) - cos(ring_bins[ri + 1]))

            if length(T) < 10
                push!(fluxes_d, E -> zero(E))
                continue
            end

            ep = bl_to_xyz(p)
            er = bl_to_xyz(r)

            abscos = abs.(sum(ep .* er; dims=2))
            ws = w ./ abscos
            total_num = sum(ws)
            println("$filename, length(T) = $(length(T)), Nsim = $Nsim, total_num = $total_num")


            dens = hist_estimate(T; ws=ws, nbins=nbins, method="log", kwargs...)
            if check
                println("delta_S = $(delta_S)\nnormalized? "
                    * "$(quad(lnT -> exp(lnT) * dens(exp(lnT)),
                         log(minimum(T)), log(maximum(T)), atol=0))")
            end

            function wrapper(E, dens=dens, total_num=total_num, delta_S=delta_S, Nsim=Nsim)
                return dens(E)*total_num*S0*norm / (16π*delta_S*Nsim)
            end

            push!(fluxes_d, wrapper)
        end
        push!(fluxes, fluxes_d)
    end
    close(f)
    return fluxes, ring_bins, depth
end

end  # module DarkProp
