module BoucheTrous

using AbstractFFTs: fftfreq, rfftfreq
using CartesianBoxes: CartesianBox
using LinOps: LinOpDFT, LinOpDiag, LinOpGrad, LinOpSelect
using OptimPackNextGen: vmlmb, AutoDiffObjectiveFunction
using Zygote: Zygote
using DifferentiationInterface: AutoZygote, prepare_gradient
using StatsBase: mad, median
using ImageFiltering: mapwindow


export bouchetrous!, redetectbadpix!

"""
    bouchetrous!(img, badpix; pupil_radius=nothing, kwargs...)

Interpolate bad pixels in each frame of a three-dimensional image stack,
modifying `img` in place. Entries where `badpix` is `true` are corrected;
other entries are left unchanged.

`img` and `badpix` must have the same size and one-based indexing. When
`pupil_radius` is `nothing`, a gradient regularizer is used; otherwise, the
value selects a Fourier-domain regularizer. Additional keyword arguments are
passed to `OptimPackNextGen.vmlmb`.

Returns `nothing`.
"""
function bouchetrous!(img::AbstractArray{T, N}, badpix::AbstractArray{Bool, N}; pupil_radius = nothing, kwds...) where {T, N}
    Base.require_one_based_indexing(img, badpix)

    sz = size(img)

    F = make_function(T, sz[1:2], pupil_radius)

    return Threads.@threads for n in axes(img, 3)
        data = img[:, :, n]
        S = LinOpSelect(badpix[:, :, n])
        p0 = S * data
        data[badpix[:, :, n]] .= 0
        cost(x) = sum(abs2, F * (S' * x + data))
        prep = prepare_gradient(cost, AutoZygote(), S * data)
        fg! = AutoDiffObjectiveFunction(cost, prep, AutoZygote())
        p = vmlmb(fg!, p0; lower = minimum(data), kwds...)
        img[badpix[:, :, n]] .= p
    end
end


function make_function(::Type{T}, sz::NTuple{2, Int}, pupil_radius::AbstractFloat) where {T}
    pupil = sqrt.(rfftfreq(sz[1]) .^ 2 .+ (fftfreq(sz[2]) .^ 2)') .< pupil_radius
    F = LinOpDFT(T, sz[1:2])
    D = LinOpDiag(.!pupil)
    return D * F
end

make_function(::Type{T}, sz::NTuple{2, Int}, ::Nothing) where {T} = LinOpGrad(sz)


"""
    get_box(badpix)

Return the smallest `CartesianBox` enclosing the unmasked pixels in the
three-dimensional Boolean array `badpix`. A pixel is unmasked when its value
is `false`; the bounds are computed across all slices in the third dimension.
At least one unmasked pixel is required.
"""
function get_box(badpix)
    bx = findfirst(n -> !all(view(badpix, n, :, :)), axes(badpix, 1))
    ex = findlast(n -> !all(view(badpix, n, :, :)), axes(badpix, 1))
    by = findfirst(n -> !all(view(badpix, :, n, :)), axes(badpix, 2))
    ey = findlast(n -> !all(view(badpix, :, n, :)), axes(badpix, 2))
    return CartesianBox(bx:ex, by:ey)
end


## Detect and interpolate spurious bad pixel using median filter
function redetectbadpix!(data, weights)
    med = mapwindow(median, data, (3, 3, 1))
    medw = mapwindow(median, weights, (3, 3, 1))
    P = (med .- data) .* (medw .+ 1.0e-6)
    md = mad(P)
    badpix = .!((-10 * md) .< P .< (10 * md))
    data[badpix] .= med[badpix]
    return weights[badpix] .= 0
end

end # module BoucheTrous
