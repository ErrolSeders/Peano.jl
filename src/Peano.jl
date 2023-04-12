module Peano

using GilbertCurves
using ComplexPaths
using Intervals

export HilbertPath, pointsonpath

function _genplane(min_coord::Complex, max_coord::Complex, width::Integer, height::Integer)
    real = range(min_coord.re, max_coord.re,length=width)
    imag = range(min_coord.im, max_coord.im,length=height)
    complexplane = zeros(Complex{Float64},(height, width))
    for (i,x) ∈ collect(enumerate(real))
        complexplane[:,i] .+= x
    end
    for (i,y) ∈ collect(enumerate(imag))
        complexplane[i,:] .+= (y * 1im)
    end
    reverse!(complexplane, dims=1)
    return complexplane
end

_linearinterplambda(startp,endp) = eval(Meta.parse("t ->  $startp + t*($(endp - startp))"))

function _build_piecewise_dict(points) 

    out = Dict{Interval, Function}()

    indexes = points |> size |> gilbertindices

    li = indexes[1:end-1]
    ri = indexes[2:end]
    @assert length(li) == length(ri)

    indexpairs = [(li[i], ri[i]) for i ∈ eachindex(li)]
    num_partitions = length(indexpairs)

    partition = range(0.0, 10.0, length = num_partitions) |> collect
    lb = partition[1:end-1]
    rb = partition[2:end]
    @assert length(lb) == length(rb)
    bounds = [(lb[i], rb[i]) for i ∈ eachindex(lb)]

    for (i,b) ∈ enumerate(bounds[1:end-1])
        startp = points[indexpairs[i][1]]
        endp = points[indexpairs[i][2]]

        pair = Interval{Closed,Open}(b[1],b[2]) => _linearinterplambda(startp,endp)
        push!(out, pair)
    end

    laststartp = points[indexpairs[end][1]]
    lastendpdiff = points[indexpairs[end][2]] - laststartp
    last_pair = bounds[end][1]..bounds[end][2] => (t -> laststartp + t*(lastendpdiff))
    push!(out,last_pair)

    @assert foldl(intersect, out |> keys) |> isempty

    return out
end

function _call_hilbertpath(t::Real, D::AbstractDict)
    for (key,value) ∈ D
        if t ∈ key
            t_map = (t - key.first) / (key.last - key.first)
            f = value
            return f(t_map)
        end
    end
end

struct HilbertPath <: AbstractPath
    parameterization::Function
    start::Real
    ending::Real
    piecewise::Dict{Interval,Function}
    depth::Int
    subdivision::Int
    
    function HilbertPath(depth,min_coord,max_coord) 
        subdivision = 4
        plane = _genplane(min_coord,max_coord,depth*subdivision,depth*subdivision)
        dict = _build_piecewise_dict(plane)
        _call_hilbert_internal(t) = _call_hilbertpath(t,dict)
        new(    _call_hilbert_internal,
                0.0,
                1.0,
                dict,
                depth,
                subdivision,
        )
    end
end

function pointsonpath(H::HilbertPath, n::Integer)
    @assert n > 0 "n must be a positive non-zero integer"
    return [P.parameterization(t,H.piecewise) for t in range(P.start, P.ending, length=n)]
end 
    
end
