module Peano

using GilbertCurves
using ComplexPaths
using Intervals

import ComplexPaths: pointsonpath

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

#! Is there a less awkward way to achieve this type of function construction? A macro potentially?

_linear_interp_lambda(startp,endp) = eval(Meta.parse("t ->  $startp + t*($(endp - startp))"))

function _build_hilbert_path(points) 

    indexes = points |> size |> gilbertindices

    li = indexes[1:end-1]
    ri = indexes[2:end]
    @assert length(li) == length(ri)

    indexpairs = [(li[i], ri[i]) for i ∈ eachindex(li)]
    
    paths = Vector{Path}()

    for ii ∈ indexpairs
        
        startp = points[ii[1]]
        endp = points[ii[2]]
    
        path = Path(_linear_interp_lambda(startp,endp),0..1)
        push!(paths,path)

    end

    out = PiecewisePath(paths...)

    return out
end

struct HilbertPath <: AbstractPath
    piecewise::PiecewisePath
    domain::Interval
    depth::Int
    subdivision::Int
    
    function HilbertPath(depth,min_coord,max_coord) 
        subdivision = 4
        plane = _genplane(min_coord,max_coord,depth*subdivision,depth*subdivision)
        piecewise = _build_hilbert_path(plane)
        new(    piecewise,
                piecewise.domain,
                depth,
                subdivision,
        )
    end
end

function ComplexPaths.pointsonpath(H::HilbertPath, n::Integer)
    pointsonpath(H.piecewise,n)
end 

end
