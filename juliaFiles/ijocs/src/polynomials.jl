using LinearAlgebra
using Polynomials
import Base.push!
import Base.in
import Base.intersect


abstract type AbstractPoly end
getbounds(p::AbstractPoly) = p.bounds
polyval(p::AbstractPoly, x::Number) = x in getbounds(p) ? _polyval(p,x) : println("raise out of bounds error")

#abstract type AbstractPolyMethod end

struct Bound{T<:Number} # possible expansion, expand to 3 bound type for 3 different endpoint states, then use traits to set lower and upper bounds, saves on storage space and speed
    lower::T
    upper::T
    includeL::Bool
    includeU::Bool
    function Bound{T}(l,u,incL, incU) where {T<:Number}
        if l==u && !(incL&&incU)
            println("raise impossible bound error")
            return
        end
        return l ≤ u ? new(l,u, incL, incU) : println("raise lower higher than upper error")
    end
end

Bound(l::T,u::T; includeLower=true, includeUpper=true) where {T<:Number} = Bound{T}(l,u, includeLower, includeUpper)
Bound(b) = Bound(b,b) # equal upper and lower bounds
Base.show(io::IO, b::Bound) = print( b.includeL ? "[" : "(", b.lower, ",", b.upper, b.includeU ? "]" : ")")
Base.show(io::IO, ::MIME"text/plain", b::Bound{T}) where{T} = print(io, "Bound{$T}\n", b) # this output is broken for arrays, fix it at some point

function in(x, b::Bound)
    b.lower < x < b.upper ||
    b.includeL && b.lower == x ||
    b.includeU && b.upper == x
end

function intersect(b1::Bound, b2::Bound) # return 1 if b1 higher than b2, -1 if b1 lower than b2
    if lower(b1) ≥ upper(b2) && !(includesLower(b1) && includesUpper(b2)) || lower(b1) > upper(b2)
        return 1 # b1 higher than b2
    end
    if lower(b2) ≥ upper(b1) && !(includesLower(b2) && includesUpper(b1)) || lower(b2) > upper(b1)
        return -1 # b1 higher than b2
    end
    l = max(lower(b1), lower(b2))
    u = min(upper(b1), upper(b2))
    return Bound(l,u,includeLower = l in b1 && l in b2, includeUpper = u in b1 && u in b2)
end

lower(b::Bound) = b.lower
upper(b::Bound) = b.upper
includesUpper(b::Bound) = upper(b) in b
includesLower(b::Bound) = lower(b) in b


mutable struct CoefPoly <: AbstractPoly
    poly::Poly
    bounds::Bound
end
create_coef_poly(coefs::Vector) = CoefPoly(Poly(coefs), Bound(-Inf, Inf))
create_coef_poly(coefs::Vector, b::Bound) = CoefPoly(Poly(coefs), b)
create_coef_poly(coefs::Vector, lower, upper) = CoefPoly(Poly(coefs), Bound(lower, upper))
_polyval(p::CoefPoly, x) = polyval(p.poly, x)


mutable struct GlobalPoly <: AbstractPoly
    poly::Vector{AbstractPoly}#array of polynomials
    bounds::Vector{Bound}
    GlobalPoly(p::AbstractPoly, b::Bound) = new([p], [b]) # new initializing globalPoly must be created with an initial poly to avoid checks further down the line
end

function push!(gloPol::GlobalPoly, locPol::AbstractPoly, gloBound::Bound) # check to make sure bounds of lower are same or less restricted than global bounds
    if !includesLower(getbounds(locPol)) #check first that local and global bound endpoints are compatible
        if includesLower(gloBound)
            return println("raise lower bound incompattible error")
        end
    elseif !includesUpper(getbounds(locPol))
        if includesUpper(gloBound)
            return println("raise upper bound incompattible error")
        end
    end
    for i in 1:length(gloPol.bounds)
        if !(typeof(intersect(gloBound, gloPol.bounds[i]))<:Number)
            return println("raise bounds intersect with existing poly error")
        end
    end
    push!(gloPol.poly, locPol)
    push!(gloPol.bounds, gloBound)
end
push!(gloPol::GlobalPoly, locPol::AbstractPoly, gloLower, gloUpper) = push!(gloPol, locPol, Bound(gloLower, gloUpper)) 
push!(gloPol::GlobalPoly, locPol::AbstractPoly) = push!(gloPol, locPol, getbounds(locPol))

function getbounds(p::GlobalPoly)
    low,up = promote(lower(p.bounds[1]), upper(p.bounds[end]))
    return Bound(low, up, includeLower=includesLower(p.bounds[1]), includeUpper=includesUpper(p.bounds[end]))
end

function _polyval(p::GlobalPoly, x)
    for i in 1:length(p.bounds)
        if x in p.bounds[i]
            return _polyval(p.poly[i], global_to_local(getbounds(p.poly[i]), p.bounds[i], x))
        end
    end
    return println("raise not found in set")
end

# following Julia convert() convention where thing being converted to comes first
local_to_global(gloB::Bound, locB::Bound, x::Number) = lower(gloB) + (upper(gloB) - lower(gloB)) * (x - lower(locB)) / ( upper(locB) - lower(locB) )
global_to_local(locB::Bound, gloB::Bound, x::Number) = lower(locB) + (upper(locB) - lower(locB)) * (x - lower(gloB)) / ( upper(gloB) - lower(gloB) )


mutable struct LagrangePoly <: AbstractPoly
    x::Vector
    y::Vector
    weights::Vector
    bounds::Bound
end

function lagrange_bary_weights(x)
    numPoints = length(x)
    weights = [prod(1/(x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints]
    return weights
end

function lagrange_eval_weights(xeval, x, y, w)
    num = mapreduce((x,y,w)-> w*y / (xeval - x),+,x,y,w)
    denom = mapreduce((x,w)-> w / (xeval - x),+,x,w)
    val = num/denom
    return isnan(val) ? y[x.==xeval] : val
end

create_lagrange_poly(x, y, lower, upper) = LagrangePoly(x, y, lagrange_bary_weights(x), Bound(lower, upper ))
create_lagrange_poly(x, y, b::Bound) = LagrangePoly(x, y, lagrange_bary_weights(x), b)
_polyval(p::LagrangePoly,x) = lagrange_eval_weights(x, p.x, p.y, p.weights)


mutable struct LagrangePolyLite <: AbstractPoly # using less space as not including bounds
    x::Vector
    y::Vector
    weights::Vector
end
create_lagrange_poly(x, y) = LagrangePolyLite(x, y, lagrange_bary_weights(x))
getbounds(p::LagrangePolyLite) = Bound(min(p.x...), max(p.x...))
_polyval(p::LagrangePolyLite,x) = lagrange_eval_weights(x, p.x, p.y, p.weights)


# LGR points

function lgr_points(order) # does not include endpoints
    if order>1
        N = 1:order-2
        a = zeros(order-1)
        b = zeros(order-2)
        a[1] = 1 / 3 
        a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
        b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
        J = SymTridiagonal(a,b)
        return pushfirst!(eigvals(J),-1)
    end
end

## temp code
g = GlobalPoly()
p = create_lagrange_poly([1,2,3],[1,2,3].^2)
push!(g,p, Bound(0,1))
push!(g,p, 1.0001,2.0)
polyval(g, 1.5)
