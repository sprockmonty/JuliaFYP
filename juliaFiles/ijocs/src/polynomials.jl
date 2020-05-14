using LinearAlgebra
using Polynomials
#using Reduce # for genvandpoly macro 
#Reduce.Preload() # required otherwise Reduce has a fit
import Base.push!
import Base.in
import Base.intersect
import Polynomials.polyval

# specifiers for the createpoly function
#abstract type AbstractPolySpecifier end
#struct LGRSpecifier <: AbstractPolySpecifier end
#struct LagrangeSpecifier <: AbstractPolySpecifier end
#struct CoefSpecifier <: AbstractPolySpecifier end
#struct VandSpecifier <: AbstractPolySpecifier end
# struct GlobalSpecifier <: AbstractPolySpecifier end

# abstract poly types and functions
#abstract type AbstractPoly end
#abstract type AbstractLagrangePoly<:AbstractPoly end
#getbounds(p::AbstractPoly) = p.bounds
#polyval(p::AbstractPoly, x::Number) = x in getbounds(p) ? _polyval(p,x) : println("raise out of bounds error")
#polyval(p::AbstractPoly, x::Number, y) = x in getbounds(p) ? _polyval(p,x,y) : println("raise out of bounds error")
#polyval!(p::AbstractPoly, x::Number, y) = x in getbounds(p) ? _polyval!(p,x,y) : println("raise out of bounds error") # get value at x and update poly with new y vals
#polyder(p::AbstractPoly, der::Number) = x in getbounds(p) ? _polyder(p,der) : println("raise out of bounds error") # returns derivate polynomial that can be evaulated with polyval

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
    if lower(b1) ≥ upper(b2) && !(includeslower(b1) && includesupper(b2)) || lower(b1) > upper(b2)
        return 1 # b1 higher than b2
    end
    if lower(b2) ≥ upper(b1) && !(includeslower(b2) && includesupper(b1)) || lower(b2) > upper(b1)
        return -1 # b1 higher than b2
    end
    l = max(lower(b1), lower(b2))
    u = min(upper(b1), upper(b2))
    return Bound(l,u,includeLower = l in b1 && l in b2, includeUpper = u in b1 && u in b2)
end

lower(b::Bound) = b.lower
upper(b::Bound) = b.upper
includesupper(b::Bound) = upper(b) in b
includeslower(b::Bound) = lower(b) in b
unbounded() = Bound(-Inf,Inf)

#mutable struct CoefPoly <: AbstractPoly
#    poly::Poly
#    bounds::Bound
#end
#
#create_coef_poly(coefs) = CoefPoly(Poly(coefs), unbounded())
#create_coef_poly(coefs, b::Bound) = CoefPoly(Poly(coefs), b)
#create_coef_poly(coefs, lower, upper) = CoefPoly(Poly(coefs), Bound(lower, upper))
#createpoly(::CoefSpecifier, coefs) = create_coef_poly(coefs) # write a macro to do this?
#createpoly(::CoefSpecifier, coefs, b::Bound) = create_coef_poly(coefs, bound)
#createpoly(::CoefSpecifier, coefs, lower, upper) = create_coef_poly(coefs, lower, upper)
#_polyval(p::CoefPoly, x) = polyval(p.poly, x)
#_polyder(p::CoefPoly, der) = polyder(p.poly, der) # compute the derth derivate of p


#mutable struct GlobalPoly <: AbstractPoly
#    poly::Vector{AbstractPoly}#array of polynomials
#    bounds::Vector{Bound}
#    GlobalPoly(p::AbstractPoly, b::Bound) = new([p], [b]) # new initializing globalPoly must be created with an initial poly to avoid checks further down the line
#end
#
#create_global_poly(p::AbstractPoly, b::Bound) = GlobalPoly(p,b)
#
#function push!(gloPol::GlobalPoly, locPol::AbstractPoly, gloBound::Bound) # check to make sure bounds of lower are same or less restricted than global bounds
#    if !includeslower(getbounds(locPol)) #check first that local and global bound endpoints are compatible
#        if includeslower(gloBound)
#            return println("raise lower bound incompattible error")
#        end
#    elseif !includesupper(getbounds(locPol))
#        if includesupper(gloBound)
#            return println("raise upper bound incompattible error")
#        end
#    end
#    for i in 1:length(gloPol.bounds)
#        if !(typeof(intersect(gloBound, gloPol.bounds[i]))<:Number)
#            return println("raise bounds intersect with existing poly error")
#        end
#    end
#    push!(gloPol.poly, locPol)
#    push!(gloPol.bounds, gloBound)
#end
#push!(gloPol::GlobalPoly, locPol::AbstractPoly, gloLower, gloUpper) = push!(gloPol, locPol, Bound(gloLower, gloUpper)) 
#push!(gloPol::GlobalPoly, locPol::AbstractPoly) = push!(gloPol, locPol, getbounds(locPol))
#
#function getbounds(p::GlobalPoly)
#    low,up = promote(lower(p.bounds[1]), upper(p.bounds[end]))
#    return Bound(low, up, includeLower=includeslower(p.bounds[1]), includeUpper=includesupper(p.bounds[end]))
#end
#
#function _polyval(p::GlobalPoly, x)
#    for i in 1:length(p.bounds)
#        if x in p.bounds[i]
#            return _polyval(p.poly[i], global_to_local(getbounds(p.poly[i]), p.bounds[i], x))
#        end
#    end
#    return println("raise not found in set")
#end
#
## following Julia convert() convention where thing being converted to comes first
#local_to_global(gloB::Bound, locB::Bound, x::Number) = lower(gloB) + (upper(gloB) - lower(gloB)) * (x - lower(locB)) / ( upper(locB) - lower(locB) )
#global_to_local(locB::Bound, gloB::Bound, x::Number) = lower(locB) + (upper(locB) - lower(locB)) * (x - lower(gloB)) / ( upper(gloB) - lower(gloB) )


#mutable struct lagrangepoly <: abstractlagrangepoly
#    x::vector{number}
#    y::Vector{Number}
#    weights::Vector
#    bounds::Bound
#end
#
#function lagrange_bary_weights(x)
#    numPoints = length(x)
#    weights = [prod(1/(x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints]
#    return weights
#end
#
#function lagrange_eval_weights(p::AbstractLagrangePoly, xeval, y)
#    x = p.x
#    w = p.w
#    num = mapreduce((x,y,w)-> w*y / (xeval - x),+,x,y,w)
#    denom = mapreduce((x,w)-> w / (xeval - x),+,x,w)
#    val = num/denom
#    return isnan(val) ? y[x.==xeval][1] : val
#end

#create_lagrange_poly(x, y, b::Bound) = LagrangePoly(x, y, lagrange_bary_weights(x), b)
#create_lagrange_poly(x, y, lower, upper) = LagrangePoly(x, y, lagrange_bary_weights(x), Bound(lower, upper))
#createpoly(::LagrangeSpecifier,x,y,b::Bound) =  create_lagrange_poly(x,y,b::Bound)
#createpoly(::LagrangeSpecifier,x,y,lower, upper) =  create_lagrange_poly(x,y,lower, upper)
#updatepoly!(p::AbstractLagrangePoly, y) = p.y = y
#_polyval(p::AbstractLagrangePoly,x) = lagrange_eval_weights(p,x, p.y)
#_polyval(p::AbstractLagrangePoly,x,y) = lagrange_eval_weights(p, x, y)
#_polyval!(p::AbstractLagrangePoly,x,y) = begin updatepoly!(y); lagrange_eval_weights(p, x, y) end
#
#mutable struct LagrangeDerPoly <: AbstractLagrangePoly # derivate of 
#    x::Vector{Number}
#    y::Vector{Number}
#    #weights::Vector
#    der::Int # the order derivate to calculate
#    bounds::Bound
#end
#
#_polyder(p::AbstractLagrangePoly, der) = LagrangeDerPoly(p.x,p.y,der,p.bounds)
#
#function lagrange_der_weights(xval, x, der) # evaluate the derivate weight of the lagrange poly at xval
#        xlen = length(x)
#    if der == 1
#        return [sum( prod((xval - x[m])  / (x[j] - x[m]) for m = 1:xlen if m ≠ j && m ≠ l) /(x[j] - x[l]) for l = 1:xlen if l ≠ j) for j = 1:xlen] 
#    end
#end
#
#function lagrange_eval_weights(p::LagrangeDerPoly, xval, y)
#    return sum(y .* lagrange_der_weights(xval, p.x, p.der))
#end


#mutable struct LagrangePolyLite <: AbstractLagrangePoly # using less space as not including bounds
#    x::Vector{Number}
#    y::Vector{Number}
#    weights::Vector
#end
#create_lagrange_poly(x, y) = LagrangePolyLite(x, y, lagrange_bary_weights(x))
#createpoly(::LagrangeSpecifier,x,y) =  create_lagrange_poly(x,y)
#getbounds(p::LagrangePolyLite) = Bound(min(p.x...), max(p.x...))


#mutable struct LGRPoly <: AbstractLagrangePoly 
#    x::Vector{Number}
#    y::Vector{Number}
#    weights::Vector
#end
#function create_LGR_poly(y)
#    lgrpoints = lgr_points(length(y))
#    LGRPoly(lgrpoints,y,lagrange_bary_weights(lgrpoints))
#end
#createpoly(::LGRSpecifier,y) =  create_LGR_poly(y)
#getbounds(p::LGRPoly) = Bound(-1, 1, includeUpper = false)

# LGR points

#function lgr_points(order) # does not include endpoints
#    if order>1
#        N = 1:order-2
#        a = zeros(order-1)
#        b = zeros(order-2)
#        a[1] = 1 / 3 
#        a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
#        b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
#        J = SymTridiagonal(a,b)
#        return pushfirst!(eigvals(J),-1)
#    end
#end

# The Vandermonde matrix method

#macro genvandpoly(index,leny,x,y) # generate poly functions using the Vandermonde matrix method
#    leny = :($leny)
#    symX = []
#    symY = []
#    A = Array{Expr}(undef, leny, leny)
#    for i = 1:leny
#        push!(symX, Meta.parse("x$i"))
#        push!(symY, Meta.parse("y$i"))
#    end
#    for i = 1:leny # Assemble A matrix
#        for j = 1:leny
#            pwr = leny-j
#            var = symX[i]
#            A[i,j] = :($var^$pwr)
#        end
#    end
#    coefs = Algebra.:*(Algebra.inv(A),symY) # invert to find coefs of polynomial
#    for i = 1:leny
#        str = repr(coefs[i])
#        for j = 1:leny
#            str = replace(str,"y"*string(j)=>"$y"*"["*string(j)*"]")
#            str = replace(str,"x"*string(j)=>"$x"*"["*string(j)*"]")
#        end
#        coefs[i] = Meta.parse(chop(str, head=2, tail=1))
#    end
#
#    return coefs[:($index)]
#end
#
#function create_vand_poly(x,y, b::Bound)
#    xLength = length(x)
#    if xLength == 2
#        c1 = @genvandpoly(1, 2,:($x),:($y))
#        c2 = @genvandpoly(2, 2,:($x),:($y))
#        return create_coef_poly([c2,c1], b)
#    elseif xLength == 3
#        c1 = @genvandpoly(1, 3,:($x),:($y))
#        c2 = @genvandpoly(2, 3,:($x),:($y))
#        c3 = @genvandpoly(3, 3,:($x),:($y))
#        return create_coef_poly([c3,c2,c1], b)
#    elseif xLength == 4
#        c1 = @genvandpoly(1, 4,:($x),:($y))
#        c2 = @genvandpoly(2, 4,:($x),:($y))
#        c3 = @genvandpoly(3, 4,:($x),:($y))
#        c4 = @genvandpoly(4, 4,:($x),:($y))
#        return create_coef_poly([c4,c3,c2,c1], b)
#    elseif xLength == 5
#        c1 = @genvandpoly(1, 5,:($x),:($y))
#        c2 = @genvandpoly(2, 5,:($x),:($y))
#        c3 = @genvandpoly(3, 5,:($x),:($y))
#        c4 = @genvandpoly(4, 5,:($x),:($y))
#        c5 = @genvandpoly(5, 5,:($x),:($y))
#        return create_coef_poly([c5,c4,c3,c2,c1], b)
#    else 
#        create_coef_poly(inv([xval^j for xval in eachindex(x), j = 0:xlength - 1]) * y, b) # create poly coefficient by inverting array of linear equations
#        println("Warning, function behaves slow when exceeding more than 5 datapoints")
#    end
#end
#create_vand_poly(x,y) = create_vand_poly(x,y,Bound(min(x...), max(x...)))
#createpoly(::VandSpecifier,x,y,b::Bound) =  create_vand_poly(x,y, b::Bound)
#createpoly(::VandSpecifier,x,y) =  create_vand_poly(x,y)

# Dispatch createpoly functions

#function createpoly(y)
#    create_LGR_poly(y)
#end
#
#function createpoly(c, b::Bound)
#    create_coef_poly(c)
#end
#
#function createpoly(x,y)
#    if length(x) < 6
#        create_vand_poly(x,y)
#    else
#        create_lagrange_poly(x,y)
#    end
#end
#
#function createpoly(x,y,b::Bound)
#    if size(x) < 6
#        create_vand_poly(x,y,b)
#    else
#        create_lagrange_poly(x,y,b)
#    end
#end
