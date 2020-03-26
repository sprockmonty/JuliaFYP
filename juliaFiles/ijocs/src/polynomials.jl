using LinearAlgebra
using Polynomials
abstract type AbstractPoly end

abstract type AbstractPolyMethod end


mutable struct CoefPoly{T} <: AbstractPoly
    poly::Poly{T}
    lowerBound::Number
    upperBound::Number # check that lower is below upper
end

mutable struct LagrangePoly <: AbstractPoly
    lowerBound::Number
    upperBound::Number # check that lower is below upper
    x::Vector # am I doing this T thing right?
    y::Vector
    weights::Vector
end


# lagrange polys
#function lagrangeEval(x::Vector, y::Vector) # returns polynomaial with coefficients of lagrange evaluation
#    numPoints = length(x)
#    lagrangeVal = sum(y[i] * prod((Poly([0,1]) - x[j]) / (x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints)
#end
function lagrangeBaryWeights(x)
    numPoints = length(x)
    weights = [prod(1/(x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints]
    return weights
end

function lagrangeEvalWeights(xeval, x, y, w)
    num = mapreduce((x,y,w)-> w*y / (xeval - x),+,x,y,w)
    denom = mapreduce((x,w)-> w / (xeval - x),+,x,w)
    val = num/denom
    return isnan(val) ? y[x.==xeval] : val
end

create_lagrange_poly(x, y) = LagrangePoly(min(x...), max(x...), x, y, lagrangeBaryWeights(x))
create_lagrange_poly(x, y, lower::Real, upper::Real) = LagrangePoly(lower, upper, x, y, lagrangeBaryWeights(x))
polyval(p::LagrangePoly,x::Number) = p.lowerBound≤x≤p.upperBound ? lagrangeEvalWeights(x, p.x, p.y, p.weights) : println("raise and error for out of bounds here")


# LGR points

function lgrPoints(order) # does not include endpoints
    if order>2
        N = 1:order-2
        a = zeros(order-1)
        b = zeros(order-2)
        a[1] = 1 / 3 
        a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
        b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
        J = SymTridiagonal(a,b)
        return eigvals(J)
    end
end




### test code
function piecewise_poly(y::Function, time, method::LagrangePolyMethod)
    lagrangeEval(time, map(y, time))
end
t = [0,1,2,3,4,5]
z = piecewise_poly((x)->x^2, t, LagrangePolyMethod())

using Plots
plotly()
scatter(lgrPoints(5),[0,0,0,0])
