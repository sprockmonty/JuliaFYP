using Polynomials
abstract type AbstractPoly end

abstract type AbstractPolyMethod end
struct LagrangePolyMethod <: AbstractPolyMethod end


mutable struct CoefPoly{T} <: AbstractPoly
    poly::Poly{T}
    lowerBound::Number
    upperBound::Number # check that lower is below upper
end

mutable struct LagrangePoly{T} <: AbstractPoly
    poly::Poly{T}
    lowerBound::Number
    upperBound::Number # check that lower is below upper
    x::Vector{T} # am I doing this T thing right?
    y::Vector{T}
end


# lagrange polys
function lagrangeEval(x::Vector, y::Vector) # returns polynomaial with coefficients of lagrange evaluation
    numPoints = length(x)
    lagrangeVal = sum(y[i] * prod((Poly([0,1]) - x[j]) / (x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints)
end

create_legrange_poly(x, y) = LagrangePoly(lagrangeEval(x,y), min(x...), max(x...), x, y)
create_legrange_poly(x, y, lower::Real, upper::Real) = LagrangePoly(lagrangeEval(x,y), lower, upper)


### test code
function piecewise_poly(y::Function, time, method::LagrangePolyMethod)
    lagrangeEval(time, map(y, time))
end
t = [0,1,2,3,4,5]
z = piecewise_poly((x)->x^2, t, LagrangePolyMethod())

