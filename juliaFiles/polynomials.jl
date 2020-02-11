struct τ{T<:Array{Any}} # parametric types
    x::T
    y::T
    function τ{T}(x::T, y::T) where T<:Array{Any}
        @assert length(x) == length(y)  
            new(x,y)
    end
end
τ([2,3],[2,3])

# unit tests
# Using function e^x
# test exact evaluation at points

using Test
funcX = [-1.0 0.0 1.0 2.0]
yLagrange = lagrangeEval(funcX, map(exp, funcX))
@test yLagrange(0) == 1
@test yLagrange(0.5) ≈ 1.606725057


# Function for evaluating lagrange polynomial basis, returns function handle with single input for point to evaluate
function lagrangeEval(points::Array, Y::Array) # define this as a type?
    function lagrangePoly(X::Real)
        numPoints = length(points)
        lagrangeVal = sum(Y[i] * prod((X - points[j]) / (points[i] - points[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints)
        return lagrangeVal
    end
    return lagrangePoly
end


