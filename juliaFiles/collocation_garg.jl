using JuMP
using Ipopt
using Polynomials, MorePolynomials



struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array # incorperate lower bound and ub
end

mutable struct TrajProblem
    stateVectorGuess::AbstractPolynomial 
    controlVectorGuess::AbstractPolynomial
    stateVectorDiff::AbstractPolynomial
    diffMat::AbstractArray # differentiation matrix for state
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
    τ::AbstractVector
    numStates::Int   # can determine these from size of guess
    numControls::Int
    weights::AbstractVector
end

# work around to make current problem globally available
mutable struct CurrentProblem
    nullableProblem::Union{TrajProblem,Nothing}
end

const CURRENT_PROBLEM = CurrentProblem(nothing)
setCurrentProblem(problem::TrajProblem) = (CURRENT_PROBLEM.nullableProblem = problem)
getCurrentProblem() = CURRENT_PROBLEM.nullableProblem

function solve(problem::TrajProblem) # remove 2.5s to get this to work, and work out why it doesn't work with transformation
    setCurrentProblem(problem)
    model = Model(with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0)) # remove cpu time at some point
    stateLen = length(problem.stateVectorGuess) # number of points on each state vector
    controlLen = length(problem.controlVectorGuess)
    numStates = problem.numStates # number of state vectors
    numControls = problem.numControls
    # assign state and control varibles
    @variable(model, x[1:numStates, 1:stateLen])
    @variable(model, u[1:numControls, 1:controlLen])
    # set state and control variables to problem guess
    [set_start_value(x[i,j], problem.stateVectorGuess[j][2]) for i = 1:numStates, j = 1:stateLen] 
    [set_start_value(u[i,j], problem.controlVectorGuess[j][2]) for i = 1:numControls, j = 1:stateLen-1] 
    
    #register custom defined functions
    register(model, :stateDerUpdate!, numStates * stateLen + 1, stateDerUpdate!, ΔstateDerUpdate!)

    # state boundary constraints
    # initial state
    [fix(x[i,1], problem.boundaryConstraints.initialState[i]) for i = 1:length(problem.boundaryConstraints.initialState)]   # adjust this for lb and ub
    # final state
    [fix(x[i,end], problem.boundaryConstraints.finalState[i]) for i = 1:length(problem.boundaryConstraints.finalState)]   # adjust this for lb and ub

    @NLexpression(model, dynamics[i=1:stateLen-1], 2*x[i] + 2*u[i]*sqrt(x[i]))
    @NLexpression(model, objecFun[i=1:stateLen-1], x[i] + u[i]^2)

    for j = 1:stateLen - 1
        @NLconstraint(model, stateDerUpdate!(j,x...) - 2.5*dynamics[j] == 0)
    end
    @NLconstraint(model, x[1] + 2.5*sum(problem.weights[i] * dynamics[i] for i = 1:stateLen-1) - x[stateLen] == 0)

    @NLobjective(model, Min, 2.5*sum(problem.weights[i] * objecFun[i] for i = 1:stateLen-1)) 
    optimize!(model)
    return value.(x), value.(u)
end

function stateDerUpdate!(ζr, stateVector...) # make sure timestep vector matches length of state vector
    problem = getCurrentProblem()
    stateArray = zeros(typeof(problem.stateVectorGuess[1][1]), length(problem.stateVectorGuess))
    [stateArray[i] = stateVector[i] for i = 1:length(stateVector)]
    ζr = convert(Int, ζr)
    τ = problem.τ
    # assemble state and control vector from tuple inputs
    if stateArray == problem.stateVectorGuess.(τ)
        return problem.stateVectorDiff(τ[ζr])
    else
        problem.stateVectorGuess[:] = stateArray
        problem.stateVectorDiff[:] = problem.diffMat*stateArray
        return problem.stateVectorDiff(τ[ζr])
    end

end

function ΔstateDerUpdate!(g, ζr,stateVector...)
    problem = getCurrentProblem()
    stateArray = zeros(typeof(problem.stateVectorGuess[1][1]), length(problem.stateVectorGuess))
    [stateArray[i] = stateVector[i] for i = 1:length(stateVector)]
    ζr = convert(Int, ζr)
    τ = problem.τ
    g[1] = 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    if stateArray[:] == problem.stateVectorGuess.(τ)
        g[2:end] = problem.diffMat[ζr,:]
        return g
    else
        problem.stateVectorGuess[:] = stateArray
        problem.stateVectorDiff[:] = problem.diffMat*stateArray
        g[2:end] = problem.diffMat[ζr,:]
        return g
    end
end





####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################
numCollocationPoints = 39

lgrPoints = lgr_points(numCollocationPoints)
yploc = LagrangePoly([lgrPoints...,1],ones(numCollocationPoints+1))
uploc = LGRPoly(ones(numCollocationPoints))
boundaryConstraints = BoundaryConstraint([2],[1])
problem = TrajProblem(yploc, uploc, derivative(yploc),derivmatrix(yploc), boundaryConstraints,[lgrPoints...,1],1,1,lgr_weights(uploc))
x,u = solve(problem)


using Plots
plotly()
transform(x) = 2.5.*x .+ 2.5 
p1 = plot(transform([lgrPoints...,1]),y[1:end],labels=["" ""], ylabel="x")
p2 = plot(transform(lgrPoints),u[1:end],labels=["" ""], ylabel="u", xlabel="t")
plot(p1,p2,layout = (2,1))
