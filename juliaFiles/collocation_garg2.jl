using JuMP
using Ipopt
using FiniteDiff
using Polynomials
using MorePolynomials


struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array # incorperate lower bound and ub
end

mutable struct TrajProblem
    objectiveFunc::Function
    dynamicsFunc::Function
    pathConstraintFunc::Function
    stateVector::Array{AbstractPolynomial} 
    stateVectorDiff::Array{AbstractPolynomial}
    stateDiffMat::AbstractArray # differentiation matrix for state
    controlVector::Array{AbstractPolynomial}
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
    numStates::Int   # can determine these from size of guess
    numControls::Int
    weights::AbstractVector
    τf
    τp
    τo
    tf::Number
    t0::Number
end

# work around to make current problem globally available
mutable struct CurrentProblem
    nullableProblem::Union{TrajProblem,Nothing}
end

const CURRENT_PROBLEM = CurrentProblem(nothing)
setCurrentProblem(problem::TrajProblem) = (CURRENT_PROBLEM.nullableProblem = problem)
getCurrentProblem() = CURRENT_PROBLEM.nullableProblem

function solve(problem::TrajProblem) # multiple dispatch on this function
    setCurrentProblem(problem)
    model = Model(with_optimizer(Ipopt.Optimizer)) # remove cpu time at some point
    stateLen = length(problem.stateVector[1]) # number of points on each state vector
    controlLen = length(problem.controlVector[1])
    numStates = problem.numStates # number of state vectors
    numControls = problem.numControls

    # assign state and control varibles
    @variable(model, x[1:numStates, 1:stateLen])
    @variable(model, u[1:numControls, 1:controlLen])
    # set state and control variables to problem guess
    [set_start_value(x[i,j], problem.stateVector[i][j][2]) for i = 1:numStates, j = 1:stateLen] 
    [set_start_value(u[i,j], problem.controlVector[i][j][2]) for i = 1:numControls, j = 1:controlLen] 
    
    #register custom defined functions
    register(model, :collocateConstraint, (problem.numStates + problem.numControls) * stateLen + 2, collocateConstraint, ΔcollocateConstraint)
    register(model, :pathConstraint, (problem.numStates + problem.numControls) * stateLen + 1, pathConstraint, ΔpathConstraint)
    register(model, :objectiveFuncInterp, (problem.numStates + problem.numControls) * stateLen, objectiveFuncInterp, ΔobjectiveFuncInterp)
    register(model, :endPointInterp, (problem.numStates + problem.numControls) * stateLen + 1, endPointInterp, ΔendPointInterp)

    # state boundary constraints
    # initial state
    [fix(x[i,1], problem.boundaryConstraints.initialState[i]) for i = 1:length(problem.boundaryConstraints.initialState)]   # adjust this for lb and ub
    # final state
    [fix(x[i,end], problem.boundaryConstraints.finalState[i]) for i = 1:length(problem.boundaryConstraints.finalState)]   # adjust this for lb and ub

    for i = 1:length(problem.pathConstraintFunc(problem.stateVector, problem.controlVector, problem.τp[1])) * length(problem.τp)
        @NLconstraint(model, pathConstraint(i,x...,u...) <=0)
    end

    for i = 1:problem.numStates
        @NLconstraint(model, endPointInterp(i,x...,u...) ==0)
        for j = 1:length(problem.τf)
            @NLconstraint(model, collocateConstraint(i,j,x...,u...) == 0)
        end
    end
    
    @NLobjective(model, Min, objectiveFuncInterp(x...,u...)) 
    optimize!(model)
    return value.(x), value.(u), model
end

function endPointInterp(ζr, stateControlVector...)  
    ζr = convert(Int, ζr)
    problem = getCurrentProblem()
    τf = problem.τf 
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.stateVector[i][:] = stateVector[i,:]
    end
    for i = 1:problem.numControls
        problem.controlVector[i][:] = controlVector[i,:]
    end
    # return interpolated integral
    f = map(t -> problem.dynamicsFunc(problem.stateVector, problem.controlVector, t), τf)
    e = zeros(1, problem.numStates)
    weights = lgr_weights(length(τf))
    for i = 1:problem.numStates
        e[i] = problem.stateVector[i](1) - problem.stateVector[i](-1) -  (tf-t0) * 0.5 * sum([weights[k]*f[k][i] for k = 1:length(τf)])
    end
    return e[ζr]
end

function endPointInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    τf = problem.τf 
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    s = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numStates)
    c = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numControls)
    for i = 1:problem.numStates
        s[i] = LagrangePoly(problem.stateVector[i].x,convert(Array{Float64,1}, stateVector[i,:]))
    end
    for i = 1:problem.numControls
        c[i] = LagrangePoly(problem.controlVector[i].x,convert(Array{Float64,1}, controlVector[i,:]))
    end
    f = map(t -> problem.dynamicsFunc(s, c, t), τf)
    e = zeros(1, problem.numStates)
    weights = lgr_weights(length(τf))
    for i = 1:problem.numStates
        e[i] = s[i](1) - s[i](-1) -  (tf-t0) * 0.5 * sum([weights[k]*f[k][i] for k = 1:length(τf)])
    end
    return e
end

function ΔendPointInterp(g, ζr, stateControlVector...)
    problem = getCurrentProblem()
    g[1] = 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = FiniteDiff.finite_difference_jacobian(endPointInterp, [stateVector; controlVector])
    g[2:end] = getStateControlFromXU(problem, jacobian[ζr,:])
    return g
end

function objectiveFuncInterp(stateControlVector...)  
    problem = getCurrentProblem()
    τo = problem.τo 
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.stateVector[i][:] = stateVector[i,:]
    end
    for i = 1:problem.numControls
        problem.controlVector[i][:] = controlVector[i,:]
    end
    # return interpolated integral
    obj = map(t -> problem.objectiveFunc(problem.stateVector, problem.controlVector, t), τo)
    return (tf-t0) * 0.5 * sum(lgr_weights(length(τo)).*obj)
end

function objectiveFuncInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    τo = problem.τo 
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    s = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numStates)
    c = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numControls)
    for i = 1:problem.numStates
        s[i] = LagrangePoly(problem.stateVector[i].x,convert(Array{Float64,1}, stateVector[i,:]))
    end
    for i = 1:problem.numControls
        c[i] = LagrangePoly(problem.controlVector[i].x,convert(Array{Float64,1}, controlVector[i,:]))
    end
    obj = map(t -> problem.objectiveFunc(s, c, t), problem.τo )
    out = (tf-t0) * 0.5 * sum(lgr_weights(length(τo)).*obj)
    return out
end

function ΔobjectiveFuncInterp(g, stateControlVector...)
    problem = getCurrentProblem()
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = FiniteDiff.finite_difference_jacobian(objectiveFuncInterp, [stateVector; controlVector])
    g[:] = getStateControlFromXU(problem, jacobian)
    return g
end

function collocateConstraint(ζr, ζc, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    τf = problem.τf
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.stateVector[i][:] = stateVector[i,:]
        problem.stateVectorDiff[i][:] = problem.stateDiffMat[i]*problem.stateVector[i][:]
    end
    for i = 1:problem.numControls
        problem.controlVector[i][:] = controlVector[i,:]
    end
    ζ = zeros(problem.numStates, length(τf))
    f = map(t -> problem.dynamicsFunc(problem.stateVector, problem.controlVector, t), τf)
    for i = 1:problem.numStates
        for j = 1:length(τf)
            ζ[i,j] = problem.stateVectorDiff[i](τf[j]) - (tf-t0) * 0.5 * f[j][i]
        end
    end
    return ζ[ζr, ζc]
end

function collocateConstraint(stateControlVector::Array) # make sure timestep vector matches length of state vector
    problem = getCurrentProblem()
    problem = getCurrentProblem()
    τf = problem.τf 
    tf = problem.tf
    t0 = problem.t0
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    s = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numStates)
    sDiff = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numStates)
    c = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numControls)
    for i = 1:problem.numStates
        s[i] = LagrangePoly(problem.stateVector[i].x,convert(Array{Float64,1}, stateVector[i,:]))
        sDiff[i] = LagrangePoly(problem.stateVector[i].x,convert(Array{Float64,1}, problem.stateDiffMat[i]*stateVector[i,:]))
    end
    for i = 1:problem.numControls
        c[i] = LagrangePoly(problem.controlVector[i].x,convert(Array{Float64,1}, controlVector[i,:]))
    end
    ζ = zeros(problem.numStates, length(τf))
    f = map(t -> problem.dynamicsFunc(s, c, t), τf)
    for i = 1:problem.numStates
        for j = 1:length(τf)
            ζ[i,j] = sDiff[i](τf[j]) - (tf-t0) * 0.5 * f[j][i]
        end
    end
    return ζ
end
#function trueJacobian(stateVector, controlVector, τf)
#    problem = getCurrentProblem()
#    jac = zeros(length(τf), length(stateVector[1]) * 2)
#    for i = 1:problem.numStates*length(τf)
#        for j = 1:length(stateVector[1]) +length(controlVector[1])
#            if mod(j,2) == 1
#                jac[i,j] += problem.stateDiffMat[1][i,convert(Int,(j+1)/2)]
#            end
#        end
#    end
#    jac = jac + FiniteDiff.finite_difference_jacobian(collocateConstraint, [stateVector; controlVector])
#end
#function dynamicsFunc2(stateControlVector) 
#    problem = getCurrentProblem()
#    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
#    x1dot = 2 .* stateVector .+ 2 .* controlVector .* sqrt.(abs.(stateVector)) 
#    return x1dot
#end

function ΔcollocateConstraint(g, ζr, ζc, stateControlVector...)
    problem = getCurrentProblem()
    g[1:2] .= 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = FiniteDiff.finite_difference_jacobian(collocateConstraint, [stateVector; controlVector])
    jacobianRow = (ζc - 1) * problem.numStates + ζr # get the row of the jacobian that matches our output function
    g[3:end] = getStateControlFromXU(problem, jacobian[jacobianRow,:])
    return g
end

function pathConstraint(ζr, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    problem = getCurrentProblem()
    τp = problem.τp 
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.stateVector[i][:] = stateVector[i,:]
    end
    for i = 1:problem.numControls
        problem.controlVector[i][:] = controlVector[i,:]
    end

    notFlatPaths = map(t -> problem.pathConstraintFunc(problem.stateVector, problem.controlVector, t), τp) # paths before array flattening
    ζ = notFlatPaths[1]
    for i = 2:length(notFlatPaths)
        push!(ζ, notFlatPaths[i]...)
    end
    return ζ[ζr]
end

function pathConstraint(stateControlVector::Array) # make sure timestep vector matches length of state vector
    problem = getCurrentProblem()
    τp = problem.τp 
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    s = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numStates)
    c = Vector{Union{Nothing, AbstractPolynomial}}(nothing, problem.numControls)
    for i = 1:problem.numStates
        s[i] = LagrangePoly(problem.stateVector[i].x,convert(Array{Float64,1}, stateVector[i,:]))
    end
    for i = 1:problem.numControls
        c[i] = LagrangePoly(problem.controlVector[i].x,convert(Array{Float64,1}, controlVector[i,:]))
    end

    notFlatPaths = map(t -> problem.pathConstraintFunc(s, c, t), τp) # paths before array flattening
    ζ = notFlatPaths[1]
    for i = 2:length(notFlatPaths)
        push!(ζ, notFlatPaths[i]...)
    end
    return ζ
end

function ΔpathConstraint(g, ζr,stateControlVector...)
    problem = getCurrentProblem()
    g[1] = 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = FiniteDiff.finite_difference_jacobian(pathConstraint, [stateVector; controlVector])
    g[2:end] = getStateControlFromXU(problem, jacobian[ζr,:])
    return g
end

function getXUFromStateControl(problem, stateControlVector::Tuple) # get separate state and control vector matricies from input tuple, add type to problem
    stateControlVector = collect(stateControlVector)
    stateVector = zeros(problem.numStates, length(problem.stateVector[1]))
    controlVector = zeros(problem.numControls, length(problem.stateVector[1]))
    for i = 1:length(stateControlVector) ÷ (problem.numStates + problem.numControls)
        jx = (i-1) * (problem.numStates) + 1
        ju = (i-1) * (problem.numControls) + length(problem.stateVector[1]) * problem.numStates + 1
        stateVector[:,i] = stateControlVector[jx:jx+problem.numStates-1]
        controlVector[:,i] = stateControlVector[ju:ju+problem.numControls-1]
    end
    return stateVector, controlVector
end

function getXUFromStateControl(problem, stateControlVector::Array) # get separate state and control vector matricies from input tuple
    stateVector = Array{Union{Missing, typeof(stateControlVector[1]) }}(missing, problem.numStates, length(problem.stateVector[1]))
    controlVector = Array{Union{Missing,typeof(stateControlVector[1]) }}(missing, problem.numControls, length(problem.stateVector[1]))
    stateVector[:,:] = stateControlVector[1:problem.numStates, :]
    controlVector[:,:] = stateControlVector[problem.numStates + 1, :]
    return stateVector, controlVector
end

function getStateControlFromXU(problem, stateControlVector::Array) # get tuple of statecontrols from a 1D array with states and controls in a column row pattern
    x = ones(1,length(problem.stateVector[1])*problem.numStates)
    u = ones(1,length(problem.stateVector[1])*problem.numControls)
    for i = 1:length(problem.stateVector[1])
        j = (i-1)*(problem.numStates + problem.numControls) + 1
        x[(i-1)*problem.numStates+1:i*problem.numStates] = stateControlVector[j:j+problem.numStates-1]
        u[(i-1)*problem.numControls+1:i*problem.numControls] = stateControlVector[j+problem.numStates:j+problem.numStates+problem.numControls-1]
    end
    return [x u]
end

####################################################################################################################################
### problem specific functions #####################################################################################################
####################################################################################################################################

function dynamicsFunc(stateVector, controlVector, t) 
    x1dot = 2*stateVector[1](t) + 2 * controlVector[1](t) * sqrt(abs(stateVector[1](t))) 
    return [x1dot]
end
objectiveFunc(stateVector, controlVector, t) = 0.5 * stateVector[1](t)+controlVector[1](t)^2

function pathConstraintFunc(stateVector, controlVector, t) # less than or equal to output value. Output must by 1 D vector of constraints
    bound = zeros(4)
    bound[1] = - stateVector[1](t)
    bound[2] = stateVector[1](t) - 2
    bound[3] = -controlVector[1](t) - 4
    bound[3] = controlVector[1](t) - 2
    return bound
end

####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################

numCP = 15
T = 5
lgrPoints = lgr_points(numCP)
controlVector = LagrangePoly([lgrPoints...,1],ones(numCP+1))
timeStep = ones(1,numCP) * T/(numCP)
time = pushfirst!(cumsum(timeStep';dims=1)[:,1],0.0)

stateVector1 = LagrangePoly([lgrPoints...,1],ones(numCP+1))
boundaryConstraints = BoundaryConstraint([2],[1])
problem = TrajProblem(objectiveFunc, dynamicsFunc, pathConstraintFunc, [stateVector1], [derivative(stateVector1)], [derivmatrix(stateVector1)], [controlVector], boundaryConstraints, 1, 1, lgr_weights(numCP),lgr_points(numCP+2), lgrPoints, lgrPoints, T,0)
x,u,model = solve(problem)



using Plots
plotly()
p1 = plot(time, x[1,:],labels=["" ""], ylabel="Position q")
p2 = plot(time, x[2,:],labels=["" ""],ylabel="Angle θ")
p3 = plot(time, u',labels=["" ""],ylabel="Force", xlabel="time t")
plot(p1,p2,p3,layout = (3,1), )


### temp
f = [78.630478282274, -10.746997196711233, -29.59005409998428, 23.958990876433276, 3.3364963223970188]
x = [2.0,0.31100483259184986,1.6478987634606341,7.357266844312955,1.0]
stateVec = LagrangePoly()
u = [26.38585863829565 -10.193162885278168 -12.80895893827755 1.7040925032356469 0.6682481611985094]
end = [-1.4325349678164079] 
