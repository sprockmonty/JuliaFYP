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
    sVector::Array{AbstractPolynomial} 
    stateVectorDiff::Array{AbstractPolynomial}
    sVectorDiff::Array{AbstractPolynomial}
    stateDiffMat::AbstractArray # differentiation matrix for state
    controlVector::Array{AbstractPolynomial}
    cVector::Array{AbstractPolynomial} 
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
    numStates::Int   # can determine these from size of guess
    numControls::Int
    τf
    τfWeights
    τp
    τpWeights
    τo
    τoWeights
    t0::Number
    tf::Number
    endJac::Array
    colJac::Array
    objJac::Array
    patJac::Array
    stateControlCache
    endJacCache
    colJacCache
    objJacCache
    patJacCache
    function TrajProblem(objectiveFunc, dynamicsFunc, pathConstraintFunc, stateVector, controlVector, boundaryConstraints, τf, τp, τo, t0,tf)
        numStates = length(stateVector)
        numControls = length(controlVector)
        sVector = [LagrangePoly(stateVector[i].x, stateVector[i][:]) for i = 1:numStates]
        cVector = [LagrangePoly(controlVector[i].x, controlVector[i][:]) for i = 1:numControls]
        stateVectorDiff = [derivative(stateVector[i]) for i = 1:numStates]
        sVectorDiff = [derivative(stateVector[i]) for i = 1:numStates]
        stateDiffMat = [derivmatrix(stateVector[i]) for i = 1:numStates]
        τfWeights = lgr_weights(length(τf))
        τpWeights = lgr_weights(length(τp))
        τoWeights = lgr_weights(length(τo))
        stateControlLen = numStates*length(stateVector[1]) + numControls * length(controlVector[1])
        endJac = zeros(numStates, stateControlLen)
        colJac = zeros(length(τf)*numStates, stateControlLen)
        objJac = zeros(stateControlLen)
        patJac = zeros(length(pathConstraintFunc(stateVector, controlVector, τp[1])) * length(τp), stateControlLen) 
        stateControlCache = tuple(zeros(numStates*length(stateVector[1]) + numControls*length(controlVector[1]))...)
        endJacCache = tuple(zeros(numStates*length(stateVector[1]) + numControls*length(controlVector[1]))...)
        colJacCache = tuple(zeros(numStates*length(stateVector[1]) + numControls*length(controlVector[1]))...)
        objJacCache = tuple(zeros(numStates*length(stateVector[1]) + numControls*length(controlVector[1]))...)
        patJacCache = tuple(zeros(numStates*length(stateVector[1]) + numControls*length(controlVector[1]))...)
        return new(objectiveFunc, dynamicsFunc, pathConstraintFunc, stateVector, sVector, stateVectorDiff, sVectorDiff, stateDiffMat, controlVector, cVector, boundaryConstraints, numStates, numControls, τf, τfWeights, τp, τpWeights, τo, τoWeights, t0, tf, endJac, colJac, objJac, patJac, stateControlCache, endJacCache,colJacCache, objJacCache, patJacCache)
    end
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
    updateStateControl(problem, stateControlVector)
    # return interpolated integral
    f = map(t -> problem.dynamicsFunc(problem.stateVector, problem.controlVector, t), τf)
    e = zeros(1, problem.numStates)
    weights = lgr_weights(length(τf))
    for i = 1:problem.numStates
        e[i] = problem.stateVector[i](1) - problem.stateVector[i](-1) -  (tf-t0) * 0.5 * sum([problem.τfWeights[k]*f[k][i] for k = 1:length(τf)])
    end
    return e[ζr]
end

function endPointInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    τf = problem.τf 
    tf = problem.tf
    t0 = problem.t0
    updateSC(problem, stateControlVector)
    s = problem.sVector
    c = problem.cVector
    f = map(t -> problem.dynamicsFunc(s, c, t), τf)
    e = zeros(1, problem.numStates)
    for i = 1:problem.numStates
        e[i] = s[i](1) - s[i](-1) -  (tf-t0) * 0.5 * sum([problem.τfWeights[k]*f[k][i] for k = 1:length(τf)])
    end
    return e
end

function ΔendPointInterp(g, ζr, stateControlVector...)
    problem = getCurrentProblem()
    g[1] = 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    if stateControlVector!== problem.endJacCache
        stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
        problem.endJac = FiniteDiff.finite_difference_jacobian(endPointInterp, [stateVector; controlVector])
        problem.endJacCache = stateControlVector
    end
    g[2:end] = getStateControlFromXU(problem, problem.endJac[ζr,:])
    return g
end

function objectiveFuncInterp(stateControlVector...)  
    problem = getCurrentProblem()
    τo = problem.τo 
    tf = problem.tf
    t0 = problem.t0
    updateStateControl(problem, stateControlVector)
    # return interpolated integral
    obj = map(t -> problem.objectiveFunc(problem.stateVector, problem.controlVector, t), τo)
    return (tf-t0) * 0.5 * sum(problem.τoWeights.*obj)
end

function objectiveFuncInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    τo = problem.τo 
    tf = problem.tf
    t0 = problem.t0
    updateSC(problem, stateControlVector)
    s = problem.sVector
    c = problem.cVector
    obj = map(t -> problem.objectiveFunc(s, c, t), problem.τo )
    return (tf-t0) * 0.5 * sum(problem.τoWeights.*obj)
end

function ΔobjectiveFuncInterp(g, stateControlVector...)
    problem = getCurrentProblem()
    if stateControlVector!== problem.objJacCache
        stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
        problem.objJac = FiniteDiff.finite_difference_jacobian(objectiveFuncInterp, [stateVector; controlVector])
        problem.objJacCache = stateControlVector
    end
    g[:] = getStateControlFromXU(problem, problem.objJac)
    return g
end

function collocateConstraint(ζr, ζc, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    τf = problem.τf
    tf = problem.tf
    t0 = problem.t0
    updateStateControl(problem, stateControlVector)
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
    updateSC(problem, stateControlVector)
    s = problem.sVector
    c = problem.cVector
    sDiff = problem.sVectorDiff
    ζ = zeros(problem.numStates, length(τf))
    f = map(t -> problem.dynamicsFunc(s, c, t), τf)
    for i = 1:problem.numStates
        for j = 1:length(τf)
            ζ[i,j] = sDiff[i](τf[j]) - (tf-t0) * 0.5 * f[j][i]
        end
    end
    return ζ
end

function ΔcollocateConstraint(g, ζr, ζc, stateControlVector...)
    problem = getCurrentProblem()
    g[1:2] .= 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    if stateControlVector!== problem.colJacCache
        stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
        problem.colJac = FiniteDiff.finite_difference_jacobian(collocateConstraint, [stateVector; controlVector])
        problem.colJacCache = stateControlVector
    end
    jacobianRow = (ζc - 1) * problem.numStates + ζr # get the row of the jacobian that matches our output function
    g[3:end] = getStateControlFromXU(problem, problem.colJac[jacobianRow,:])
    return g
end

function pathConstraint(ζr, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    problem = getCurrentProblem()
    τp = problem.τp 
    updateStateControl(problem, stateControlVector)
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
    updateSC(problem, stateControlVector)
    s = problem.sVector
    c = problem.cVector

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
    if stateControlVector!== problem.patJacCache
        stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
        problem.patJac = FiniteDiff.finite_difference_jacobian(pathConstraint, [stateVector; controlVector])
        problem.patJacCache = stateControlVector
    end
    g[2:end] = getStateControlFromXU(problem, problem.patJac[ζr,:])
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

function updateStateControl(problem, stateControlVector)
    stateControlVector == problem.stateControlCache && return 1
    problem.stateControlCache = stateControlVector
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.stateVector[i][:] = stateVector[i,:]
        problem.stateVectorDiff[i][:] = problem.stateDiffMat[i]*problem.stateVector[i][:]
    end
    for i = 1:problem.numControls
        problem.controlVector[i][:] = controlVector[i,:]
    end
end

function updateSC(problem, stateControlVector)
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    for i = 1:problem.numStates
        problem.sVector[i][:] = convert(Array{Float64,1}, stateVector[i,:])
        problem.sVectorDiff[i][:] = convert(Array{Float64,1}, problem.stateDiffMat[i]*stateVector[i,:])
    end
    for i = 1:problem.numControls
        problem.cVector[i][:] = convert(Array{Float64,1}, controlVector[i,:])
    end
end
####################################################################################################################################
### problem specific functions #####################################################################################################
####################################################################################################################################

function dynamicsFunc(stateVector, controlVector, t) 
    l = 0.5
    m1 = 1
    m2 = 0.3
    g = 9.81
    x2 = stateVector[2](t)
    u = controlVector[1](t)
    x1dot = stateVector[3](t)  
    x2dot = stateVector[4](t)
    x3dot = (l * m2 * sin(x2)*x2dot^2 + u + m2 * g * cos(x2)*sin(x2) )/ (m1 + m2 *(1 - cos(x2)^2))
    x4dot = - (l*m2*cos(x2)*sin(x2)*x2dot^2+u*cos(x2)+(m1+m2)*g*sin(x2)) / (l*m1+l*m2*(1-cos(x2)^2))
    return [x1dot;x2dot;x3dot;x4dot]
end
objectiveFunc(stateVector, controlVector, t) = controlVector[1](t)^2

function pathConstraintFunc(stateVector, controlVector, t) # less than or equal to output value. Output must by 1 D vector of constraints
    dmax = 2
    umax = 20
    bounds = zeros(4)
    bounds[1] = -dmax .- stateVector[1](t)
    bounds[2] = -dmax .+ stateVector[1](t)
    bounds[3] = -umax .- controlVector[1](t)
    bounds[4] = -umax .+ controlVector[1](t)
    return bounds
end

####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################

numCP = 20
T = 2
lgrPoints = lgr_points(numCP)
controlVector = LagrangePoly([lgrPoints...,1],zeros(numCP+1))
timeStep = ones(1,numCP) * T/(numCP)
time = pushfirst!(cumsum(timeStep';dims=1)[:,1],0.0)

stateVector1 = LagrangePoly([lgrPoints...,1],ones(numCP+1) .* time ./ T)
stateVector2 = LagrangePoly([lgrPoints...,1],ones(numCP+1) .* π .* time ./ T)
stateVector3 = LagrangePoly([lgrPoints...,1],zeros(numCP+1))
stateVector4 = LagrangePoly([lgrPoints...,1],zeros(numCP+1))
boundaryConstraints = BoundaryConstraint([0;0;0;0],[1;π;0;0])
τf = lgr_points(numCP)
τp = lgr_points(numCP)
τo = lgr_points(numCP)
problem = TrajProblem(objectiveFunc, dynamicsFunc, pathConstraintFunc, [stateVector1, stateVector2, stateVector3, stateVector4], [controlVector], boundaryConstraints,τf, τp, τo, 0, T)
x,u,model = solve(problem)



using Plots
plotly()
p1 = plot(problem.stateVector[1],labels=["" ""], ylabel="Position q")
p2 = plot(problem.stateVector[2],labels=["" ""],ylabel="Angle θ")
p3 = plot(problem.controlVector[1],labels=["" ""],ylabel="Force", xlabel="time t")
plot(p1,p2,p3,layout = (3,1), )


