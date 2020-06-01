using JuMP
using Ipopt
using ForwardDiff


struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array # incorperate lower bound and ub
end

mutable struct TrajProblem
    objectiveFunc::Function
    dynamicsFunc::Function
    pathConstraintFunc::Function
    controlVectorGuess::Array
    stateVectorGuess::Array # each row represents a state, maybe use an add state guess function which adds each state guess
    timeStep::Union{Array, Float64} # if only one value input, turn this into an array of that one timestep, must be one less than state length
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
    numStates::Int   # can determine these from size of guess
    numControls::Int
    numCollocationPoints::Int # user / program specified
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
    model = Model(with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0)) # remove cpu time at some point

    # assign state and control varibles
    @variable(model, x[1:problem.numStates, 1:problem.numCollocationPoints])
    @variable(model, u[1:problem.numControls, 1:problem.numCollocationPoints])
    @constraint(model,x[1,:] >= 1.0)
    @constraint(model,x[3,:] >= 0.6)
    # set state and control variables to problem guess
    [set_start_value(x[i,j], problem.stateVectorGuess[i,j]) for i = 1:problem.numStates, j = 1:problem.numCollocationPoints] 
    [set_start_value(u[i,j], problem.controlVectorGuess[i,j]) for i = 1:problem.numControls, j = 1:problem.numCollocationPoints] 
    
    #register custom defined functions
    register(model, :collocateConstraint, (problem.numStates + problem.numControls) * problem.numCollocationPoints + 2, collocateConstraint, ΔcollocateConstraint)
    register(model, :pathConstraint, (problem.numStates + problem.numControls) * problem.numCollocationPoints + 1, pathConstraint, ΔpathConstraint)
    register(model, :objectiveFuncInterp, (problem.numStates + problem.numControls) * problem.numCollocationPoints, objectiveFuncInterp, ΔobjectiveFuncInterp)

    # state boundary constraints
    # initial state
    [fix(x[i,1], problem.boundaryConstraints.initialState[i]) for i = 1:length(problem.boundaryConstraints.initialState)]   # adjust this for lb and ub
    
    for i = 1:length(problem.boundaryConstraints.finalState)
        if !ismissing(problem.boundaryConstraints.finalState[i])
            fix(x[i,end], problem.boundaryConstraints.finalState[i])
        end
    end

    for i = 1:length(problem.pathConstraintFunc(problem.stateVectorGuess, problem.controlVectorGuess))
        @NLconstraint(model, pathConstraint(i,x...,u...) <=0)
    end

    for i = 1:problem.numStates
        for j = 1:problem.numCollocationPoints - 1
            @NLconstraint(model, collocateConstraint(i,j,x...,u...) == 0)
        end
    end
    
    @NLobjective(model, Min, objectiveFuncInterp(x...,u...)) 
    optimize!(model)
    return value.(x), value.(u)
end

function objectiveFuncInterp(stateControlVector...)  
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)

    # return interpolated integral
    return sum(0.5 .* problem.timeStep[1] .* (problem.objectiveFunc(stateVector[2:end,:], controlVector[2:end]) .+ problem.objectiveFunc(stateVector[2:end,:], controlVector[2:end])) ) # should we use a map function here? test with a time, also needs to be rewritten for control vector with multiple dimensions
end

function objectiveFuncInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    # return interpolated integral
    return [sum(0.5 .* problem.timeStep[1] .* (problem.objectiveFunc(stateVector[2:end,:], controlVector[2:end]) .+ problem.objectiveFunc(stateVector[1:end-1,:], controlVector[1:end-1])))]  # should we use a map function here? test with a time, also needs to be rewritten for control vector with multiple dimensions
end

function ΔobjectiveFuncInterp(g, stateControlVector...)
    problem = getCurrentProblem()
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = ForwardDiff.jacobian(objectiveFuncInterp, [stateVector; controlVector])
    g[:] = getStateControlFromXU(problem, jacobian)
    return g
end

function collocateConstraint(ζr, ζc, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64

    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)

    ΔstateVector = 0.5 .* problem.timeStep[1] .* (problem.dynamicsFunc(stateVector[:,2:end], controlVector[:,2:end]) + problem.dynamicsFunc(stateVector[:,1:end-1], controlVector[:,1:end-1])) # how do we ensure dynamicsFunc has the right inputs and outputs if it is user defined? also have a look at making timestep dynamic for each value
    ζ = stateVector[:, 2:end] - stateVector[:,1:end-1] - ΔstateVector
    return ζ[ζr, ζc]
end

function collocateConstraint(stateControlVector::Array) # make sure timestep vector matches length of state vector
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)

    ΔstateVector = 0.5 .* problem.timeStep[1] .* (problem.dynamicsFunc(stateVector[:,2:end], controlVector[:,2:end]) + problem.dynamicsFunc(stateVector[:,1:end-1], controlVector[:,1:end-1])) # how do we ensure dynamicsFunc has the right inputs and outputs if it is user defined? also have a look at making timestep dynamic for each value
    ζ = stateVector[:, 2:end] - stateVector[:,1:end-1] - ΔstateVector
    return ζ
end

function ΔcollocateConstraint(g, ζr, ζc, stateControlVector...)
    problem = getCurrentProblem()
    g[1:2] .= 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = ForwardDiff.jacobian(collocateConstraint, [stateVector; controlVector])
    jacobianRow = (ζc - 1) * problem.numStates + ζr # get the row of the jacobian that matches our output function
    g[3:end] = getStateControlFromXU(problem, jacobian[jacobianRow,:])
    return g
end

function pathConstraint(ζr, stateControlVector...) # make sure timestep vector matches length of state vector
    ζr = convert(Int, ζr)
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)

    ζ = problem.pathConstraintFunc(stateVector, controlVector)
    return ζ[ζr]
end

function pathConstraint(stateControlVector::Array) # make sure timestep vector matches length of state vector
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    ζ = problem.pathConstraintFunc(stateVector, controlVector)
    return ζ
end

function ΔpathConstraint(g, ζr,stateControlVector...)
    problem = getCurrentProblem()
    g[1] = 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = ForwardDiff.jacobian(pathConstraint, [stateVector; controlVector])
    g[2:end] = getStateControlFromXU(problem, jacobian[ζr,:])
    return g
end

function getXUFromStateControl(problem, stateControlVector::Tuple) # get separate state and control vector matricies from input tuple, add type to problem
    stateControlVector = collect(stateControlVector)
    stateVector = zeros(problem.numStates, problem.numCollocationPoints)
    controlVector = zeros(problem.numControls, problem.numCollocationPoints)
    for i = 1:length(stateControlVector) ÷ (problem.numStates + problem.numControls)
        jx = (i-1) * (problem.numStates) + 1
        ju = (i-1) * (problem.numControls) + problem.numCollocationPoints * problem.numStates + 1
        stateVector[:,i] = stateControlVector[jx:jx+problem.numStates-1]
        controlVector[:,i] = stateControlVector[ju:ju+problem.numControls-1]
    end
    return stateVector, controlVector
end

function getXUFromStateControl(problem, stateControlVector::Array) # get separate state and control vector matricies from input tuple
    stateVector = Array{Union{Missing, typeof(stateControlVector[1]) }}(missing, problem.numStates, problem.numCollocationPoints)
    controlVector = Array{Union{Missing,typeof(stateControlVector[1]) }}(missing, problem.numControls, problem.numCollocationPoints)
    stateVector[:,:] = stateControlVector[1:problem.numStates, :]
    controlVector[:,:] = stateControlVector[problem.numStates + 1, :]
    return stateVector, controlVector
end

function getStateControlFromXU(problem, stateControlVector::Array) # get tuple of statecontrols from a 1D array with states and controls in a column row pattern
    x = ones(1,problem.numCollocationPoints*problem.numStates)
    u = ones(1,problem.numCollocationPoints*problem.numControls)
    for i = 1:problem.numCollocationPoints
        j = (i-1)*(problem.numStates + problem.numControls) + 1
        x[(i-1)*problem.numStates+1:i*problem.numStates] = stateControlVector[j:j+problem.numStates-1]
        u[(i-1)*problem.numControls+1:i*problem.numControls] = stateControlVector[j+problem.numStates:j+problem.numStates+problem.numControls-1]
    end
    return [x u]
end

####################################################################################################################################
### problem specific functions #####################################################################################################
####################################################################################################################################

function dynamicsFunc(stateVector, controlVector) 
    h = stateVector[1,:]
    v = stateVector[2,:]
    m = stateVector[3,:]
    T = controlVector[1,:]
    D(h0,h,v) = 320 * v^2 * exp(-500*((h-h0)/h0))
    g(h0,h) = (h0/h)^2

    x1dot = h
    x2dot = map( (T,hi,v,m) -> (T-D(h[1],hi,v)) / m - g(h[1],hi), T,h,v,m)
    x3dot = -2T
    return [x1dot';x2dot';x3dot']
end
objectiveFunc(stateVector, controlVector) = -stateVector[1,end]

function pathConstraintFunc(stateVector, controlVector) # less than or equal to output value. Output must by 1 D vector of constraints
    T_max = 3.5
    T_min = 0
    h_0 = 1
    m_0 = 1
    m_f = 0.6

    bounds = zeros(typeof(stateVector[1,1]),6,size(stateVector,2))
    bounds[1,:] = -m_f .- stateVector[3,:] # weight constraints
    bounds[2,:] = -m_0 .+ stateVector[3,:]
    bounds[3,:] = -T_min .- controlVector[1,:] # thurst constraints
    bounds[4,:] = -T_max .+ controlVector[1,:]
    bounds[5,:] = -stateVector[2,:] # velocity constraint
    bounds[6,:] = h_0 .- stateVector[1,:] # velocity constraint
    return [bounds...]
end

####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################

numCP = 30
T = 0.2
controlVectorGuess = zeros(1,numCP)
timeStep = ones(1,numCP-1) * T/(numCP-1)
time = pushfirst!(cumsum(timeStep';dims=1)[:,1],0.0)
stateVectorGuess = ones((3,numCP))
stateVectorGuess[2,:] = time .* (1 .- time)
stateVectorGuess[3,:] = -0.4.*time.+1

boundaryConstraints = BoundaryConstraint([0;1;1],[missing;missing;0.6])
problem = TrajProblem(objectiveFunc, dynamicsFunc, pathConstraintFunc, controlVectorGuess, stateVectorGuess, timeStep, boundaryConstraints, 3, 1, numCP)
x,u = solve(problem)


using Plots
plotly()
plot(time, x[1,:])
plot(time, x[2,:])
plot(time, u')

