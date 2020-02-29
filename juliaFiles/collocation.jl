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
    model = Model(Ipopt.Optimizer)

    # assign state and control varibles
    @variable(model, x[1:problem.numStates, 1:problem.numCollocationPoints])
    @variable(model, u[1:problem.numControls, 1:problem.numCollocationPoints])
    # set state and control variables to problem guess
    [set_start_value(x[i,j], problem.stateVectorGuess[i,j]) for i = 1:problem.numStates, j = 1:problem.numCollocationPoints] 
    [set_start_value(u[i,j], problem.controlVectorGuess[i,j]) for i = 1:problem.numControls, j = 1:problem.numCollocationPoints] 
    
    #register custom defined functions
    register(model, :collocateConstraint, (problem.numStates + problem.numControls) * problem.numCollocationPoints + 2, collocateConstraint, ΔcollocateConstraint)
    register(model, :objectiveFuncInterp, (problem.numStates + problem.numControls) * problem.numCollocationPoints, objectiveFuncInterp, ΔobjectiveFuncInterp)

    # state boundary constraints
    # initial state
    [fix(x[i,1], problem.boundaryConstraints.initialState[i]) for i = 1:length(problem.boundaryConstraints.initialState)]   # adjust this for lb and ub
    # final state
    [fix(x[i,end], problem.boundaryConstraints.finalState[i]) for i = 1:length(problem.boundaryConstraints.finalState)]   # adjust this for lb and ub

    for i = 1:problem.numStates
        for j = 1:problem.numCollocationPoints - 1
            @NLconstraint(model, collocateConstraint(i,j,x...,u...) == 0)
        end
    end
    
    @NLobjective(model, Min, objectiveFuncInterp(x...,u...)) 
    optimize!(model)
end

function objectiveFuncInterp(stateControlVector...)  
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)

    # return interpolated integral
    return sum(0.5 .* problem.timeStep[1] .* (problem.objectiveFunc(controlVector[2:end]) .+ problem.objectiveFunc(controlVector[1:end-1])) ) # should we use a map function here? test with a time, also needs to be rewritten for control vector with multiple dimensions
end

function objectiveFuncInterp(stateControlVector::Array)  # can we type union this? so we don't have to define two functions, only one
    problem = getCurrentProblem()
    # assemble state and control vector from tuple inputs
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    # return interpolated integral
    return [sum(0.5 .* problem.timeStep[1] .* (problem.objectiveFunc(controlVector[2:end]) .+ problem.objectiveFunc(controlVector[1:end-1])))]  # should we use a map function here? test with a time, also needs to be rewritten for control vector with multiple dimensions
end

function ΔobjectiveFuncInterp(g, stateControlVector...)
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = ForwardDiff.jacobian(objectiveFuncInterp, [stateVector; controlVector])
    g[:] = jacobian[:]
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
    g[1:2] .= 0 # gradient of xr and xc is zero as they're just indexes and will not be changed
    ζr = convert(Int, ζr)
    ζc = convert(Int, ζc) # convert from jumps conversion to float64
    stateVector, controlVector = getXUFromStateControl(problem, stateControlVector)
    jacobian = ForwardDiff.jacobian(collocateConstraint, [stateVectorGuess; controlVectorGuess])
    jacobianRow = (ζc - 1) * ζr + ζr # get the row of the jacobian that matches our output function
    g[3:end] = jacobian[jacobianRow,:]
    return g
end

function getXUFromStateControl(problem, stateControlVector::Tuple) # get separate state and control vector matricies from input tuple
    stateControlVector = collect(stateControlVector)
    stateVector = zeros(problem.numStates, problem.numCollocationPoints)
    controlVector = zeros(problem.numControls, problem.numCollocationPoints)
    for i = 1:length(stateControlVector) ÷ (problem.numStates + problem.numControls)
        j = i * (problem.numStates + problem.numControls) - (problem.numStates + problem.numControls) + 1
        stateVector[:,i] = stateControlVector[j:j+problem.numStates-1]
        controlVector[:,i] = stateControlVector[j+problem.numStates:j+problem.numStates+problem.numControls-1]
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

####################################################################################################################################
### problem specific functions #####################################################################################################
####################################################################################################################################

dynamicsFunc(stateVector, controlVector) = [stateVector[2,:]' ; controlVector]# think of a better way of doing this to test that this function has the correct inputs and outputs, maybe do a check in the type definition too
objectiveFunc(controlVector) = controlVector.^2


####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################

controlVectorGuess = zeros(1,30)
stateVectorGuess = [transpose(0:29) ; ones(1,30)]
timeStep = ones(1,29)
boundaryConstraints = BoundaryConstraint([0;0],[0;30])
problem = TrajProblem(objectiveFunc, dynamicsFunc, controlVectorGuess, stateVectorGuess, timeStep, boundaryConstraints, 2, 1, 30)
solve(problem)

using Plots
plotly()
