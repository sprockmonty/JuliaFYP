using JuMP
using Ipopt

const DIFFH = 0.000001 # finite difference difference constant

struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array
end

mutable struct TrajProblem
    objectiveFunc::Function
    dynamicsFunc::Function
    controlVectorGuess::Array
    stateVectorGuess::Array # each row represents a state, maybe use an add state guess function which adds each state guess
    timeStep::Union{Array, Float64} # if only one value input, turn this into an array of that one timestep, must be one less than state length
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
    stateVector::Array # Not touched by user, need to find way of making this private and set it to guess as default
end

# work around to make current problem globally available
mutable struct CurrentProblem
    nullableProblem::Union{TrajProblem,Nothing}
end

const CURRENT_PROBLEM = CurrentProblem(nothing)
setCurrentProblem(problem::TrajProblem) = (CURRENT_PROBLEM.nullableProblem = problem)
getCurrentProblem() = CURRENT_PROBLEM.nullableProblem

function solve(problem::TrajProblem) # multiple dispatch on this function
    problem.stateVector = problem.stateVectorGuess # to ensure if we rerun the problem the state vector will start at guess
    setCurrentProblem(problem)
    model = Model(with_optimizer(Ipopt.Optimizer))
    uIndex = 1:size(problem.controlVectorGuess,2)
    numFinalBounds = length(problem.boundaryConstraints.finalState)
    # assign control varibles
    @variable(model, u[uIndex]) # make this work when u has more than one dimention
    map(i->set_start_value(u[i], problem.controlVectorGuess[i]), uIndex)
    register(model, :boundConstrainFunc!, length(uIndex) + 1, boundConstrainFunc!, ΔboundConstrainFunc)
    register(model, :objectiveFuncInterp, length(uIndex), objectiveFuncInterp, ΔobjectiveFuncInterp)
    @NLconstraint(model,[i= 1:numFinalBounds], boundConstrainFunc!(i,u...)==0.0) # maybe add error tolerance to this?
    @NLobjective(model, Min, objectiveFuncInterp(u...)) 
    optimize!(model)
end

function objectiveFuncInterp(controlVector...)  # this will need to be rewritten when controlVector has more than one dimension
    controlVector = collect(controlVector)'
    problem = getCurrentProblem()
    return sum(0.5 .* problem.timeStep .* (problem.objectiveFunc(controlVector[2:end]) .+ problem.objectiveFunc(controlVector[1:end-1])) ) # should we use a map function here? test with a time
end

function ΔobjectiveFuncInterp(points...)
    # add some stuff here
end

function boundConstrainFunc!(selectedBound,controlVector...) # updates state vector
    controlVector = collect(controlVector)'
    selectedBound = convert(Integer, selectedBound) # for some reason jump turns this into a float, must convert back to int
    problem = getCurrentProblem()
    if problem.stateVectorGuess == problem.stateVector
        problem.stateVector = collocate(problem, problem.stateVectorGuess, problem.controlVectorGuess)
    elseif selectedBound == 1
        problem.stateVector = collocate(problem, problem.stateVector, controlVector)
    end
    return (problem.stateVector[:,end] .- problem.boundaryConstraints.finalState)[selectedBound]
end

function boundConstrainFunc(selectedBound,controlVector...) # does not update state vector
    controlVector = collect(controlVector)'
    selectedBound = convert(Integer, selectedBound) # for some reason jump turns this into a float, must convert back to int
    problem = getCurrentProblem()
    if problem.stateVectorGuess == problem.stateVector
        stateVector = collocate(problem, problem.stateVectorGuess, problem.controlVectorGuess)
    elseif selectedBound == 1
        stateVector = collocate(problem, problem.stateVector, controlVector)
    end
    return (stateVector[:,end] .- problem.boundaryConstraints.finalState)[selectedBound]
end

function ΔboundConstrainFunc(points...)
    # f(boundConstrainFunc)
end

function collocate(problem::TrajProblem, stateVector, controlVector) # make sure timestep vector matches length of state vector
    # stateVector[:,1] = initialState this line might not be needed due to later line, but idk what difference this makes, best to ask
    ΔstateVector = 0.5 * problem.timeStep .* (problem.dynamicsFunc(stateVector[:,2:end], controlVector[:,2:end]) + problem.dynamicsFunc(stateVector[:,1:end-1], controlVector[:,1:end-1])) # how do we ensure dynamicsFunc has the right inputs and outputs if it is user defined?
    return stateVectorOut = cumsum([problem.boundaryConstraints.initialState ΔstateVector], dims=2)
end


####################################################################################################################################
### problem specific functions #####################################################################################################
####################################################################################################################################

dynamicsFunc(stateVector, controlVector) = [stateVector[2,:]' ; controlVector]# think of a better way of doing this to test that this function has the correct inputs and outputs, maybe do a check in the type definition too
objectiveFunc(controlVector) = controlVector.^2


####################################################################################################################################
### Forward difference code - move this to separate file / module later ############################################################
####################################################################################################################################
function forwardDiff(point::Int, func::Function, h::Float64)
    return ( func(point + h) - func(point) ) / h   
end

function forwardDiff(points::Tuple, func::Function, h::Float64)
    println(points) ###
    return ( func(points .+ h) - func(points) ) / h   
end
####################################################################################################################################
### example code  ##################################################################################################################
####################################################################################################################################

controlVectorGuess = zeros(1,30)
stateVectorGuess = [transpose(0:29) ; ones(1,30)]
timeStep = ones(1,29)
boundaryConstraints = BoundaryConstraint([0;0],[0;30])
problem = TrajProblem(objectiveFunc, dynamicsFunc, controlVectorGuess, stateVectorGuess, timeStep, boundaryConstraints, stateVectorGuess)
solve(problem)
