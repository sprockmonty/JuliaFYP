using JuMP
using Ipopt

struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array
end

mutable struct TrajProblem
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
    setCurrentProblem(problem)
    model = Model(with_optimizer(Ipopt.Optimizer))
    uIndex = 1:size(problem.controlVectorGuess,2)
    numFinalBounds = length(problem.boundaryConstraints.finalState)
    # assign control varibles
    @variable(model, u[uIndex]) # make this work when u has more than one dimention
    map(i->set_start_value(u[i], problem.controlVectorGuess[i]), uIndex)
    
    register(model, :boundConstrainFunc!, numFinalBounds + 1, boundConstrainFunc!, autodiff=true)
    for i = 1:numFinalBounds 
        @NLconstraint(model,[i], boundConstrainFunc!(i,u...)==0) # maybe add error tolerance to this?
    end
end


function boundConstrainFunc!(selectedBound,controlVector...)
    problem = getCurrentProblem()
    if problem.stateVectorGuess == problem.stateVector
        problem.stateVector = collocate(problem, problem.stateVectorGuess, problem.controlVectorGuess)
    elseif selectedBound == 1
        problem.stateVector = collocate(problem, problem.stateVector, controlVector)
    end
    return [problem.stateVector[:,end] .- problem.boundaryConstraints.finalState][selectedBound]
end

function collocate(problem::TrajProblem, stateVector, controlVector) # make sure timestep vector matches length of state vector
    # stateVector[:,1] = initialState this line might not be needed due to later line, but idk what difference this makes, best to ask
    ΔstateVector = 0.5 * problem.timeStep .* (problem.dynamicsFunc(stateVector[:,2:end], controlVector[:,2:end]) + problem.dynamicsFunc(stateVector[:,1:end-1], controlVector[:,1:end-1])) # how do we ensure dynamicsFunc has the right inputs and outputs if it is user defined?
    return stateVectorOut = cumsum([problem.boundaryConstraints.initialState ΔstateVector], dims=2)
end

function dynamicsFunc(stateVector, controlVector)
    return ∂stateVector = [stateVector[2,:]' ; controlVector]
end
####################################################################################################################################

# example code
controlVectorGuess = zeros(1,30)
stateVectorGuess = [transpose(0:29) ; ones(1,30)]
timeStep = ones(1,29)
boundaryConstraints = BoundaryConstraint([0;0],[0;30])

problem = TrajProblem(dynamicsFunc, controlVectorGuess, stateVectorGuess, timeStep, boundaryConstraints, stateVectorGuess)
