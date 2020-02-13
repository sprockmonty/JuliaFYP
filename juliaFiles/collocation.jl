struct TrajProblem
    dynamicsFunc::Function
    controlVectorGuess::Array
    stateVectorGuess::Array # each row represents a state, maybe use an add state guess function which adds each state guess
    timeStep::Union{Array, Float64} # if only one value input, turn this into an array of that one timestep, must be one less than state length
    boundaryConstraints::BoundaryConstraint # add checks that these are dimensionally correct
end

struct BoundaryConstraint # make a default setting if unbounded
    initialState::Array # make sure these are Nx1 size matricies
    finalState::Array
end

function solve(problem::TrajProblem) # multiple dispatch on this function
    
end


function collocate(dynamicsFunc, stateVector, controlVector, initialState, timeStep) # make sure timestep vector matches length of state vector
    # stateVector[:,1] = initialState this line might not be needed due to later line, but idk what difference this makes, best to ask
    ΔstateVector = 0.5 * timeStep .* (dynamicsFunc(stateVector[:,2:end], controlVector[:,2:end]) + dynamicsFunc(stateVector[:,1:end-1], controlVector[:,1:end-1])) # how do we ensure dynamicsFunc has the right inputs and outputs if it is user defined?
    stateVectorOut = cumsum([initialState ΔstateVector], dims=2)
end



####################################################################################################################################

# example code

controlVectorGuess = zeros(1,30)
stateVectorGuess = [transpose(0:29) ; ones(1,30)]
timeStep = ones(1,29)

boundayConstraints = BoundaryConstraint([0;0],[0;30])

TrajProblem(dynamicsFunc, controlVectorGuess, stateVectorGuess, timeStep, boundaryConstraints)

solve(TrajProblem)


