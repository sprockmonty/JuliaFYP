abstract type AbstractCollocation end # store the type of collocations

# collocation types
struct TrapizoidalCollocation <: AbstractCollocation end
struct EulerCollocation <: AbstractCollocation end
struct RungeKuttaCollocation <: AbstractCollocation end
struct HermiteCollocation <: AbstractCollocation end


mutable struct CollocationProblem
    stateGuess::Array{AbstractState}
    controlGuess::Array{AbstractState}
    time::AbstractState # doesn't have to be same lenght as states and controls, but should handle this if that's not the case (i.e. duplicate the time state if only one time state is supplied)
    collocationType::AbstractCollocation
    CollocationProblem(collocationType::AbstractCollocation) = new(Array{AbstractState}[], Array{AbstractState}[], Array{AbstractState},collocationType) # check if states have consistent dimensions, if they do, replace array of states with single state?
end

# if dimentions are consitent with previous, convert all states to consistent states, else, leave as abstract

# CollocationProblem methods
add_state!(probem::CollocationProblem, state::AbstractState) =  append!(problem.states)
add_control!(probem::CollocationProblem, state::AbstractState) =  append!(problem.controls)
change_collocation!(problem::CollocationProblem, collocationType::AbstractCollocation) =  problem.collocationType = collocationType
get_collocation_method(problem::CollocationProblem) = problem.collocationType
function clean_problem(problem::CollocationProblem)
    problem.stateGuess[time] 
end
