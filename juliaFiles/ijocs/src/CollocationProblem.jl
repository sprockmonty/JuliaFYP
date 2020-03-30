abstract type AbstractDiscretization end # store the type of collocations

# collocation types
struct TrapezoidalDiscretization <: AbstractDiscretization end
struct EulerDiscretization <: AbstractDiscretization end
struct RungeKuttaDiscretization <: AbstractDiscretization end
struct HermiteDiscretization <: AbstractDiscretization end


mutable struct CollocationProblem
    states::Array{AbstractState}
    controls::Array{AbstractState}
    time::AbstractState # doesn't have to be same lenght as states and controls, but should handle this if that's not the case (i.e. duplicate the time state if only one time state is supplied)
    collocationType::AbstractDiscretization
    CollocationProblem(collocationType::AbstractDiscretization) = new(Array{AbstractState}[], Array{AbstractState}[], Array{AbstractState},collocationType) # check if states have consistent dimensions, if they do, replace array of states with single state?
end

# if dimentions are consitent with previous, convert all states to consistent states, else, leave as abstract

# CollocationProblem methods
add_state!(probem::CollocationProblem, state::AbstractState) =  append!(problem.states)
add_control!(probem::CollocationProblem, state::AbstractState) =  append!(problem.controls)
change_collocation!(problem::CollocationProblem, collocationType::AbstractDiscretization) =  problem.collocationType = collocationType
get_collocation_method(problem::CollocationProblem) = problem.collocationType
get_state_with_name(name::String) = # needs defining


function clean_problem(problem::CollocationProblem) # needs finishing
    problem.stateGuess[time] 
end
