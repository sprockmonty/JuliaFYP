abstract type AbstractState end # rename this to something else that's not state?

abstract type AbstractStateMethod end

struct MultiPolyMethod <: AbstractStateMethod end
struct RadauMethod <: AbstractStateMethod end


mutable struct PolyState{T<:AbstractPoly} <: AbstractState
    poly::T
    name::String
end


# State methods
create_state(p::PolyState, name::String = "") = PolyState(p, name)
create_state(time::Vector, ::RadauMethod, name::String ="" ) # globally using radau
create_state(::MultiPolyMethod, name::String ="" ) = PolyState(GlobalPoly, name)

add_interval(state::PolyState{GlobalPoly}, interval) = # needs finishing 
add_interval(state::PolyState, interval) = # convert to global poly if possible 

set_name(state::AbstractState, name::String) = state.name = name
#update_state!(state::AbstractState, update) = state.stateDefinition = update # check type of update is same as abstract state
    

