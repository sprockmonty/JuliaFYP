abstract type AbstractState end

mutable struct ArrayState <: AbstractState
    stateDefinition::Array
    stateName::String
end

mutable struct ConsistentArrayState <: AbstractState # array dimensions are consistent
    stateDefinition::Array #add check that dimensions are consistent
    stateName::Union{Array{String}, String}
end

mutable struct FuncState <: AbstractState # need to define indexing for this so that indexing is dependant on time states
    stateDefinition::Function # how do we ensure this function has correct inputs?
    stateName::String
end

#mutable struct Time # do we need this?
#    data
#end

# State methods

create_state(state::Array, name::String = "") = ArrayState(state, name)
create_state(state::Array{Number, 2}, name::Union{Array{String}, String} = "") = ConsistentArrayState(state, name)
create_state(state::Function, name::String = "") = FuncState(state, name)

update_state!(state::AbstractState, update) = state.stateDefinition = update # check type of update is same as abstract state
    

