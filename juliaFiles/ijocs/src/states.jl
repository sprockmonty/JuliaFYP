abstract type AbstractState end # rename this to something else that's not state?


#mutable struct ArrayState <: AbstractState
#    stateDefinition::Array
#    time::Union{Array, Missing} # need to specify the time at each array point, make it so that this is automatically set as the problem time 
#    stateName::String
#end
#
#mutable struct ConsistentArrayState <: AbstractState # array dimensions are consistent
#    stateDefinition::Array #add check that dimensions are consistent
#    stateName::Union{Array{String}, String}
#end
#
#mutable struct FuncState <: AbstractState # need to define indexing for this so that indexing is dependant on time states
#    stateDefinition::Function # how do we ensure this function has correct inputs?
#    stateName::String
#end

mutable struct TimeState <: AbstractState
    stateDefinition::Array
    stateName::String
end



# State methods

create_state(state::Array, time=missing, name::String = "") = ArrayState(state, time, name)
create_state(state::Function, name::String = "") = FuncState(state, name)

create_time(time::Array, name::String = "") = TimeState(time, name)
update_state!(state::AbstractState, update) = state.stateDefinition = update # check type of update is same as abstract state
    

