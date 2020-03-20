function collocation(::TrapizoidalCollocation, stateArray::Array{AbstractState}, controlArray::Array{AbstractState}, timeArray::AbstractState, Δfunc) # trapizoidal collocation
    Δstate = Δfunc(stateArray, controlArray, timeArray) # return tuple of changing states for time
    map((state, array, time, func) -> , stateArray, timeArray, Δfunc) 
end

