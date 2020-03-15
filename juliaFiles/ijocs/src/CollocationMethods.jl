function collocation(::TrapizoidalCollocation, stateArray::Array{AbstractState}, controlArray::Array{AbstractState}, timeArray::Array{AbstractState}, Δfunc) # trapizoidal collocation
    Δstate = Δfunc(stateArray, controlArray, timeArray) # return tuple of changing states for time
    map((state, array, time, func) -> , stateArray, timeArray, Δfunc) 
end

function collocation(::TrapizoidalCollocation, stateArray::Array{ConsistentArrayState}, controlArray::Array{ConsistentArrayState}, timeArray::Array{ConsistentArrayState}, Δfunc) # trapizoidal collocation
    
end
