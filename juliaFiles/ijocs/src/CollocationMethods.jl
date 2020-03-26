function collocation(::TrapezoidalDiscretization, stateArray::Array{AbstractState}, controlArray::Array{AbstractState}, timeArray::AbstractState, Δfunc) # trapezoidal collocation
    Δstate = Δfunc(stateArray, controlArray, timeArray) # return tuple of changing states for time
    map((state, array, time, func) -> , stateArray, timeArray, Δfunc) 
end

