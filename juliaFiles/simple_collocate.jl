using JuMP, Ipopt, Plots

function blockMoving(numCollocationPoints)
    model = Model(Ipopt.Optimizer)

    # Timestep
    hk = 29/(numCollocationPoints-1) # for now set as constant, later add variable timestep
    time = 0:hk:29

    # Variables and start values
    @variable(model, x[i=1:numCollocationPoints])
    set_start_value.(x,time)
    @variable(model, v[1:numCollocationPoints], start = 1)
    @variable(model, u[1:numCollocationPoints], start = 0)

    # Boundary constraints
    fix(x[1], 0)
    fix(x[end], 29)
    fix(v[1], 0)
    fix(v[end], 0)

    # System dynamics
    @NLexpression(model, xDot[i=1:numCollocationPoints], v[i])
    @NLexpression(model, vDot[i=1:numCollocationPoints], u[i])
    
    # Objective function
    @NLexpression(model, objective[i=1:numCollocationPoints], u[i]^2)
    for i in 1:numCollocationPoints - 1
        # Collocation points
        @NLconstraint(model, x[i+1] - x[i] - 0.5 * hk * (xDot[i+1] + xDot[i]) == 0)
        @NLconstraint(model, v[i+1] - v[i] - 0.5 * hk * (vDot[i+1] + vDot[i]) == 0)
    end

    # Objective function interpolation
    @NLexpression(model, objectiveInterp, sum(0.5 * hk * (objective[i] + objective[i+1]) for i in 1:numCollocationPoints - 1))
    @NLobjective(model, Min, objectiveInterp)
    optimize!(model)
    return value.(x), value.(v), value.(u), time
end

x,v,u,times = blockMoving(3000)

# Plotting
plotly()
plot(times, x, xlabel="time (s)", ylabel="x (m)", title="Location plot")
plot(times, v, xlabel="time (s)", ylabel="velocity (m/s)", title="Velocity plot")
plot(times, u, xlabel="time (s)", ylabel="Force (N)", title="Force plot")
