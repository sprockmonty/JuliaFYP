using JuMP, Ipopt


function blockMoving()
    model = Model(Ipopt.Optimizer)

    numCollocationPoints = 30
    hk = 1/29 # for now set as constant, later add variable timestep

    # Variables and start values
    @variable(model, x[i=1:numCollocationPoints])
    set_start_value.(x,0:1/29:1)
    @variable(model, v[1:numCollocationPoints], start = 1)
    @variable(model, u[1:numCollocationPoints], start = 0)

    # Boundary constraints
    fix(x[1], 0)
    fix(x[30], 1)
    fix(v[1], 0)
    fix(v[30], 0)

    # System dynamics
    @NLexpression(model, xDot[i=1:numCollocationPoints], v[i])
    @NLexpression(model, vDot[i=1:numCollocationPoints], u[i])
    
    # objective function
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

    println(value.(u))
    println(value.(x))
end

blockMoving()

