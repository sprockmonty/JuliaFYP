function createpoly(x,y, b::Bound)
    xLength = length(x)
    if length(x) == 2
        m = (y[1] - y[2]) / (x[1] - x[2])
        c = (-y[1]x[2] + y[2]*x[1]) / (x[1] - x[2])
        return create_coef_poly([c,m], b)
    end
end

