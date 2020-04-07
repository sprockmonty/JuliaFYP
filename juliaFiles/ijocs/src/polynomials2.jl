using ForceImport
using Reduce

macro genPolyFit(index,leny,x,y)
    leny = :($leny)
    symX = []
    symY = []
    A = Array{Expr}(undef, leny, leny)
    for i = 1:leny
        push!(symX, Meta.parse("x$i"))
        push!(symY, Meta.parse("y$i"))
    end
    for i = 1:leny
        for j = 1:leny
            pwr = leny-j
            var = symX[i]
            A[i,j] = :($var^$pwr)
        end
    end
    coefs = Algebra.:*(Algebra.inv(A),symY)
    for i = 1:leny
        str = repr(coefs[i])
        for j = 1:leny
            str = replace(str,"y"*string(j)=>"y"*"["*string(j)*"]")
            str = replace(str,"x"*string(j)=>"x"*"["*string(j)*"]")
        end
        coefs[i] = Meta.parse(chop(str, head=2, tail=1))
    end

    return coefs[:($index)]
end
### macro test code
f = [1,2]
t = [0,4]
@macroexpand @genPolyFit(1, 2,x,y)


function createpoly2(x,y, b::Bound)
    xLength = length(x)
    if xLength == 2
        c1 = @genPolyFit(1, 2,:x,:y)
        c2 = @genPolyFit(2, 2,:x,:y)
        return create_coef_poly([c2,c1], b)
    end
end
@benchmark createpoly(f,t,Bound(-100,100))
###



function createpoly(x,y, b::Bound)
    xLength = length(x)
    if xLength == 2
        pa = (-y[1]x[2] + y[2]*x[1]) / (x[1] - x[2])
        pb = (y[1] - y[2]) / (x[1] - x[2])
        return create_coef_poly([pa,pb], b)
    elseif xLength == 3
        pa = -y[1]*(x[2]*x[3])/(x[1]*x[2] + x[1]*x[3] - x[2]*x[3] - x[1]^2) -y[2]*(x[1]*x[3])/(x[1]*x[2] - x[1]*x[3] + x[2]*x[3] - x[2]^2) + y[3]*(x[1]*x[2])/(x[1]*x[2] - x[1]*x[3] - x[2]*x[3] + x[3]^2)
        pb = y[1]*(x[2] + x[3])/(x[1]*x[2] + x[1]*x[3] - x[2]*x[3] - x[1]^2) + y[2]*(x[1] + x[3])/(x[1]*x[2] - x[1]*x[3] + x[2]*x[3] - x[2]^2) -y[3]*(x[1] + x[2])/(x[1]*x[2] - x[1]*x[3] - x[2]*x[3] + x[3]^2)
        pc =  -y[1]/(x[1]*x[2] + x[1]*x[3] - x[2]*x[3] - x[1]^2) -y[2]/(x[1]*x[2] - x[1]*x[3] + x[2]*x[3] - x[2]^2) + y[3]/(x[1]*x[2] - x[1]*x[3] - x[2]*x[3] + x[3]^2)
        return create_coef_poly([pa,pb,pc], b)
    elseif xLength == 4
        pa = y[1]*(x[2]*x[3]*x[4])/(x[1]^2*x[2] + x[1]^2*x[3] + x[1]^2*x[4] - x[1]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] + x[2]*x[3]*x[4])+ y[2]*(x[1]*x[3]*x[4])/(x[1]*x[2]^2 + x[2]^2*x[3] + x[2]^2*x[4] - x[2]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[3]*x[4] - x[2]*x[3]*x[4])+ y[3]*(x[1]*x[2]*x[4])/(x[1]*x[3]^2 + x[2]*x[3]^2 + x[3]^2*x[4] - x[3]^3 - x[1]*x[2]*x[3] + x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])+ y[4]*(x[1]*x[2]*x[3])/(x[1]*x[4]^2 + x[2]*x[4]^2 + x[3]*x[4]^2 - x[4]^3 + x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])
        pb = -y[1]*(x[2]*x[3] + x[2]*x[4] + x[3]*x[4])/(x[1]^2*x[2] + x[1]^2*x[3] + x[1]^2*x[4] - x[1]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] + x[2]*x[3]*x[4]) -y[2]*(x[1]*x[3] + x[1]*x[4] + x[3]*x[4])/(x[1]*x[2]^2 + x[2]^2*x[3] + x[2]^2*x[4] - x[2]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[3]*x[4] - x[2]*x[3]*x[4]) -y[3]*(x[1]*x[2] + x[1]*x[4] + x[2]*x[4])/(x[1]*x[3]^2 + x[2]*x[3]^2 + x[3]^2*x[4] - x[3]^3 - x[1]*x[2]*x[3] + x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4]) -y[4]*(x[1]*x[2] + x[1]*x[3] + x[2]*x[3])/(x[1]*x[4]^2 + x[2]*x[4]^2 + x[3]*x[4]^2 - x[4]^3 + x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])
        pc = y[1]*(x[2] + x[3] + x[4])/(x[1]^2*x[2] + x[1]^2*x[3] + x[1]^2*x[4] - x[1]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] + x[2]*x[3]*x[4])+ y[2]*(x[1] + x[3] + x[4])/(x[1]*x[2]^2 + x[2]^2*x[3] + x[2]^2*x[4] - x[2]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[3]*x[4] - x[2]*x[3]*x[4])+ y[3]*(x[1] + x[2] + x[4])/(x[1]*x[3]^2 + x[2]*x[3]^2 + x[3]^2*x[4] - x[3]^3 - x[1]*x[2]*x[3] + x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])+y[4]*(x[1] + x[2] + x[3])/(x[1]*x[4]^2 + x[2]*x[4]^2 + x[3]*x[4]^2 - x[4]^3 + x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])
        pd = -y[1]/(x[1]^2*x[2] + x[1]^2*x[3] + x[1]^2*x[4] - x[1]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] + x[2]*x[3]*x[4]) -y[2]/(x[1]*x[2]^2 + x[2]^2*x[3] + x[2]^2*x[4] - x[2]^3 - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[3]*x[4] - x[2]*x[3]*x[4]) -y[3]/(x[1]*x[3]^2 + x[2]*x[3]^2 + x[3]^2*x[4] - x[3]^3 - x[1]*x[2]*x[3] + x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4]) -y[4]/(x[1]*x[4]^2 + x[2]*x[4]^2 + x[3]*x[4]^2 - x[4]^3 + x[1]*x[2]*x[3] - x[1]*x[2]*x[4] - x[1]*x[3]*x[4] - x[2]*x[3]*x[4])
        return create_coef_poly([pa,pb,pc,pd],b)

    end
end
createpoly(x,y) = createpoly(x,y,Bound(-Inf, Inf))


### test code

@benchmark polyfit([-1,0,1,2],[1,2,4,20])
@benchmark createpoly([-1,0,1,2],[1,2,4,20])
@benchmark create_lagrange_poly([-1,0,1,2],[1,2,4,20])

@benchmark polyfit(   [-1,0,1],[1,2,4])
@benchmark createpoly([-1,0,1],[1,2,4])
