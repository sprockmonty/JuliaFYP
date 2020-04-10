
function createpoly(y)
    create_LGR_poly(y)
end

function createpoly(c, b::Bound)
    create_coef_poly(c)
end

function createpoly(x,y)
    if length(x) < 6
        create_vand_poly(x,y)
    else
        create_lagrange_poly(x,y)
    end
end

function createpoly(x,y,b::Bound)
    if size(x) < 6
        create_vand_poly(x,y,b)
    else
        create_lagrange_poly(x,y,b)
    end
end

### example code

createpoly([5,6,7,8]) # create LGR poly
createpoly([0,1,2], unbounded()) # create unbounded coefficient poly
createpoly([0,1,2],[4,25,36]) # create low order polynomial
createpoly(1:10,[4,5,2,3,4,5,6,7,8,10]) # create high order polynomial
createpoly(LGRSpecifier(), [1,2,3]) # create LGR poly with LGR specifier
