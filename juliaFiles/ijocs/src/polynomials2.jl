







### test code

k = [1,5]
l = [1,6]
m = [1,2,3]
n = [1,2,8]
o = [1,2,6,20,25]
p = [1,2,3,4,5]

#first order poly
@benchmark polyfit(k,l) # 1e-6 s
@benchmark create_lagrange_poly(k,l) # 3e-7 s
@benchmark createpolyOld(k,l) # 1.4 e-7 s
@benchmark createpoly(k,l) # 1.4 e-7 s

#second order poly
@benchmark polyfit(m,n) # 1.6 e-6 s
@benchmark create_lagrange_poly(m,n) 
@benchmark createpolyOld(m,n)
@benchmark createpoly(m,n)


#third order poly
@polyfit(o,p)
@benchmark create_lagrange_poly(o,p)
@benchmark createpolyOld(o,p)
@createpoly(o,p)
