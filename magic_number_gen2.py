import sympy as s
x,y,z=s.symbols("x,y,z")
n=5
ar=[s.log(x+y)]+[1/(x+y)**i for i in range(1,n)]
ar2=[-s.log(1/x)]+[1/(x)**i for i in range(1,n)]
for i in range(2,n):
    print("Complex y%s{y%s*y1};"%(i,i-1))
dic={i:[] for i in range(n)}
for li,l in enumerate(ar):
    se=s.series(l,x,s.oo,6)
    for mi,m in enumerate(ar2):
        co=se.coeff(m).as_coeff_exponent(y)
        #print(l,m,co)
        a,b=co
        a=float(a)
        if a!=0:
            stb="" if b==0 else "*y%s"%b
            #sta="" if a==1 else "%s"%a
            dic[mi].append("%s%s"%(a,stb)+"*coef[%s]"%li)
for i in range(n):

    print("newcoef[%s]="%i,"+".join(dic[i]).replace("+-","-").replace("1.0*","")+";")

