from math import floor, log10, hypot
from math import pi as 𝜋, gamma as 𝛤

def sqp_d0(a,b,c,m,n):
    # Diameter of bounding sphere of superquadric
    # See model/des/sq_properties.f
    #
    # test case1, ellipsoidal, sqp_d0=2*a=0.002
    #    a=0.001; b=0.0005; c=0.0002; m=2.0; n=2.0
    # test case2, cylindrical, sqp_d0=0.0020396078 (this is approximately a cylinder)
    #    a=0.001; b=0.0005; c=0.0002; m=200000; n=2.0
    # test case3, cubic, sqp_d0=2*(0.0005^2+0.001^2+0.0002^2)^0.5=0.0022715634 (this is approximately a cube)
    #     a=0.001; b=0.0005; c=0.0002; m=200000; n=200000
    a,b = max(a,b), min(a,b)
    if m == n == 2:
        return 2*max(a,c)
    if n == 2:
        return 2*hypot(a,c)
    if n < 2.1:
        n = 2.1
    if m == 2:
        𝛼 = 0
    else:
        if abs(m-2) < 0.1:
            m = 2.1 if m>2 else 1.9
        𝛼 = (b/a)**(2/(m-2))

    𝛾 = (1+𝛼**m)**(n/m-1)
    𝛽 = (𝛾*c**2 / a**2)**(1/(n-2))
    𝜁 = 1 / ((1+𝛼**m)**(n/m)+𝛽**n)**(1/n)
    x = a*𝜁
    y = b*𝛼*𝜁
    z = c*𝛽*𝜁
    return 2*hypot(x,y,z)

def sqp_volume(a,b,c,m,n):
    e1, e2 = 2/n, 2/m
    def B(a,b): return  𝛤(a)*𝛤(b)/𝛤(a+b) # Beta function
    #Jaklič, A., Leonardis, A., Solina, F. (2000). Superquadrics and Their Geometric Properties.
    # formula 2.60 https://cse.buffalo.edu/~jryde/cse673/files/superquadrics.pdf
    V = 2 * a*b*c * e1*e2 * B(e1/2 + 1, e1) * B(e2/2, e2/2)
    return V
