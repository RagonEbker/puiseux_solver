from math import sqrt
import random
import time
from uplifting_beta_recursive import Info,calculate_smallest_root,calculate_alpha
from sympy import Poly, Rational, random_poly,symbols,nsimplify,roots
from itertools import *

x = symbols("x")
y = symbols("y")





#Types of singularites
crunode = Poly(x**2 - y**2,y) 
cusp = Poly(x**3 - y**2,y)
tacnode = Poly(x**4 - y**2,y)
rhampoid_cusp = Poly(x**5 - y**2,y)
wiki = Poly(y**2-(x**3+x**2),y)
#domain = K[x]
##n_s number of singularities needs to be smaller than d_x
#a and b are ranges for random numbers

def random_poly_puiseux(d_x,d_y,a,b,frac_exp,n_s):
    p = Poly(1,y,domain='ZZ[x**(1/10)]')
    lin_facs = []
    roots = []
    for _ in range(d_y):
        root_p_c = Poly(random_poly(x,d_x,a,b),x).all_coeffs()
        root_r_c = [root_p_c[i]*x**(Rational(i,random.randint(1,frac_exp))) for i in range(len(root_p_c))]
        roots.append(root_r_c)
       
    for i in range(n_s):
        sing_val = random.randint(a,b)
        for j in range(random.randint(1,d_y)):
            roots[j][i] = sing_val*x**Rational(i,frac_exp)
    for root in roots:
        lin_fac = (y - sum(root))

        lin_facs.append(lin_fac)
        p *= lin_fac
    print("d_x: " + str(d_x))
    print("d_y: " + str(d_y))
    print("Roots: ")
    print(repr(roots))
    print("Poly: ")
    print(repr(p))
    return p,lin_facs


#this is example 3 from improving complexity bounds
def pouteaux_poly(N,domain = False):
   if (domain):
      p = Poly(1,y,domain=domain)
   else:
      p = Poly(1,y)
   lin_facs = []

   for i in range(N):
      sum_k = 0
      for j in range(i):
         sum_k += x**j
      s_k = (2*x**i + sum_k)
      lin_fac = (y - s_k)

      lin_facs.append(lin_fac)
      p *= lin_fac
   return p,lin_facs

#this is example 1 from improving complexity bounds
def pouteaux_poly_2(N,domain = False):
   if domain:
      p = Poly((y - x**N )*(y**2 - x**3 ),y,domain)
   else:
      p = Poly((y - x**N )*(y**2 - x**3 ),y)
   return p

def pouteaux_poly_3(N,domain = False):
   print("")

def swinnerton_dyer(N):
   combos = list(product([0,1],repeat=N))
   lin_facs = []
   p = Poly(1,y)
   for combo in combos:
      lin_fac = y
      for i in range(len(combo)):
         lin_fac += (-1)**combo[i] * sqrt(x+i)
      lin_facs.append(lin_fac)
      p*=lin_fac
   return p,lin_facs

#precision
d = 20
x_lift = 0
root_dict_list = []
info = Info(d,root_dict_list,x_lift)
d_x = 2
d_y = 3
n_s = 2
a = 0
b = 5
frac_exp = 5
p_p_N = 5

#p_r, lin_facs = pouteaux_poly_2(p_p_N)

#p_r = pouteaux_poly_2(p_p_N)
#p_s,lin_facs = swinnerton_dyer(3)
#p_t_ = p_t*p_t*p_t
#p_r = Poly(y**3 + (-5*x**(1/14) - 2*x**(1/10) - 4*x**(2/9) - 2*x**(1/7) - x**(1/6) - 4*x**(1/4) - 10)*y**2 + (8*x**(29/90) + 8*x**(23/63) + 2*x**(13/42) + 16*x**(17/36) + 10*x**(6/35) + 8*x**(11/28) + 20*x**(9/28) + 5*x**(5/21) + 8*x**(7/20) + 4*x**(7/18) + 2*x**(4/15) + 10*x**(3/14) + 30*x**(1/14) + 18*x**(1/10) + 24*x**(2/9) + 18*x**(1/7) + 5*x**(1/6) + 20*x**(1/4) + 29)*y - 32*x**(155/252) - 10*x**(71/210) - 32*x**(103/180) - 40*x**(59/140) - 8*x**(67/126) - 40*x**(29/90) - 40*x**(23/63) - 8*x**(22/45) - 8*x**(13/42) - 16*x**(17/36) - 50*x**(6/35) - 40*x**(13/28) - 32*x**(11/28) - 20*x**(9/28) - 10*x**(8/21) - 5*x**(5/21) - 32*x**(7/20) - 4*x**(7/18) - 8*x**(4/15) - 50*x**(3/14) - 25*x**(1/14) - 40*x**(1/10) - 20*x**(2/9) - 40*x**(1/7) - 4*x**(1/6) - 16*x**(1/4) - 20, y, domain='ZZ[x**(1/10)]')

#frac_exp = 1
p_r, lin_facs = random_poly_puiseux(d_x,d_y,a,b,frac_exp,n_s)
#p_r = Poly((y-(5+ 2*x**(1/5)+x))*(y-(5+ 2*x**(1/5)+x))*(y-(5+ 4*x**(1/4)+5*x**(2/3)+ 3*x**(3/2))),y)
p_r = Poly((y-(5+3*x+ x**(1/2)))*(y-(5+3*x)),y)
p_r = Poly(nsimplify(p_r),y)


#p_r = Poly((y-(1+3*x**(Rational(1,2))+5*x**2))*(y-(1+3*x**(Rational(4,10))+3*x**2))*(y-(4 + 2*x)),y)
#*(y-(2+4*x)*
#p_r(100)


#p_r = Poly((y-(3+3*x))*(y-(2+3*x))*(y-(4+3*x)),y)
start_time = time.time()
#p_r = Poly(y**3 + (-4*x**Rational(3,5) - 10*x**Rational(1,5) - 5*x**Rational(2,3) - 10*x - 7)*y**2 + (12*x**Rational(19,15) + 35*x**Rational(13,15) + 36*x**Rational(8,5) + 55*x**Rational(6,5) + 40*x**Rational(4,5) + 16*x**(3/5) + 25*x**(2/5) + 50*x**(1/5) + 48*x**(5/3) + 6*x**(4/3) + 23*x**(2/3) + 9*x**2 + 49*x + 16)*y - 108*x**(34/15) - 105*x**(28/15) - 30*x**(23/15) - 60*x**(22/15) - 24*x**(19/15) - 50*x**(16/15) - 85*x**(13/15) - 45*x**(11/5) - 180*x**(9/5) - 72*x**(8/5) - 25*x**(7/5) - 155*x**(6/5) - 80*x**(4/5) - 16*x**(3/5) - 75*x**(2/5) - 60*x**(1/5) - 27*x**(8/3) - 54*x**(7/3) - 123*x**(5/3) - 12*x**(4/3) - 26*x**(2/3) - 18*x**2 - 158*x - 12, y, domain='EX')
#p_r = Poly(nsimplify(p_r),y)
#p_r = Poly((y - (2+5*x**(1/5)+9*x))*(y - (2+5*x**(1/5) + 3*x**(2/3)))*(y - (3+x+2*x**(2/3) + 4*x**0.6)),y)
p_r = Poly((y-(1+x+x**2+x**3))*(y-(0.5+2*x)),y)
f = Poly(2*x**4 + x**2*y+4*x*y**2 + 4*y**3,y)
g = Poly((y-(1+x+x**2+x**3+x**4+x**5+x**6))*(y-(1+x+x**2+2*x**3)),y)
g = Poly((y-(1+x))*(y-(1+x)),y)
calculate_smallest_root(g,info)
print("--- %s seconds ---" % (time.time() - start_time))

alpha = calculate_alpha(info)
print("The first root is: " + str(alpha))
print("")
