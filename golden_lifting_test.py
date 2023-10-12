from math import sqrt
import random
import time
from uplifting_beta_recursive import Info,calculate_smallest_root,calculate_alpha
from sympy import Poly, Rational, random_poly,symbols,nsimplify,roots
from itertools import *

x = symbols("x")
y = symbols("y")





#Different types of singularites
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


#this is example 3 from DOI: 10.1145/2755996.2756650
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

#this is example 1 from DOI: 10.1145/2755996.2756650
def pouteaux_poly_2(N,domain = False):
   if domain:
      p = Poly((y - x**N )*(y**2 - x**3 ),y,domain)
   else:
      p = Poly((y - x**N )*(y**2 - x**3 ),y)
   return p

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

#Here we create the initial info object
#precision
d = 20
x_lift = 0
root_dict_list = []
info = Info(d,root_dict_list,x_lift)

#These are the parameters to create a random polynomial
d_x = 2
d_y = 3
n_s = 2
a = 0
b = 5
frac_exp = 5


p_r, lin_facs = random_poly_puiseux(d_x,d_y,a,b,frac_exp,n_s)

#A few more example polynomials
#p_r = Poly((y-(5+ 2*x**(1/5)+x))*(y-(5+ 2*x**(1/5)+x))*(y-(5+ 4*x**(1/4)+5*x**(2/3)+ 3*x**(3/2))),y)
#p_r = Poly((y-(1+3*x**(Rational(1,2))+5*x**2))*(y-(1+3*x**(Rational(4,10))+3*x**2))*(y-(4 + 2*x)),y)
#p_r = Poly((y-(3+3*x))*(y-(2+3*x))*(y-(4+3*x)),y)
# p_r = Poly((y-(1+x+x**2+x**3))*(y-(0.5+2*x)),y)
# f = Poly(2*x**4 + x**2*y+4*x*y**2 + 4*y**3,y)
# g = Poly((y-(1+x+x**2+x**3+x**4+x**5+x**6))*(y-(1+x+x**2+2*x**3)),y)
# g = Poly((y-(1+x))*(y-(1+x)),y)

start_time = time.time()

calculate_smallest_root(g,info)
print("--- %s seconds ---" % (time.time() - start_time))

alpha = calculate_alpha(info)
print("The first root is: " + str(alpha))
