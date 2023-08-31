from math import sqrt
import random
from uplifting_beta_recursive import Info,calculate_smallest_root,calculate_alpha
from sympy import Poly, Rational, random_poly,symbols
from itertools import *

x = symbols("x")
y = symbols("y")


d = 10
x_lift = 0
root_dict_list = []
info = Info(d,root_dict_list,x_lift)
d_x = 3
d_y = 3
a = 0
b = 5
frac_exp = 20
p_p_N = 4


#Types of singularites
crunode = Poly(x**2 - y**2,y) 
cusp = Poly(x**3 - y**2,y)
tacnode = Poly(x**4 - y**2,y)
rhampoid_cusp = Poly(x**5 - y**2,y)
#domain = K[x]

def random_poly_puiseux(d_x,d_y,a,b,frac_exp):
   p = Poly(1,y,domain='ZZ[x**(1/10)]')
   lin_facs = []
   for _ in range(d_y):
      root_p_c = Poly(random_poly(x,d_x,a,b),x).all_coeffs()
      root_r_c = [root_p_c[i]*x**(Rational(i,random.randint(1,frac_exp))) for i in range(len(root_p_c))]

      lin_fac = (y - sum(root_r_c))

      lin_facs.append(lin_fac)
      p *= lin_fac
   #p = Poly(lin_facs,y)
   return p,lin_facs

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


p_p, lin_facs = pouteaux_poly(p_p_N)
#p_p = pouteaux_poly_2(p_p_N)
#p_s,lin_facs = swinnerton_dyer(3)
#p_t_ = p_t*p_t*p_t
#p_r = Poly(y**3 + (-5*x**(1/14) - 2*x**(1/10) - 4*x**(2/9) - 2*x**(1/7) - x**(1/6) - 4*x**(1/4) - 10)*y**2 + (8*x**(29/90) + 8*x**(23/63) + 2*x**(13/42) + 16*x**(17/36) + 10*x**(6/35) + 8*x**(11/28) + 20*x**(9/28) + 5*x**(5/21) + 8*x**(7/20) + 4*x**(7/18) + 2*x**(4/15) + 10*x**(3/14) + 30*x**(1/14) + 18*x**(1/10) + 24*x**(2/9) + 18*x**(1/7) + 5*x**(1/6) + 20*x**(1/4) + 29)*y - 32*x**(155/252) - 10*x**(71/210) - 32*x**(103/180) - 40*x**(59/140) - 8*x**(67/126) - 40*x**(29/90) - 40*x**(23/63) - 8*x**(22/45) - 8*x**(13/42) - 16*x**(17/36) - 50*x**(6/35) - 40*x**(13/28) - 32*x**(11/28) - 20*x**(9/28) - 10*x**(8/21) - 5*x**(5/21) - 32*x**(7/20) - 4*x**(7/18) - 8*x**(4/15) - 50*x**(3/14) - 25*x**(1/14) - 40*x**(1/10) - 20*x**(2/9) - 40*x**(1/7) - 4*x**(1/6) - 16*x**(1/4) - 20, y, domain='ZZ[x**(1/10)]')

#frac_exp = 1
p_r, lin_facs = random_poly_puiseux(d_x,d_y,a,b,frac_exp)
#p_r = Poly((y-(1+3*x**(Rational(1,2))+5*x**2))*(y-(1+3*x**(Rational(4,10))+3*x**2))*(y-(4 + 2*x)),y)
#*(y-(2+4*x)*
#p_r(100)




p_r = Poly((y-(3+3*x))*(y-(2+3*x))*(y-(4+3*x)),y)
calculate_smallest_root(p_r,info)
alpha = calculate_alpha(info)
print("The first root is: " + str(alpha))
print("")
