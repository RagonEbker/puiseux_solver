import math
from sympy import Pow, sqrt
import numpy as np
from sympy.polys.ring_series import ring,PolyElement,rs_series_inversion,rs_series,_rs_series,rs_trunc
from sympy import RR, parse_expr,Matrix,simplify,Mod,Add,Mul,Pow,real_roots,symbols,Poly,roots,nan,re,eye,Identity,series,solve,Rational,div,degree,integer_nthroot,binomial,GF
import random
from sympy.simplify.simplify import nsimplify
import collections
from sympy.core.numbers import One,sympify,Float,S
from sympy.core.rules import Transform
import signal
from sympy.polys.specialpolys import random_poly
from itertools import *
import scipy.spatial


class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)

z = symbols('z')
y = symbols('y')
x = symbols('x')
a = symbols('a')
b = symbols('b')
c = symbols("c")
d = symbols("d")
e = symbols("e")
f = symbols("f")
g = symbols("g")
# p = 
K = GF(131)
##TODO PRECISION EINBAUEN
class Info:
  #Precision d
  #x_lifts
  #all the info about the roots in root_dict_list
  #alpha is our first root
  def __init__(self, d,root_dict_list,x_lift):
    self.d = d
    self.root_dict_list = root_dict_list
    self.x_lift = x_lift




#Helper functions
def intify(expr):
   if(is_expr(expr)):

      floats = S(expr).atoms(Float)
      ints = [i for i in floats if int(i) == i]
      return expr.xreplace(dict(zip(ints, [int(i) for i in ints])))
   else:
      if int(expr) == expr:
         return int(expr)
      else:
         return expr

def is_expr(expr):
   expr = sympify(expr)
   symbols = expr.free_symbols
   if not symbols:
      return False
   else:
      return True
   
def calculate_alpha(info):
   return(sum([list(x)[0] for x in info.root_dict_list]))

def calculate_h_shift(info):
   return list(info.root_dict_list[-1])[0]

def shift_vertically(p,mtp):   
   coeffs = p.all_coeffs()
   l_coeffs =len(coeffs)
   new_poly = sum([coeffs[i]*mtp**i*y**(l_coeffs-i) for i in range(l_coeffs)])
   print("We shift vertically by: " + str(mtp))
   return Poly(new_poly,y)

def order_roots(p_r):
   p_r = collections.OrderedDict(sorted(p_r.items(), key=lambda _:_[0].as_coeff_exponent(x)[0]))
   return p_r
#Checks if our Polynomial consist of exactly one root
##TODO For the lift we need the multiplicity
##REWRITE THIS
##DO THIS WITH THE ROOTS FUNCTION SOOO MUCH EASIER


def check_one_root(p):
   r = roots(p)
   # if(len(list(r)) == 1):
   #    return (r,True)
   # else:
      
   return (r, False)

##Returns whole dict of root, smallest root that is not zero, if it exists, otherwise 0
##and mtpcty of that root
def get_sub_x_root(p,x_lift):
   p_x = Poly(p.subs(x,0),y)
   p_r_old = roots(p_x)
   p_r = dict()
   for key, value in p_r_old.items():
    # do something with value
    ##maybe say if value != 0
    if (key): p_r[key*x**x_lift] = value

   if (p_r):
      p_r = collections.OrderedDict(sorted(p_r.items(), key=lambda _:_[0].as_coeff_exponent(x)[0]))
      s_r = list(p_r)[0]
      #s_r = p_r_list[0]
      mtpcty = p_r[s_r]
   #smallest root if not zero
   else:
      s_r = 0
      mtpcty = p_r_old[0]
      
   return p_r,s_r,mtpcty




def shift_horizontally(p,shift):
   print("We shift horizontally by: " + str(shift))
   return Poly((p(z+shift)).subs(z,y),y)


def golden_lifting(p_shift,mtpcty,info):
   shifted_coeffs = list(reversed(p_shift.coeffs()[1:]))
   cutoff_coeffs = []
   #SET MULTIPLICTIY TO 1
   for i in range(min(len(shifted_coeffs),mtpcty+1)):
      coeff = shifted_coeffs[i]
      terms = coeff.as_ordered_terms()
      lowest_term = min(terms, key=lambda term: term.as_coeff_exponent(x)[1])
      cutoff_coeffs.append(lowest_term)
   cutoff_coeffs = list(reversed(cutoff_coeffs))
   new_poly = Poly.monic(Poly(cutoff_coeffs,y))
   
   if (mtpcty > 1):
         #here we test if new_poly consists of just a single root
         check_test_mtpcty = False
         try:
            timeout = 10
            #signal.alarm(timeout)
            try: 
               test_mtpcty = check_one_root(new_poly)
               check_test_mtpcty = test_mtpcty[1]
            except TimeoutException:
               pass
         except:
            pass
         if (check_test_mtpcty):
            r = test_mtpcty[0]
            info.root_dict_list.append(order_roots(r))
         else:
            calculate_smallest_root_q_x(new_poly,info)
   else:
      r = roots(new_poly)
      info.root_dict_list.append(order_roots(r))


def calculate_smallest_root_q_x(p,info):

   #TOODO ADD ROOT DICT TO GLOBAL OBJECT
   root_dict,shift_number,_ = get_sub_x_root(p,info.x_lift)
   if(shift_number):
      info.root_dict_list.append(root_dict)
      return
   if (check_done(p,info)):
      return

   
   #lift the polynomial horizontally when all its roots are zeroes
   else:
      
      slopes = get_newton_slopes(p)
      slopes = [x for x in slopes if x!= 0]
      min_slope = nsimplify(min(slopes))
      p_shift = shift_vertically(p,x**(min_slope)) 
      info.x_lift = -min_slope
      info.d = info.d-1
      calculate_smallest_root_q_x(p_shift,info)

def check_done(p,info):
   return (p(nsimplify(calculate_alpha(info), tolerance=0.001, rational=True)) == 0)

def calculate_initial_shift(p,info):
   root_dict,shift_number, _ = get_sub_x_root(p,info.x_lift)
   if (shift_number):
      info.root_dict_list.append(root_dict)
      return p

   #lift the polynomial horizontally when all its roots are zeroes
   else:
      slopes = get_newton_slopes(p)
      slopes = [x for x in slopes if x != 0]
      min_slope = min(np.abs(slopes))
      p_shift = shift_vertically(p,x**(-min_slope)) 
      info.x_lift = min_slope
      info.d = info.d-1
      info = calculate_smallest_root_q_x(p_shift,info)
      return p_shift


#Calculates the "smallest" root of a polynomial over the Puiseux field
def calculate_smallest_root(p,info):
   
   p_shift = calculate_initial_shift(p,info)
   if (check_done(p,info)):
      return    
   print(calculate_alpha(info))

   for _ in range(1,info.d):

      info.d = info.d-1
      p_shift = shift_horizontally(p_shift,calculate_h_shift(info))
      mtpcty = list(info.root_dict_list[-1].values())[0]
      
      golden_lifting(p_shift,mtpcty,info)
      if (check_done(p,info)):
         return 
      if (p_shift == 0):
         return 
   return 

def get_lower(polygon):
   minx = np.argmin(polygon[:, 0])
   maxx = np.argmax(polygon[:, 0]) + 1
   if minx >= maxx:
      lower_curve = np.concatenate([polygon[minx:], polygon[:maxx]])
   else:
      lower_curve = polygon[minx:maxx]
   return lower_curve

def slope(point1,point2):
   result = (point2[1]-point1[1])/(point2[0]-point1[0])
   return result

def get_newton_slopes(f):

   coeff_dict = f.as_expr().as_coefficients_dict(y)


   newton_points = []
   for k,v in coeff_dict.items():
      if (v.atoms(Pow)):
         if (isinstance(k,One)):
               newton_points.append([0,v.atoms(Pow).pop().as_base_exp()[1]])
         else:
               newton_points.append([k.as_base_exp()[1], v.atoms(Pow).pop().as_base_exp()[1]])
      else:
         if (x in v.free_symbols):
            if (k.as_base_exp()[0] == 1):
               newton_points.append([0, v.atoms().pop().coeff(x)])
            else:
               newton_points.append([k.as_base_exp()[1],1])
         else:
            if (k.as_base_exp()[0] == 1):
                  newton_points.append([0, v.atoms().pop().coeff(x)])
            else:
                  newton_points.append([k.as_base_exp()[1],0])


   points = np.asarray(newton_points)
   terminated = False
   try:
      hull = scipy.spatial.ConvexHull(points)

      lower_curve = get_lower(points[hull.vertices])
      slopes = [slope(lower_curve[i],lower_curve[i+1]) for i in range(len(lower_curve)-1)]
      terminated = True
   except:
      pass
   if (not terminated):
      slopes = [slope(points[0],points[-1])]
   return slopes

def random_poly_puiseux(d_x,d_y,a,b,frac_exp,domain):
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
##TODO EXAMPLES
##Examples that have the properties
##for each combination we want an example
#different exponent, same exponent, different coefficient with same exponent, same coefficient with same exponent
#integer exponents, rational exponents
#Finite Fields, ZZ,QQ,RR,CC
#Different precisions
#textbook examples
#competitve examples (with times)

#def random_poly()

d = 10
x_lift = 0
root_dict_list = []
info = Info(d,root_dict_list,x_lift)
d_x = 5
d_y = 5
a = 0
b = 5
frac_exp = 20
#p_r, lin_facs = random_poly_puiseux(d_x,d_y,a,b,frac_exp)
p_p_N = 10
#domain = K[x]
#p_p, lin_facs = pouteaux_poly(p_p_N)
#p_p = pouteaux_poly_2(p_p_N)
p_s,lin_facs = swinnerton_dyer(2)
#p_t = Poly((y-(1+3*x**(Rational(1,2))+5*x**2))*(y-(1+3*x**(Rational(4,10))+3*x**2))*(y-(4 + 2*x)),y)
#p_t_ = p_t*p_t*p_t


calculate_smallest_root(p_s,info)
alpha = calculate_alpha(info)
print("The first root is: " + str(alpha))
print("")
