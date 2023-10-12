from sympy import Pow
import numpy as np
from sympy.polys.ring_series import PolyElement,rs_series_inversion,rs_series,_rs_series,rs_trunc
from sympy import RR, parse_expr,Matrix,simplify,Mod,Add,Mul,Pow,real_roots,symbols,Poly,roots,nan,re,eye,Identity,series,solve,Rational,div,degree,integer_nthroot,binomial,GF,parse_expr
from sympy.simplify.simplify import nsimplify
import collections
from sympy.core.numbers import One,sympify,Float,S
import signal
from sympy.polys.specialpolys import random_poly
import scipy.spatial
from sympy.polys.polyfuncs import horner


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

def shift_horizontally(p,shift):
   #print("We shift horizontally by: " + str(shift))
   q = horner(p,wrt=y)
   return_poly_1 = Poly(q.subs(y,y+shift),y)
   return return_poly_1

def shift_horizontally_2(p,shift,mtpcty):
   #print("We shift horizontally by: " + str(shift))
   coeffs = p.all_coeffs()
   return_poly = Poly(0,y)
   for i in range(mtpcty):
      a = coeffs[-(i+1)]
      new_term = a*(y+shift)**i
      return_poly += new_term
   return return_poly

def shift_vertically(p,mtp):   
   coeffs = p.all_coeffs()
   l_coeffs =len(coeffs)
   new_poly = sum([coeffs[i]*mtp**i*y**(l_coeffs-i) for i in range(l_coeffs)])
   #print("We shift vertically by: " + str(mtp))
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
   #smallest root if not zero
   else:
      s_r = 0
      
   return p_r,s_r

def golden_lifting(p_shift,mtpcty,info):
   shifted_coeffs = list(reversed(p_shift.coeffs()))
   cutoff_coeffs = []
   #SET MULTIPLICTIY TO 1
   for i in range(mtpcty+1):
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
   root_dict,shift_number = get_sub_x_root(p,info.x_lift)
   if(shift_number):
      info.root_dict_list.append(root_dict)
      return
   if (check_done(p,info)):
      return

   
   #lift the polynomial horizontally when all its roots are zeroes
   else:
      
      slopes = get_newton_slopes(p)
      slopes = [x for x in slopes if x!= 0]
      min_slope = nsimplify(min(slopes,key=abs))
      p_shift = shift_vertically(p,x**(min_slope)) 
      info.x_lift = -min_slope
      info.d = info.d-1
      calculate_smallest_root_q_x(p_shift,info)

def check_done(p,info):
   return (p(nsimplify(calculate_alpha(info), tolerance=0.001, rational=True)) == 0)

def calculate_initial_shift(p,info):
   root_dict,shift_number = get_sub_x_root(p,info.x_lift)
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

   for i in range(1,info.d):

      info.d = info.d-1
      p_shift = shift_horizontally(p_shift,calculate_h_shift(info))
      mtpcty = list(info.root_dict_list[-1].values())[0]
      print(i)
      print(p_shift)
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


