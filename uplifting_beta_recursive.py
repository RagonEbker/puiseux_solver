from sympy import Pow
import numpy as np
from sympy.polys.ring_series import ring,PolyElement,rs_series_inversion,rs_series,_rs_series,rs_trunc
from sympy import RR, parse_expr,Matrix,simplify,Mod,Add,Mul,Pow,real_roots,symbols,Poly,roots,nan,re,eye,Identity,series,solve,Rational,div,degree,integer_nthroot
from sympy.simplify.simplify import nsimplify

from sympy.core.numbers import One,sympify,Float,S

import scipy.spatial
import cProfile
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
# K = GF(5)

class Info:
  #Precision d
  #x_lifts
  #all the info about the roots in root_dict_list
  #alpha is our first root
  def __init__(self, d, x_lifts,root_dict_list,alpha):
    self.d = d
    self.x_lifts = x_lifts
    self.root_dict_list = root_dict_list
    self.alpha = alpha
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
   
def shift_vertically(p,mtp):   
   A = Matrix.companion(p)
   A = A * mtp
   p = A.charpoly(y)
   return p

#Checks if our Polynomial consist of exactly one root
def check_one_root(p):
   lift = (p.coeffs()[-1]).as_coeff_exponent(x)[1]**(1/2)
   p_c = shift_vertically(p,x**(-lift))
   constant_term = Poly(p_c.subs(x,0),y).coeffs()[-1]
   r_constant = integer_nthroot(constant_term,p_c.degree())
   if (r_constant[1]): 
      r_dict = dict({r_constant[0] * x**(lift) : p_c.degree()})     
   return r_constant[1],r_dict

##Returns whole dict of root, smallest root that is not zero, if it exists, otherwise 0
##and mtpcty of that root
def get_sub_x_root(p):
   p_x = Poly(p.subs(x,0),y)
   p_r_old = roots(p_x)
   p_r = dict()
   for key, value in p_r_old.items():
    # do something with value
    p_r[key] = value
   p_r_list = np.asarray(list(p_r))
   
   #smallest root if not zero
   if (not p_r_list.any()):
      s_r = 0
      mtpcty = p_r[0]
   else:
      s_r = np.min(p_r_list[np.nonzero(p_r_list)])
      #s_r = p_r_list[0]
      mtpcty = p_r[s_r]
   return p_r,s_r,mtpcty

def shift_horizontally(p,shift):
   A = Matrix.companion(p)
   A = A - shift*eye(A.shape[0])
   p_shift = A.charpoly(y)
   print("We shift horizontally by: " + str(shift))
   return p_shift

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
            test_mtpcty = check_one_root(new_poly)
            check_test_mtpcty = test_mtpcty[0]
         except:
            pass
         if (check_test_mtpcty):
            r = test_mtpcty[1]
            r_number = list(r)[0]
            info.alpha += r_number
            info.root_dict_list.append(r)
         else:
            info = calculate_smallest_root_q_x(new_poly,info)
   else:
      r = roots(new_poly)
      r_number = list(r)[0]
      info.alpha += r_number
      info.root_dict_list.append(r)
   return info


def calculate_smallest_root_q_x(p,info):

   #TOODO ADD ROOT DICT TO GLOBAL OBJECT
   root_dict,shift_number = get_sub_x_root(p)
   info.root_dict_list.append(root_dict)
   #lift the polynomial horizontally when all its roots are zeroes
   if (shift_number == 0):
      
      slopes = get_newton_slopes(p)
      min_slope = min(np.abs(slopes))
      p_shift = shift_vertically(p,x**(-min_slope)) 
      info.x_lifts += min_slope
      info.d = info.d-1
      shift_number = calculate_smallest_root_q_x(p_shift,info)
      return info
   info.alpha += shift_number*x**info.x_lifts
   print(info.alpha)
   p_shift = p
   return info

#Calculates the "smallest" root of a polynomial over the Puiseux field
def calculate_smallest_root(p,info):
   
   #TOODO ADD ROOT DICT TO GLOBAL OBJECT
   root_dict,shift_number,mtpcty = get_sub_x_root(p)
   info.root_dict_list.append(root_dict)
   #lift the polynomial horizontally when all its roots are zeroes
   if (shift_number == 0):
      
      slopes = get_newton_slopes(p)
      min_slope = min(np.abs(slopes))
      p_shift = shift_vertically(p,x**(-min_slope)) 
      info.x_lifts += min_slope
      info.d = info.d-1
      shift_number = calculate_smallest_root(p_shift,info)
      
   info.alpha += shift_number*x**info.x_lifts
   print(info.alpha)
   p_shift = p
   for i in range(1,info.d):
      info.d = info.d-1
      p_shift = shift_horizontally(p_shift,shift_number)
      info = golden_lifting(p_shift,mtpcty,info)
      last_root = info.root_dict_list[-1]
      shift_number = list(last_root)[0]*x**info.x_lifts
      print(info.alpha)
      mtpcty = list(last_root.values())[0]
      if (p(intify(info.alpha)) == 0):
         return info
   return info

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
         if (k.as_base_exp()[0] == 1):
               newton_points.append([0, v.atoms().pop().coeff(x)])
         else:
               newton_points.append([k.as_base_exp()[1], v.atoms().pop().coeff(x)])

   points = np.asarray(newton_points)
   hull = scipy.spatial.ConvexHull(points)

   lower_curve = get_lower(points[hull.vertices])
   slopes = [slope(lower_curve[i],lower_curve[i+1]) for i in range(len(lower_curve)-1)]
   return slopes
  


##TODO EXAMPLES
##Examples that have the properties
##for each combination we want an example
#different exponent, same exponent, different coefficient with same exponent, same coefficient with same exponent
#integer exponents, rational exponents
#Finite Fields, ZZ,QQ,RR,CC
#Different precisions
#textbook examples
#competitve examples (with times)

# p_t = Poly(1,y)
# n = 3
#Beispiel 2
#**0.5
# for i in range(1,n):
#    lf = (y- (5*x**(1/2) + 3*x**2))
#    p_t*= lf
#    print(lf)

# n = 2
# ##**06
# for i in range(1,n):
#    lf = (y- (3))
#    p_t*= lf
#    print(lf)

# p_t *= (y-(5*x**(1/2) + 2*x**2))

#p_t = Poly(y**6 + x**6 + 3*x**2*y**4+3*x**4*y**2 - 4*x**2*y**2,y)
#p_t = Poly(4*y**3 + 4*x*y**2+x**2*y+2*x**4,y)



# p_t = Poly(x**5 + 8*x**4-2*x**2*y**2 - y**3+2*y**4,y)
##TODO!!
p_t = Poly((y-(1+3*x**(0.5)+x**2))*(y-(1+3*x+3*x**2))*(y-(4 + 3*x)),y)
p_t *= p_t


#For this one 
#roots(shift_multiply(new_poly,x**(-0.5)))
#With s multiplicity ICH GLAUBE HIER WAR EINS MIT selbem exponenten überall und unterschiedlichen coefficients


#DAS HIER IST DAS WICHTIGSTE BSP bisher 
#roots(shift_multiply(new_poly,x**(-0.4)).subs(x,0))
#p_t = Poly((y-(1+2*x**(0.5)+x**2))*(y-(1+4*x**(0.4)+x**2))*(y-(1+3*x**(0.5)+3*x**2))*(y-(4 + 3*x)),y)

#roots(shift_multiply(new_poly,x**(-0.5)).subs(x,0))
#functioniert auch
#IMMER DEN NIEDRIGSTEN EXPO NEHMEN
#p_t = Poly((y-(1+2*x**(0.4)+x**2))*(y-(1+3*x**(0.5)+x**2))*(y-(1+3*x**(0.6)+3*x**2))*(y-(4 + 3*x)),y)



#p_t = Poly((y-(1+3*x**3+x**4))*(y-(4 + 3*x)),y)

#p_t = Poly((y-(1+2*x**(2)+x**2))*(y-(1+3*x**(1)+x**2))*(y-(1+3*x**(2)+3*x**2))*(y-(4 + 3*x)),y)


#p_t = Poly((y-(x+x**(1/2)))*(y-x),y)
#p_t = Poly((y-(x**(1/2))),y,domain="ZZ[x**0.5]")
#p_t = Poly()
#p_t = Poly((y-(a+b*x)*(y-(c+d*x))*(y-(f+g*x))),y)
#Precision
d = 5
x_lifts = 0
root_dict_list = []
alpha = 0
info = Info(d,x_lifts,root_dict_list,alpha)
#cProfile.run('calculate_smallest_root(p_t,info)',"restats")

info = calculate_smallest_root(p_t,info)
print("")
print("The first root is: " + str(info.alpha))

#p*= (y-(2 + x))
##TODO LOOK AT THIS POLY IT HAS SMALLER THAN X^1 LOWEST B_1_min
##Could be possible that this is not repearable and we have to work with convex hull
#p_t = Poly(x**5 + 8*x**4-2*x**2*y**2 - y**3+2*y**4,y)

#x shift test
#p = Poly((y - (x + x**2))*(y-(x + x**3)),y)
#p = shift_multiply(p,x**(-1))




##CORRECT THE THING WITH THE COMPANION MATRIX
#IT MAY CAUSE ROUNDING ERRORS

##MULTIPLICITIES!!!!

#This should print a power series with a valuation 
#higher than 0.5 and probably even higher than that
#print(p(a+shift_number))
print("")