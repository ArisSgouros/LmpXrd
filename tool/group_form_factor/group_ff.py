import numpy as np
import argparse
from scipy.special import spherical_jn
from scipy.optimize import minimize

parser = argparse.ArgumentParser(description='Calculate group scattering factors')
parser.add_argument('group_file', type=str, help='Path of the xyz mol file')
parser.add_argument('formfactorfile', type=str, default="", help='Path of the form factor file')
parser.add_argument('qmax', type=float, help='qmax')
parser.add_argument('dq', type=float, help='dq')
parser.add_argument('-orig', '--orig', type=str, default="0,0,0", help='Origin coordinates.')
parser.add_argument('-scheme', '--scheme', type=str, default="sum", help='Select the weighting scheme (sum/narten).')
parser.add_argument('-ff_group_file', '--ff_group_file', type=str, default="", help='Path of the group form factor as a function of q.')
parser.add_argument('-init_guess', '--init_guess', type=str, default="1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0", help='Initial guess of group coeffs.')

def CartesianToSpherical(r_cart, r_orig):
   xx, yy, zz = (r_cart - r_orig)
   rr = np.sqrt(xx*xx + yy*yy + zz*zz)
   if rr == 0.0:
      return np.zeros(3)
   else:
      th = np.arccos(zz/rr)
      ph = np.arctan2(yy,xx)
      return np.array([rr, th, ph])

class FormFactCoeff:
   def __init__(self):
      self.a = [0.0]*4
      self.b = [0.0]*4
      self.c = 0.0

def FormFact(coeff, qq):
   fa = coeff.c
   for ii in range(4):
      fa += coeff.a[ii] * np.exp(-coeff.b[ii]*pow(qq/(4*np.pi),2))
   return fa

if __name__ == "__main__":
   args = parser.parse_args()
   file_group = args.group_file
   file_ff_atom   = args.formfactorfile
   qmax = args.qmax
   dq = args.dq
   scheme = args.scheme
   r_orig = np.array([float(item) for item in args.orig.split(',')])
   file_ff_group = args.ff_group_file

   print("group file : ", file_group)
   print("ff file    : ", file_ff_atom)
   if scheme not in ["narten", "sum"]:
      print("Unknown weighting scheme " + str(scheme))
      sys.exit()
   print("scheme     : ", scheme)
   print("origin     : ", r_orig)
   print("qmax       : ", qmax)
   print("dq         : ", dq)

   #             a1   b1   a2   b2   a3   b3   a4   b4   c
   init_guess = [float(item) for item in args.init_guess.split(',')]

   nbin = int(qmax / dq)
   qq_list = [(ii+0.5)*dq for ii in range(nbin)]

   natom = 0
   rcart_list = []
   element_list = []
   element_types = set()

   # parse group
   with open(file_group) as foo:
      line_split = foo.readline()
      natom = int(line_split[0])
      foo.readline()
      for iatom in range(natom):
         line_split = foo.readline().split()
         el, x, y, z = line_split
         element_list.append(el)
         element_types.add(el)
         rcart_list.append(np.array([float(x), float(y), float(z)]))

   #
   # This section deals with the atomic form factors
   #
   # The structure factors were retrieved from :
   # http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
   print( "Reading the atomic form factors from",file_ff_atom,"..\n")

   form_fact_coeff_el = {}

   for element in element_types:
      form_fact_coeff_el[element] = FormFactCoeff()

   with open(file_ff_atom,"r") as openFileObject:
      for cur_line in openFileObject:
         cur_line = cur_line.split()
         for element in element_types:
            if element == cur_line[0]:
               cur_line.pop(0)
               coeff_list = [float(item) for item in cur_line]
               form_fact_coeff_el[element].a = [coeff_list[0], coeff_list[2], coeff_list[4], coeff_list[6]]
               form_fact_coeff_el[element].b = [coeff_list[1], coeff_list[3], coeff_list[5], coeff_list[7]]
               form_fact_coeff_el[element].c = coeff_list[8]

   print( "Atomic Form Factors")
   for element in form_fact_coeff_el:
      print( "element :", element)
      print( "a", form_fact_coeff_el[element].a)
      print( "b", form_fact_coeff_el[element].b)
      print( "c", form_fact_coeff_el[element].c)
   print()

   rsph_list = []
   for rcart in rcart_list:
      rsph_list.append(CartesianToSpherical(rcart, r_orig))

   # compute the molecular scattering factor for m=0, n=0 [eq. 2 in https://doi.org/10.1063/1.435096]
   ff_group_list = []

   if scheme == "narten":
      mm = 0
      nn = 0
      for qq in qq_list:
         a00 = 0.0
         for ii in range(natom):
            el = element_list[ii]
            [rr, th, phi] = rsph_list[ii]
            ff = FormFact(form_fact_coeff_el[el], qq)
            jn = spherical_jn(mm, qq*rr)
            a00 += ff*spherical_jn(mm, qq*rr)
            #print("   el: ", el, ", ff: ", ff, ", jn: ", jn, ", rr: ", rr)
         ff_group_list.append(a00)
   elif scheme == "sum":
      for qq in qq_list:
         ff_group = 0.0
         for ii in range(natom):
            el = element_list[ii]
            ff_atom = FormFact(form_fact_coeff_el[el], qq)
            ff_group += ff_atom
         ff_group_list.append(ff_group)

   #
   # fit effective form factor coefficients for group
   #
   def ObjFunc(coeff_list, qq_list, ff_group_list):
      # unpack the coefficients
      ff_g = FormFactCoeff()
      ff_g.a = [coeff_list[0], coeff_list[2], coeff_list[4], coeff_list[6]]
      ff_g.b = [coeff_list[1], coeff_list[3], coeff_list[5], coeff_list[7]]
      ff_g.c = coeff_list[8]

      merit = 0.0
      n_eval = len(qq_list)
      for ii in range(n_eval):
         qq = qq_list[ii]
         ff_model = FormFact(ff_g, qq)
         merit += pow(ff_model - ff_group_list[ii], 2)
      merit /= n_eval
      #print(merit)
      return merit

   res = minimize(ObjFunc, init_guess, args=(qq_list, ff_group_list),  method='nelder-mead',options={'xtol': 1e-9, 'disp': True})
   res = minimize(ObjFunc, res.x, args=(qq_list, ff_group_list),  method='BFGS',options={'disp': True})

   ff_g = FormFactCoeff()
   ff_g.a = [res.x[0], res.x[2], res.x[4], res.x[6]]
   ff_g.b = [res.x[1], res.x[3], res.x[5], res.x[7]]
   ff_g.c = res.x[8]

   print("Fitting completed.")
   print("   error:")
   print("      ", ObjFunc(res.x, qq_list, ff_group_list))
   print("   coeffs:")
   print("      %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s %-16s" % ("a1","b1","a2","b2","a3","b3","a4","b4","c"))
   print("      %-16.9f %-16.9f %-16.9f %-16.9f %-16.9f %-16.9f %-16.9f %-16.9f %-16.9f" % \
                (ff_g.a[0], ff_g.b[0], ff_g.a[1], ff_g.b[1], ff_g.a[2], ff_g.b[2], ff_g.a[3], ff_g.b[3], ff_g.c))

   # export the group form factor
   if not file_ff_group == "":
      foo = open(file_ff_group, "w")
      foo.write("%-16s %-16s %-16s " % ("q", "ff_group_dir", "ff_group_fit"))
      for el in element_list:
         foo.write("%-16s " % (el))
      foo.write("\n")
      for ii in range(nbin):
         qq = qq_list[ii]
         ff_group_dir = ff_group_list[ii]
         ff_group_fit = FormFact(ff_g, qq)
         foo.write("%-16.9f %-16.9f %-16.9f " % (qq, ff_group_dir, ff_group_fit))
         for el in element_list:
            foo.write("%-16.9f " % (FormFact(form_fact_coeff_el[el], qq_list[ii])))
         foo.write("\n")
      foo.close()
