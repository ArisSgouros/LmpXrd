import numpy as np
import argparse
from scipy.special import spherical_jn

parser = argparse.ArgumentParser(description='Calculate group scattering factors')
parser.add_argument('group_file', type=str, help='Path of the xyz mol file')
parser.add_argument('formfactorfile', type=str, default="", help='Path of the form factor file')
parser.add_argument('qmax', type=float, help='qmax')
parser.add_argument('dq', type=float, help='dq')
parser.add_argument('-orig', '--orig', type=str, default="0,0,0", help='Origin coordinates.')
parser.add_argument('-scheme', '--scheme', type=str, default="sum", help='Select the weighting scheme (sum/narten).')
parser.add_argument('-ff_group_file', '--ff_group_file', type=str, default="", help='Path of the group form factor as a function of q')

def CartesianToSpherical(r_cart, r_orig):
   xx, yy, zz = (r_cart - r_orig)
   rr = np.sqrt(xx*xx + yy*yy + zz*zz)
   if rr == 0.0:
      return np.zeros(3)
   else:
      th = np.arccos(zz/rr)
      ph = np.arctan2(yy,xx)
      return np.array([rr, th, ph])

def FormFact(p, qq):
   fa = p["c"]
   for ii in range(4):
      fa += p["a"][ii] * np.exp(-p["b"][ii]*pow(qq/(4*np.pi),2))
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

   FFp = {}
   for element in element_types:
      FFp[element] = {"a" : [0.0]*4, "b" : [0.0]*4, "c" : 0.0}

   with open(file_ff_atom,"r") as openFileObject:
      for cur_line in openFileObject:
         cur_line = cur_line.split()
         for element in element_types:
            if element == cur_line[0]:
               FFp[element]["a"][0] = float(cur_line[1])
               FFp[element]["a"][1] = float(cur_line[3])
               FFp[element]["a"][2] = float(cur_line[5])
               FFp[element]["a"][3] = float(cur_line[7])
               FFp[element]["b"][0] = float(cur_line[2])
               FFp[element]["b"][1] = float(cur_line[4])
               FFp[element]["b"][2] = float(cur_line[6])
               FFp[element]["b"][3] = float(cur_line[8])
               FFp[element]["c"]    = float(cur_line[9])

   print( "Atomic Form Factors")

   for element in FFp:
      print( "element :", element)
      print( "a", FFp[element]["a"])
      print( "b", FFp[element]["b"])
      print( "c", FFp[element]["c"])
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
            ff = FormFact(FFp[el], qq)
            jn = spherical_jn(mm, qq*rr)
            a00 += ff*spherical_jn(mm, qq*rr)
            #print("   el: ", el, ", ff: ", ff, ", jn: ", jn, ", rr: ", rr)
         ff_group_list.append(a00)
   elif scheme == "sum":
      for qq in qq_list:
         ff_group = 0.0
         for ii in range(natom):
            el = element_list[ii]
            ff_atom = FormFact(FFp[el], qq)
            ff_group += ff_atom
         ff_group_list.append(ff_group)

   # export the group form factor
   if not file_ff_group == "":
      foo = open(file_ff_group, "w")
      foo.write("%-16s %-16s " % ("q", "ff_group"))
      for el in element_list:
         foo.write("%-16s " % (el))
      foo.write("\n")
      for ii in range(nbin):
         foo.write("%-16.9f %-16.9f " % (qq_list[ii], ff_group_list[ii]))
         for el in element_list:
            foo.write("%-16.9f " % (FormFact(FFp[el], qq_list[ii])))
         foo.write("\n")
      foo.close()
