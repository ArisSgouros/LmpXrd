import os
import sys
import ast
import numpy as np
import math as math

n_command_line_var = len(sys.argv)

if n_command_line_var != 9 and n_command_line_var != 10:
   print()
   print( "The required formats are the following:")
   print()
   print( "python RDF_PER_TYPE.py \"LMP_DATA\" \"LMP_CONF\" \"rmax\" \"dr\" \"qmax\" \"dq\" \"NFRAME\" \"NEVERY\"")
   print( "example: python RDF_PER_TYPE.py 92pc.data IO.CONF.lammpstrj 10 0.25 15 0.1 10 2")
   print()
   print( "OR if you wish to read the partial gr from a previous calculation")
   print()
   print( "python RDF_PER_TYPE.py \"LMP_DATA\" \"LMP_CONF\" \"rmax\" \"dr\" \"qmax\" \"dq\" \"NFRAME\" \"NEVERY\" \"Gab_FILE\"")
   print( "example: python RDF_PER_TYPE.py 92pc.data IO.CONF.lammpstrj 10 0.25 15 0.1 10 2 o.AllRDFs.dat")
   print()
   print( "*The format of the Masses section of the data file should be like the following:")
   print()
   print( " Masses")
   print()
   print( " 1    15.9994   # O")
   print( " 3    12.0107   # C")
   print( " 2    12.0107   # C")
   print( " ..")
   print()
   print( "*The format of the atom section of the data file should be like the following:")
   print( " atomId molId type charge x y z")
   print()
   print( "*The formal of the configuration file should be:")
   print( " id type x y z")
   print()
   print( "exiting..")
   sys.exit()

dataFileName = sys.argv[1]
confFileName = sys.argv[2]
rmax = float(sys.argv[3])
dr = float(sys.argv[4])
qmax = float(sys.argv[5])
dq = float(sys.argv[6])
nFrame = int(sys.argv[7])
nEvery = int(sys.argv[8])

print( "*Parameters of the calculation*")
print()
print( "dataFileName   :",dataFileName)
print( "confFileName   :",confFileName)
print( "rmax           :",rmax)
print( "dr             :",dr)
print( "qmax           :",qmax)
print( "dq             :",dq)
print( "NFrames        :",nFrame)
print( "NEvery         :",nEvery)

computeRDFs = False
readRDFs = False
if n_command_line_var == 10:
   RDFFileName = sys.argv[9]
   readRDFs = True
   print( "\nThe partial gr will be read from:",RDFFileName)
else:
   computeRDFs = True
   print( "\nThe partial gr will be computed from scratch")

#
# Initialize the histograms
#
rbins = int(rmax / dr)
r_range = [(ii+0.5)*dr for ii in range(rbins)]

qmin = 1e-5
qbins = int(qmax / dq)
q_range = [(ii+0.5)*dq for ii in range(qbins)]

#
# Read the DATA file and get the atom types
#
print( "Reading DATAFILE from:",dataFileName,"..")

fileLines = []
iLine = 0


with open(dataFileName,"r") as openFileObject:
   for curLine in openFileObject:

      if "atom types" in curLine:
         nAtomTypes = int(curLine.split()[0])
      if "atoms" in curLine:
         nAtoms = int(curLine.split()[0])
      if "Masses" in curLine:
         typesLineStart = iLine + 2
         typesLineEnd = typesLineStart + nAtomTypes
      if "Atoms" in curLine:
         atomLineStart = iLine + 2
         atomLineEnd = atomLineStart + nAtoms

      fileLines.append(curLine)
      iLine += 1

typesOfSpecies = {}
speciesOfType = {}
for curLine in range(typesLineStart,typesLineEnd):
   curLineSplit = fileLines[curLine].split()
   species = curLineSplit[3]
   type = int(curLineSplit[0])

   if species in typesOfSpecies.keys():
      typesOfSpecies[species].append(type)
   else:
      typesOfSpecies[species] = [type]

   speciesOfType[type] = species

for species in typesOfSpecies:
   print( species, typesOfSpecies[species])

IdsOfSpecies = {}
for species in typesOfSpecies:
   IdsOfSpecies[species] = []

for curLine in range(atomLineStart,atomLineEnd):
   curLineSplit = fileLines[curLine].split()
   ID = int(curLineSplit[0])
   type = int(curLineSplit[2])
   species = speciesOfType[type]
   IdsOfSpecies[species].append(ID)

N = {}
C = {}
N["all"] = 0
C["all"] = 1.0
for species in typesOfSpecies:
   N[species] = len(IdsOfSpecies[species])
   N["all"] += N[species]

for species in typesOfSpecies:
   C[species] = float(N[species]) / float(N["all"])
   print( "There are",N[species]," ",species,"atoms (",C[species],"%)")

print( "There are",N["all"], "atoms in total..")

#Write the atom IDs for each species in a file
print( "Printing the atom IDs for each species..")
for species in typesOfSpecies:
   f = open('o.IDS.'+species+'.dat', 'w')
   f.write(species+"\n")
   for ID in IdsOfSpecies[species]:
      f.write(str(ID)+" ")
f.close()

# Track the computed RDFs so that we wont have to recompute them again
gab_all = {}

#
# The section computes the RDFs from the trajectory files
#

def GetDumpFormat(filename):
   col_id = col_xx = col_yy = col_zz = -1
   f = open(confFileName,"r")
   line = ""
   for ii in range(9):
      line = f.readline()
   f.close()
   header = line.split()
   header = header[2:] # pop first two elements ('ITEM:', 'ATOMS')
   for icol in range(len(header)):
      attrib = header[icol]
      if attrib == 'id':
         col_id = icol
      if attrib in ['x', 'xu', 'xs']:
         col_xx = icol
      if attrib in ['y', 'yu', 'ys']:
         col_yy = icol
      if attrib in ['z', 'zu', 'zs']:
         col_zz = icol
   if -1 in [col_id, col_xx, col_yy, col_zz]:
      print("error with column ids in dump file")
      print("header: ", line)
      print("Id    : ", col_id)
      print("x     : ", col_xx)
      print("y     : ", col_yy)
      print("z     : ", col_zz)
      sys.exit()
   return col_id, col_xx, col_yy, col_zz

# Check the format of the dump file
col_id, col_xx, col_yy, col_zz = GetDumpFormat(confFileName)

if computeRDFs == True:
   # Start iterating over all species
   for species_A in typesOfSpecies:
      for species_B in typesOfSpecies:

         if species_A+"_"+species_B in gab_all.keys() or species_B+"_"+species_A in gab_all.keys():
            continue

         print( "RDF computation of",species_A,"->",species_B,"..")

         # Initialize the array of vectors with dimensions:
         #
         # Load the atom trajectories
         #
         f = open(confFileName,"r")

         RDF = [0] * rbins
         rAB = [0.0] * 3

         nSamples = 0
         Volume = 0.0
         for iFrame in range(nFrame):

            if iFrame % nEvery != 0:
               for ii in range(nAtoms + 9):
                  f.readline()
               continue
            else:
               nSamples += 1

            print( iFrame, "/", nFrame, "frames..")

            f.readline()                # ITEM: TIMESTEP
            Timestep = int(f.readline())

            f.readline()                # ITEM: NUMBER OF ATOMS
            f.readline()
            f.readline()                # ITEM: BOX BOUNDS

            X_BOX = f.readline().split()  # xlo xhi
            xlo = float(X_BOX[0])
            xhi = float(X_BOX[1])
            Y_BOX = f.readline().split()  # xlo xhi
            ylo = float(Y_BOX[0])
            yhi = float(Y_BOX[1])
            Z_BOX = f.readline().split()  # xlo xhi
            zlo = float(Z_BOX[0])
            zhi = float(Z_BOX[1])

            L0 = xhi - xlo
            L1 = yhi - ylo
            L2 = zhi - zlo

            Volume += L0 * L1 * L2

            f.readline()                # ITEM: FORMAT
            #
            # Loop over all atoms and create a
            # dictionary for atom IDs and coordinates
            #
            pos = {}
            for ii in range(nAtoms):
               line = f.readline().split()
               id = int(line[col_id])
               xx = float(line[col_xx])
               yy = float(line[col_yy])
               zz = float(line[col_zz])
               pos.update({id:[xx,yy,zz]})

            for id_A in IdsOfSpecies[species_A]:

               # Start of optimizations
               pos_A0 = pos[id_A][0]
               pos_A1 = pos[id_A][1]
               pos_A2 = pos[id_A][2]
               iL0 = 1.0 / L0
               iL1 = 1.0 / L1
               iL2 = 1.0 / L2
               rmax_sq = rmax * rmax
               i_dr = 1.0 / dr
               # End of optimizations

               for id_B in IdsOfSpecies[species_B]:

                  rAB_0 = pos[id_B][0] - pos_A0
                  rAB_0 -= L0 * int(round( rAB_0 * iL0 ))
                  rAB_1 = pos[id_B][1] - pos_A1
                  rAB_1 -= L1 * int(round( rAB_1 * iL1 ))
                  rAB_2 = pos[id_B][2] - pos_A2
                  rAB_2 -= L2 * int(round( rAB_2 * iL2 ))

                  rAB_sq = rAB_0 * rAB_0 + rAB_1 * rAB_1 + rAB_2 * rAB_2

                  if( rAB_sq > rmax_sq):
                     continue

                  rAB = np.sqrt(rAB_sq)
                  m = int(rAB * i_dr)
                  RDF[m] += 1

         Volume /= nSamples
         RDF[0] = 0

         f.close()

         g = open('o.'+species_A+'-'+species_B+'.dat', 'w')
         for m in range(rbins):

            rho = float(N["all"]) / Volume
            shell_vol = 4.0 / 3.0 * np.pi * (pow(m+1, 3) - pow(m, 3)) * pow(dr, 3)
            nid = shell_vol * rho

            RDF[m] = float(RDF[m]) / float(nSamples)
            RDF[m] *= float(N["all"]) / (float(N[species_A]) * float(N[species_B]) * nid)

            g.write("%10d %10.6f %10.6f %10.6f\n" % (m, r_range[m], RDF[m], shell_vol))
         g.close()

         # Set the completed flags to True in order to skip a similar computationi
         gab_all[species_A+"_"+species_B] = RDF
         gab_all[species_B+"_"+species_A] = RDF

   g = open('o.AllRDFs.dat', 'w')
   g.write("Nall = %-20s\n" % (str(N["all"])))
   g.write("rho = %-.16f\n" % (rho))
   g.write("%-20s" % ("type"))
   for key in gab_all:
      g.write("%-20s" % (key))
   g.write("\n")

   g.write("%-20s" % ("ca"))
   for key in gab_all:
      species_A = key.split("_")[0]
      g.write("%-20s" % (str(C[species_A])))
   g.write("\n")

   g.write("%-20s" % ("cb"))
   for key in gab_all:
      species_B = key.split("_")[1]
      g.write("%-20s" % (str(C[species_B])))
   g.write("\n")

   for ir in range(rbins):
      g.write("%-20.6f" %(r_range[ir]))
      for key in gab_all:
         g.write("%-20.6f" % (gab_all[key][ir]))
      g.write("\n")
   g.close()

#
# The section reads the RDFs made from previous calculations
#
if readRDFs:
   print( "Reading the partial gr from file",RDFFileName,"..")
   f = open(RDFFileName,"r")
   f.readline()
   rho = float(f.readline().split()[2])
   print( rho)
   types = f.readline().split()
   types.pop(0)
   for type in types:
      gab_all[type] = [0] * rbins

   f.readline() # CA
   f.readline() # CB
   for ir in range(rbins):
      currentLine = f.readline().split()
      currentLine.pop(0)
      for it in range(len(types)):
         type = types[it]
         gab_all[type][ir] = float(currentLine[it])

#
# This section deals with the atomic form factors
#
# The structure factors were retrieved from :
# http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

fileName = "FORM_FACTORS.dat"

print( "Reading the atomic form factors from",fileName,"..\n")

FFp = {}
for species in typesOfSpecies:
   FFp[species] = {"a" : [0.0]*4, "b" : [0.0]*4, "c" : 0.0}

with open(fileName,"r") as openFileObject:
   for curLine in openFileObject:
      curLine = curLine.split()
      for species in typesOfSpecies:
         if species == curLine[0]:
            FFp[species]["a"][0] = float(curLine[1])
            FFp[species]["a"][1] = float(curLine[3])
            FFp[species]["a"][2] = float(curLine[5])
            FFp[species]["a"][3] = float(curLine[7])
            FFp[species]["b"][0] = float(curLine[2])
            FFp[species]["b"][1] = float(curLine[4])
            FFp[species]["b"][2] = float(curLine[6])
            FFp[species]["b"][3] = float(curLine[8])
            FFp[species]["c"]    = float(curLine[9])

print( "Atomic Form Factors")

for species in FFp:
   print( "species :", species)
   print( "a", FFp[species]["a"])
   print( "b", FFp[species]["b"])
   print( "c", FFp[species]["c"])
print()

def FormFact(p, qq):
   fa = p["c"]
   for ii in range(4):
      fa += p["a"][ii] * np.exp(-p["b"][ii]*pow(qq/(4*np.pi),2))
   return fa

print( "Computation of the Faber-Ziman structure factor..\n")

SFZ_all = {}
for species_A in typesOfSpecies:
   for species_B in typesOfSpecies:

      SFZ = [0.0] * qbins
      #print( "SFZ computation of",species_A,"->",species_B,"..")

      type = species_A+"_"+species_B
      gab = gab_all[type]
      for iq in range(qbins):
         qq = q_range[iq]
         r_integral = 0.0
         for ir in range(rbins):
            rr = r_range[ir]
            r_integral += dr * pow(rr,2) * math.sin(qq*rr)/(qq*rr) * (gab[ir]-1)

         SFZ[iq] = 1.0 + 4 * np.pi * rho * r_integral
      SFZ_all[type] = SFZ

g = open('o.All_SFZ.dat', 'w')
g.write("%-20s" % ("type"))
for key in gab_all:
   g.write("%-20s" % (key))
g.write("\n")

g.write("%-20s" % ("ca"))
for key in gab_all:
   species_A = key.split("_")[0]
   g.write("%-20s" % (str(C[species_A])))
g.write("\n")

g.write("%-20s" % ("cb"))
for key in gab_all:
   species_B = key.split("_")[1]
   g.write("%-20s" % (str(C[species_B])))
g.write("\n")

for iq in range(qbins):
   g.write("%-20.6f" %(q_range[iq]))
   for key in SFZ_all:
      g.write("%-20.6f" % (SFZ_all[key][iq]))
   g.write("\n")
g.close()

print( "Computation of the X-ray weighted structure factor..\n")

Fx = [0.0] * qbins
for iq in range(qbins):
   qq = q_range[iq]

   denominator = 0.0
   for species_A in typesOfSpecies:
      FA = FormFact(FFp[species_A],qq)
      denominator += C[species_A] * FA
   denominator = pow(denominator,2)

   numerator = 0.0
   for species_A in typesOfSpecies:
      for species_B in typesOfSpecies:
         type = species_A+"_"+species_B
         FA = FormFact(FFp[species_A],qq)
         FB = FormFact(FFp[species_B],qq)

         numerator += FA * FB * C[species_A] * C[species_B] * (SFZ_all[type][iq] - 1.0)

   Fx[iq] = numerator / denominator


print( "Computation of the coherent scattering intensity..\n")

Icoh = [0.0] * qbins
for iq in range(qbins):
   qq = q_range[iq]

   denominator = 0.0
   for species_A in typesOfSpecies:
      FA = FormFact(FFp[species_A],qq)
      denominator += C[species_A] * FA
   denominator = pow(denominator,2)

   sum_c_fsq = 0.0
   for species_A in typesOfSpecies:
      FA = FormFact(FFp[species_A],qq)
      sum_c_fsq += C[species_A] * FA * FA

   Icoh[iq] = Fx[iq] * denominator + sum_c_fsq

g = open('o.Fx_Icoh.dat', 'w')
g.write("%-19s %-19s %-19s %-19s\n" % ("bin", "q", "Fx", "Icoh"))
for iq in range(qbins):
   g.write( "%-19d %-19.6f %-19.6f %-19.6f\n"  %(iq, q_range[iq], Fx[iq], Icoh[iq]))
g.close()
print( "Done!")

#
# Direct calculation of the coherent scattering intensity
#
S_total = [0.0] * qbins

for iq in range(qbins):
   qq = q_range[iq]

   term_B = 0.0

   for species_A in typesOfSpecies:
      for species_B in typesOfSpecies:

         r_integral = 0.0
         gab = gab_all[species_A+"_"+species_B]
         for ir in range(rbins):
            rr = r_range[ir]
            r_integral += dr * pow(rr,2) * math.sin(qq*rr)/(qq*rr) * (gab[ir]-1)

         FA = FormFact(FFp[species_A],qq)
         FB = FormFact(FFp[species_B],qq)

         sum = 4 * np.pi * rho * FA * FB * C[species_A] * C[species_B] * r_integral

         term_B += sum

   term_A = 0.0
   for species_A in typesOfSpecies:
      FA = FormFact(FFp[species_A],qq)
      term_A += C[species_A] * FA * FA

   S_total[iq] = term_B + term_A

g = open('o.S_XRD.dat', 'w')
g.write("%-19s %-19s %-19s\n" % ("bin", "q", "S"))
for iq in range(qbins):
   g.write( "%-19d %-19.6f %-19.6f\n"  %(iq, q_range[iq], S_total[iq]))
g.close()
print( "Done!")
