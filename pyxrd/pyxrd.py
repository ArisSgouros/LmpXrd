import os
import sys
import ast
import numpy as np
import math as math
import argparse

kTol = 1e-6

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('datafile', type=str, help='Path of the lammps data file')
parser.add_argument('rdffile', type=str, default="", help='Path of the RDF file')
parser.add_argument('rmax', type=float, help='Max radious')
parser.add_argument('dr', type=float, help='Step interval')
parser.add_argument('qmax', type=float, help='Max wavevector')
parser.add_argument('dq', type=float, help='wavevector interval')
parser.add_argument('-coltype', type=int, default=2, help='Atom type column; Full=2, Atomic=1')

args = parser.parse_args()
dataFileName = args.datafile
RDFFileName  = args.rdffile
rmax         = args.rmax
dr           = args.dr
qmax         = args.qmax
dq           = args.dq
col_type     = args.coltype

print( "*Parameters of the calculation*")
print()
print( "dataFileName   :",dataFileName)
print( "RDFFileName    :",RDFFileName)
print( "rmax           :",rmax)
print( "dr             :",dr)
print( "qmax           :",qmax)
print( "dq             :",dq)

def FormFact(p, qq):
   fa = p["c"]
   for ii in range(4):
      fa += p["a"][ii] * np.exp(-p["b"][ii]*pow(qq/(4*np.pi),2))
   return fa

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
   type = int(curLineSplit[col_type])
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
# The section reads the RDFs made from previous calculations
#
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
