import os
import sys
import ast
import numpy as np
import math as math
import argparse

kTol = 1e-6

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('datafile', type=str, help='Path of the lammps data file')
parser.add_argument('dumpfile', type=str, help='Path of the lammps dump file')
parser.add_argument('rmax', type=float, help='Max radious')
parser.add_argument('dr', type=float, help='Step interval')
parser.add_argument('nframe', type=int, help='number of frames')
parser.add_argument('nevery', type=int, help='read frame every')
parser.add_argument('-coltype', type=int, default=2, help='Atom type column; Full=2, Atomic=1')

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

if __name__ == "__main__":
   args = parser.parse_args()
   dataFileName = args.datafile
   confFileName = args.dumpfile
   rmax         = args.rmax
   dr           = args.dr
   nFrame       = args.nframe
   nEvery       = args.nevery
   col_type     = args.coltype

   print( "*Parameters of the calculation*")
   print()
   print( "dataFileName   :",dataFileName)
   print( "confFileName   :",confFileName)
   print( "rmax           :",rmax)
   print( "dr             :",dr)
   print( "NFrames        :",nFrame)
   print( "NEvery         :",nEvery)

   #
   # Initialize the histograms
   #
   rbins = int(rmax / dr)
   r_range = [(ii+0.5)*dr for ii in range(rbins)]

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
   # The section computes the RDFs from the trajectory files
   #

   # Check the format of the dump file
   col_id, col_xx, col_yy, col_zz = GetDumpFormat(confFileName)

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
               rmax_sq = rmax * rmax - kTol
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

                  if( rAB_sq >= rmax_sq):
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
