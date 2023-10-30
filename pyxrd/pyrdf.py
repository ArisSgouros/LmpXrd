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
parser.add_argument('-rdffile', type=str, default='o.AllRDFs.dat', help='Path of the rdf file')

def GetDumpFormat(filename):
   col_id = col_xx = col_yy = col_zz = -1
   f = open(filename,"r")
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
   file_data = args.datafile
   file_dump = args.dumpfile
   rmax         = args.rmax
   dr           = args.dr
   nframe       = args.nframe
   nevery       = args.nevery
   col_type     = args.coltype
   file_rdf     = args.rdffile

   print( "*Parameters of the calculation*")
   print()
   print( "file_data   :",file_data)
   print( "file_dump   :",file_dump)
   print( "file_rdf    :",file_rdf)
   print( "rmax        :",rmax)
   print( "dr          :",dr)
   print( "nframe      :",nframe)
   print( "nevery      :",nevery)

   # Initialize the binning
   rbins = int(rmax / dr)
   r_range = [(ii+0.5)*dr for ii in range(rbins)]

   # Fetch the atom types from the Lammps data file
   print( "Reading Lammps datafile from:",file_data,"..")
   lines = []
   iline = 0
   with open(file_data,"r") as openFileObject:
      for cur_line in openFileObject:
         if "atom types" in cur_line:
            n_atom_type = int(cur_line.split()[0])
         if "atoms" in cur_line:
            n_atom = int(cur_line.split()[0])
         if "Masses" in cur_line:
            line_start_type = iline + 2
            line_end_type = line_start_type + n_atom_type
         if "Atoms" in cur_line:
            line_start_atom = iline + 2
            line_end_atom = line_start_atom + n_atom
         lines.append(cur_line)
         iline += 1

   type_of_species = {}
   species_of_type = {}
   for cur_line in range(line_start_type,line_end_type):
      cur_lineSplit = lines[cur_line].split()
      species = cur_lineSplit[3]
      type = int(cur_lineSplit[0])
      if species in type_of_species.keys():
         type_of_species[species].append(type)
      else:
         type_of_species[species] = [type]
      species_of_type[type] = species

   for species in type_of_species:
      print( species, type_of_species[species])

   ids_of_species = {}
   for species in type_of_species:
      ids_of_species[species] = []

   for cur_line in range(line_start_atom,line_end_atom):
      cur_lineSplit = lines[cur_line].split()
      id = int(cur_lineSplit[0])
      type = int(cur_lineSplit[col_type])
      species = species_of_type[type]
      ids_of_species[species].append(id)

   N = {}
   C = {}
   N["all"] = 0
   C["all"] = 1.0
   for species in type_of_species:
      N[species] = len(ids_of_species[species])
      N["all"] += N[species]

   for species in type_of_species:
      C[species] = float(N[species]) / float(N["all"])
      print( "There are",N[species]," ",species,"atoms (",C[species],"%)")

   print( "There are",N["all"], "atoms in total..")
   print( "Printing the atom IDs for each species..")
   for species in type_of_species:
      f = open('o.IDS.'+species+'.dat', 'w')
      f.write(species+"\n")
      for id in ids_of_species[species]:
         f.write(str(id)+" ")
   f.close()

   # Compute the rdf
   gab_all = {}

   # Check the format of the dump file
   col_id, col_xx, col_yy, col_zz = GetDumpFormat(file_dump)

   # Start iterating over all species
   for species_a in type_of_species:
      for species_b in type_of_species:

         if species_a+"_"+species_b in gab_all.keys() or species_b+"_"+species_a in gab_all.keys():
            continue

         print( "RDF computation of",species_a,"->",species_b,"..")

         # load the atom trajectories
         f = open(file_dump,"r")

         rdf = [0] * rbins
         rAB = [0.0] * 3
         n_sample = 0
         volume = 0.0
         for iFrame in range(nframe):
            if iFrame % nevery != 0:
               for ii in range(n_atom + 9):
                  f.readline()
               continue

            n_sample += 1

            print( iFrame, "/", nframe, "frames..")

            f.readline()                # ITEM: TIMESTEP
            f.readline()
            f.readline()                # ITEM: NUMBER OF ATOMS
            f.readline()
            f.readline()                # ITEM: BOX BOUNDS

            box_x = f.readline().split()  # xlo xhi
            xlo = float(box_x[0])
            xhi = float(box_x[1])
            box_y = f.readline().split()  # xlo xhi
            ylo = float(box_y[0])
            yhi = float(box_y[1])
            box_z = f.readline().split()  # xlo xhi
            zlo = float(box_z[0])
            zhi = float(box_z[1])

            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo

            volume += lx*ly*lz

            f.readline()                # ITEM: FORMAT
            #
            # Loop over all atoms and create a
            # dictionary for atom IDs and coordinates
            #
            pos = {}
            for ii in range(n_atom):
               line = f.readline().split()
               id = int(line[col_id])
               xx = float(line[col_xx])
               yy = float(line[col_yy])
               zz = float(line[col_zz])
               pos.update({id:[xx,yy,zz]})

            for id_A in ids_of_species[species_a]:
               pos_A0 = pos[id_A][0]
               pos_A1 = pos[id_A][1]
               pos_A2 = pos[id_A][2]
               ilx = 1.0 / lx
               ily = 1.0 / ly
               ilz = 1.0 / lz
               rmax_sq = rmax * rmax - kTol
               i_dr = 1.0 / dr
               for id_B in ids_of_species[species_b]:
                  rAB_0 = pos[id_B][0] - pos_A0
                  rAB_0 -= lx * int(round( rAB_0 * ilx ))
                  rAB_1 = pos[id_B][1] - pos_A1
                  rAB_1 -= ly * int(round( rAB_1 * ily ))
                  rAB_2 = pos[id_B][2] - pos_A2
                  rAB_2 -= lz * int(round( rAB_2 * ilz ))
                  rAB_sq = rAB_0 * rAB_0 + rAB_1 * rAB_1 + rAB_2 * rAB_2

                  if rAB_sq >= rmax_sq: continue

                  rAB = np.sqrt(rAB_sq)
                  m = int(rAB*i_dr)
                  rdf[m] += 1

         volume /= n_sample
         rdf[0] = 0

         f.close()

         g = open('o.'+species_a+'-'+species_b+'.dat', 'w')
         for m in range(rbins):
            rho = float(N["all"]) / volume
            shell_vol = 4.0 / 3.0 * np.pi * (pow(m+1, 3) - pow(m, 3)) * pow(dr, 3)
            nid = shell_vol * rho
            rdf[m] = float(rdf[m]) / float(n_sample)
            rdf[m] *= float(N["all"]) / (float(N[species_a]) * float(N[species_b]) * nid)

            g.write("%10d %10.6f %10.6f %10.6f\n" % (m, r_range[m], rdf[m], shell_vol))
         g.close()

         # Set the completed flags to True in order to skip a similar computationi
         gab_all[species_a+"_"+species_b] = rdf
         gab_all[species_b+"_"+species_a] = rdf

   g = open(file_rdf, 'w')
   g.write("Nall = %-20s\n" % (str(N["all"])))
   g.write("rho = %-.16f\n" % (rho))
   g.write("%-20s" % ("type"))
   for key in gab_all:
      g.write("%-20s" % (key))
   g.write("\n")

   g.write("%-20s" % ("ca"))
   for key in gab_all:
      species_a = key.split("_")[0]
      g.write("%-20s" % (str(C[species_a])))
   g.write("\n")

   g.write("%-20s" % ("cb"))
   for key in gab_all:
      species_b = key.split("_")[1]
      g.write("%-20s" % (str(C[species_b])))
   g.write("\n")

   for ir in range(rbins):
      g.write("%-20.6f" %(r_range[ir]))
      for key in gab_all:
         g.write("%-20.6f" % (gab_all[key][ir]))
      g.write("\n")
   g.close()
