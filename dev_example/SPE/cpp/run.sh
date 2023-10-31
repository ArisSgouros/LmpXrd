#!/bin/bash
path_gij="../../../cpp/gij_partial"

# compute intra (intra)
$path_gij 15 0.5 1 INTRA ../CG.dat ../CG.lammpstrj gij_pairs.dat > o.log_rdf_intra
mv o.AllRDFs.dat o.rdf_intra.dat

# compute inter (inter)
$path_gij 15 0.5 1 INTER ../CG.dat ../CG.lammpstrj gij_pairs.dat > o.log_rdf_inter
mv o.AllRDFs.dat o.rdf_inter.dat

# compute rdf (full)
$path_gij 15 0.5 1 FULL  ../CG.dat ../CG.lammpstrj gij_pairs.dat > o.log_rdf_full
mv o.AllRDFs.dat o.rdf_full.dat
