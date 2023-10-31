#!/bin/bash
path_gij="../../../cpp/gij_partial"

# compute intra (intra)
$path_gij 15 0.5 1 INTRA ../CG.dat ../CG.lammpstrj gij_pairs.dat fmt.dat o.rdf_intra > o.log_rdf_intra

# compute inter (inter)
$path_gij 15 0.5 1 INTER ../CG.dat ../CG.lammpstrj gij_pairs.dat fmt.dat o.rdf_inter > o.log_rdf_inter

# compute rdf (full)
$path_gij 15 0.5 1 FULL  ../CG.dat ../CG.lammpstrj gij_pairs.dat fmt.dat o.rdf_full > o.log_rdf_full
