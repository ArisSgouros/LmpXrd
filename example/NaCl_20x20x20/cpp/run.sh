#!/bin/bash
path_gij="../../../cpp/gij_partial"
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
$path_gij 25.0 0.1 1 FULL $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj gij_pairs.dat full o.rdf > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/o.nacl.dat o.rdf.dat ../form_factors.dat 25.0 0.1 20.0 0.01 -atomtype='full' > o.log_xrd
