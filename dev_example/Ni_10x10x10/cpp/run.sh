#!/bin/bash
path_gij="../../../cpp/gij_partial"
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
$path_gij 17.6 0.1 1 FULL $path_data/o.ni.dat $path_data/o.ni.lammpstrj gij_pairs.dat fmt.dat o.rdf > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/o.ni.dat o.rdf.dat ../form_factors.dat 17.6 0.1 20.0 0.05 -atomtype='full' > o.log_xrd
