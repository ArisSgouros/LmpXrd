#!/bin/bash
path_gij="../../../cpp/gij_partial"
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
$path_gij 6.5 0.065 1 FULL $path_data/pos_formatted.dat $path_data/movie.lammpstrj gij_pairs.dat fmt.dat o.rdf > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/pos_formatted.dat o.rdf.dat ../form_factors.dat 6.5 0.065 20.0 0.05 -atomtype='atomic' > o.log_xrd
