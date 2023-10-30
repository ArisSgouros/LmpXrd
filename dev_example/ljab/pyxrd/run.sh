#!/bin/bash
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
python $path_python/pyrdf.py $path_data/pos_formatted.dat $path_data/movie.lammpstrj 6.5 0.065 10 1 -coltype=1 > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/pos_formatted.dat o.AllRDFs.dat 6.5 0.065 20.0 0.05 -coltype=1 > o.log_xrd
