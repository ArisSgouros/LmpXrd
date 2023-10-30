#!/bin/bash
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
python $path_python/pyrdf.py $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj 25.0 0.1 1 1 -coltype=2 > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/o.nacl.dat o.AllRDFs.dat 25.0 0.1 20.0 0.01 -coltype=2 > o.log_xrd
