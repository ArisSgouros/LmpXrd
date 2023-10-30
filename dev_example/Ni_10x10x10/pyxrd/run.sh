#!/bin/bash
path_python="../../../pyxrd/"
path_data="../gen/"

# from scratch
python $path_python/pyrdf.py $path_data/o.ni.dat $path_data/o.ni.lammpstrj 17.6 0.1 1 1 -coltype=2 > o.log_rdf

# from restart
python $path_python/pyxrd.py $path_data/o.ni.dat o.AllRDFs.dat 17.6 0.1 20.0 0.05 -coltype=2 > o.log_xrd
