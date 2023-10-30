#!/bin/bash
path_python="../../..//pyxrd/xrd.py"
path_data="../gen/"

# call
#python $path_python

# from scratch
python $path_python $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj 25.0 0.1 20.0 0.05 1 1 -coltype=2 > o.log_rdf

# from restart
python $path_python $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj 25.0 0.1 20.0 0.01 1 1 -coltype=2 -rdffile=o.AllRDFs.dat > o.log_xrd
