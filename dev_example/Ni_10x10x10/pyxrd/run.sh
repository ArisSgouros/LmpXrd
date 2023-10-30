#!/bin/bash
path_python="../../../pyxrd/xrd.py"
path_data="../gen/"

# from scratch
python $path_python $path_data/o.ni.dat $path_data/o.ni.lammpstrj 17.6 0.1 20.0 0.05 1 1 -coltype=2 > o.log_rdf

# from restart
python $path_python $path_data/o.ni.dat $path_data/o.ni.lammpstrj 17.6 0.1 20.0 0.05 1 1 -coltype=2 -rdffile=o.AllRDFs.dat > o.log_xrd
