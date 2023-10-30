#!/bin/bash
path_python="/home/aps/Programs/XrdDev/pyxrd/xrd.py"
path_data="../gen/"

#
#python $path_python

# from scratch
#python $path_python $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj 25.0 0.1 20.0 0.05 1 1

# from restart
python $path_python $path_data/o.nacl.dat $path_data/o.nacl.lammpstrj 25.0 0.1 20.0 0.01 1 1 o.AllRDFs.dat
