#!/bin/bash
#path_python="/home/aps/Programs/XrdDev/pyxrd/xrd.py"
path_python="xrd.py"
path_data="../gen/"

# simple call
#python $path_python

# from scratch 1 frame
#python $path_python $path_data/pos_formatted.dat $path_data/movie.lammpstrj 6.5 0.065 20.0 0.05 1 1

# from scratch
python $path_python $path_data/pos_formatted.dat $path_data/movie.lammpstrj 6.5 0.065 20.0 0.05 10 1

# from restart
#python $path_python $path_data/o.ni.dat $path_data/o.ni.lammpstrj 17.6 0.1 20.0 0.05 1 1 o.AllRDFs.dat
