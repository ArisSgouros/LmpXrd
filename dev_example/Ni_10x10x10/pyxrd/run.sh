#!/bin/bash
path_python="../../../pyxrd/"
path_data="../gen/"

# compute rdf
python $path_python/pyrdf.py $path_data/o.ni.dat $path_data/o.ni.lammpstrj 17.6 0.1 1 1 -atomtype='full' -rdffile='o.rdf' > o.log_rdf

# compute xrd
python $path_python/pyxrd.py $path_data/o.ni.dat o.rdf.dat form_factors.dat 17.6 0.1 20.0 0.05 -atomtype='full' > o.log_xrd
