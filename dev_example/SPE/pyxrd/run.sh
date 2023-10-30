#!/bin/bash
path_python="../../../pyxrd/"

# compute rdf (full)
python $path_python/pyrdf.py ../CG.dat ../CG.lammpstrj 15.0 0.5 20 1 -atomtype='full' -exclude_interaction='W_W' -neigh='full' -rdffile='o.rdf_full' > o.log_rdf_full

# compute intra (intra)
python $path_python/pyrdf.py ../CG.dat ../CG.lammpstrj 15.0 0.5 20 1 -atomtype='full' -exclude_interaction='W_W' -neigh='intra' -rdffile='o.rdf_intra' > o.log_rdf_intra

# compute inter (inter)
python $path_python/pyrdf.py ../CG.dat ../CG.lammpstrj 15.0 0.5 20 1 -atomtype='full' -exclude_interaction='W_W' -neigh='inter' -rdffile='o.rdf_inter' > o.log_rdf_inter
