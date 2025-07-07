#!/bin/bash

path_python="../../"
path_crystal_builder="/home/asgouros/Programs/CrystalBuilder/crystal_builder.py"

nn=25
r_max=70.0
dr=0.005
q_max=6.0
dq=0.01
path_gij="../../cpp/gij_partial"
path_python="../../pyxrd/"

system_name="o.nacl_$nn"
rdf_name="$system_name"_"$r_max"_"$dr"
xrd_name="$rdf_name"_"$q_max"_"$dq"

# generate
time python $path_crystal_builder in.basis_nacl.dat "$nn,$nn,$nn" --file_pos="$system_name.dat" --file_dump="$system_name.lammpstrj" --file_xyz="$system_name.xyz" --rc="0"

# compute rdf
time $path_gij $r_max $dr 1 FULL "$system_name.dat" "$system_name.lammpstrj" gij_pairs.dat full $rdf_name > $rdf_name.log

# compute xrd
time python $path_python/pyxrd.py "$system_name.dat" $rdf_name.dat form_factors.dat $r_max $dr $q_max $dq -atomtype='full' > $xrd_name.log

cp o.Fx_Icoh.dat "$xrd_name".Icoh.dat
