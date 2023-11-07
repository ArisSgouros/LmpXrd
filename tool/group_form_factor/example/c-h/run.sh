#!/bin/bash
python ../../group_ff.py group_c-h.xyz form_factors.dat 20.0 0.1 --orig="0.0,0.0,0.0" --scheme="narten" --ff_group_file="o.ff_group_narten.dat" > o.log_narten

python ../../group_ff.py group_c-h.xyz form_factors.dat 20.0 0.1 --orig="0.0,0.0,0.0" --scheme="sum" --ff_group_file="o.ff_group_sum.dat" > o.log_sum
