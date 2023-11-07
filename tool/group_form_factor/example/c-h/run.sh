#!/bin/bash

# group form factor w/ narten scheme
python ../../group_ff.py group_c-h.xyz form_factors.dat 20.0 0.1 --ff_group_file="o.ff_group_narten.dat" \
    --scheme="narten" --init_guess="1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0" --orig="0.0,0.0,0.0" > o.log_narten

# group form factor w/ sum scheme
python ../../group_ff.py group_c-h.xyz form_factors.dat 20.0 0.1 --ff_group_file="o.ff_group_sum.dat" \
    --scheme="sum" --init_guess="1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0" > o.log_sum
