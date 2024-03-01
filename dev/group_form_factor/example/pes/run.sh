#!/bin/bash

# sum
python ../../group_ff.py group_pes_ch3.xyz form_factors.dat 20.0 0.1 --ff_group_file="o.ff_pes_ch3_sum.dat" --scheme="sum" -weights=1 \
  -init_guess="2.3100,20.8439,1.0200,10.2075,1.5886,0.5687,0.8650,51.6512,0.2156"

# narten, orig=0,2,0
python ../../group_ff.py group_pes_ch3.xyz form_factors.dat 20.0 0.1 --ff_group_file="o.ff_pes_ch3_narten_0_2_0.dat" --scheme="narten" --orig="0.0,2.0,0.0" -weights=1 \
  -init_guess="2.3100,20.8439,1.0200,10.2075,1.5886,0.5687,0.8650,51.6512,0.2156"

# narten, orig=0,1,0
python ../../group_ff.py group_pes_ch3.xyz form_factors.dat 20.0 0.1 --ff_group_file="o.ff_pes_ch3_narten_0_1_0.dat" --scheme="narten" --orig="0.0,1.0,0.0" -weights=1 \
  -init_guess="2.3100,20.8439,1.0200,10.2075,1.5886,0.5687,0.8650,51.6512,0.2156"
