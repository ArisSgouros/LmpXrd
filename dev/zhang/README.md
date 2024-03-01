# gr
Codes for evaluating the Radial Distribution Function

# Author
 - Dr. Kai Zhang

# Notes
 - The code was employed for extra validation purposes.
 - The original code was retrieved from:
   https://github.com/statisticalmechanics/gr
 - The code has been improved slightly:
     *fix issue with comment section in .xyz files (3ac673b)
     *improve IO format (ce95785)

# Organization
The folder includes the following files and directories:
 - README   -> current file
 - gr.c     -> C code to compute the RDF
 - sq.c     -> C code to compute the scattering factor w/ the direct method
 - make.sh  -> Bash script for make
 - clean.sh -> Bash script for cleaning binary files
