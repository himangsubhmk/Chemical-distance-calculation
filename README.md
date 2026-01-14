To detect the strand breaking events one needs to compute the change of chemical distance of a pair of neighbouring particles. 

Calculation of chemical distance:


The code l-chem-dist.c will calculate the chemical distance. 

compilation: gcc l-chem-dist.c -lm

This code needs three arguments:  ./a.out N filepath-dump-series nframe

Where N is the number of particle(here 10000), filepath for example Data_full_e0-10_S0-1_B-0.01_steps-100000.dump

and nframe is the frame number for which it will calculate chemical-distance with respect to (nframe-1)th frame.

See the run-script_l-dist.sh that will calculate of the above dump-file series for nframe 1-99, note it starts with 0th-frame
