================================================================================
                        ./CDCS/packages/+cdcs_sos README
================================================================================

This folder is a MATLAB package and contains the functions used by the 
homogeneous self-dual embedding solvers for SOS problems in CDCS. The package must contain the 
functions

makeADMM.m: constructs the operators for the ADMM steps and to check the
            convergence of the ADMM algorithm.
            
makeVariables.m: initialize the variables in the ADMM algorithm.

printHeader.m : a function to set headers & lines for fprintf.

All other package-specific functions can be put in a private subfolder for better
code organization, although it is not necessary.