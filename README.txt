================================================================================
                            ADMM-PDCP README.TXT
================================================================================

An open-source MATLAB ADMM solver for partially decomposable conic optimization 
programs.

Contents
--------
* Introduction
* Quick start
* Planned features
* How to cite
* Contact us



================================================================================
                                 INTRODUCTION
================================================================================

ADMM-PDCP is a MATLAB&reg; implementation of the alternating direction method of multipliers (ADMM) 
for conic programs in the standard primal and dual vectorized forms

    minimize    c'x                     maximize    b'y
(1) subject to  Ax = b,             (2) subject to  c - A'y = z,    
                x \in K                             z \in K*

where the conic constraint x \in K can be replaced by p smaller conic constraints
x_1 \in K_1, ..., x_p \in K_p, where x_1, ..., x_p are (possibly not-disjoint) 
subsets of the original optimization variable x. ADMM-PDCP supports cartesian 
products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

At present, only chordal decomposition techniques for semidefinite cones are 
implemented, while the other supported cone types are not decomposed. 
Consequently, ADMM-PDCP is most suitable for large semidefinite programs (SDPs) 
characterized by chordal sparsity, although it can be used for any conic program
over the supported cones.



================================================================================
                                 QUICK START
================================================================================

To install ADMM-PDCP, simply run the installer script at the MATLAB command line:

	>> admmPDCPInstall;

To test your installation, run

	>> admmPDCPTest;

ADMM-PDCP is called with the syntax

	>> [x,y,z,info] = admmPDCP(At,b,c,K);
	
where At is the transpose of the data matrix in problem (1) above (the inputs 
and outputs are in the same format used by SeDuMi). Type

	>> help admmPDCP
	
for information on how to specify additional solver options.
	
NOTE: this is a research code, and is under active development. You may find 
some undocumented inputs and options that are being used for development 
purposes, in the hope that they will become part of the "official" release. If 
you have any suggestions for improvement, or find any bugs, feel free to contact
us (see the Contact Us section below).


================================================================================
                              PLANNED FEATURES
================================================================================

* YALMIP interface


================================================================================
                                 HOW TO CITE
================================================================================

If you find ADMM-PDCP useful, please cite:

1. Y. Zheng, G. Fantuzzi, A. Papachristodoulou, P. J. Goulart, A. Wynn, "Fast ADMM
   for Semidefinite Programs with Chordal Sparsity", arXiv:1609.06068 [math.OC].
   Available from https://arxiv.org/abs/1609.06068


================================================================================
                                 CONTACT US
================================================================================

To contact us about ADMM-PDCP, suggest improvements and report bugs, email

* gf910@ic.ac.uk (Giovanni Fantuzzi)
* yang.zheng@eng.ox.ac.uk	(Yang Zheng)



================================================================================
                              END OF README.TXT
================================================================================
