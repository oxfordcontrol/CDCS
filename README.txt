================================================================================
                                CDCS README.TXT
================================================================================

CDCS (Cone Decomposition Conic Solver) is an open-source MATLAB solver for sparse
conic programs with partially decomposable conic constraints. CDCS implements the
alternating direction method of multipliers (ADMM) described in on our papers
 
* "Fast ADMM for Semidefinite Programs with Chordal Sparsity" available from
(https://arxiv.org/pdf/1609.06068v2.pdf)

* "Fast ADMM for homogeneous self-dual embeddings of sparse SDPs" available from
(https://arxiv.org/abs/1611.01828)

Current version: 1.1
Release notes: CDCS is based on a temporary research code called ADMM-PDCP, which
               is no longer maintained. If you downloaded ADMM-PDCP, please 
               replace it with CDCS.


================================================================================
                                   CONTENTS
================================================================================

* Description
* Quick start
* How to cite
* Contact us
* Licence


================================================================================
                                   DESCRIPTION
================================================================================

CDCS solves in the standard primal and dual vectorized forms

		minimize 	c'x						maximize 	b'y
	(1)	subject to	Ax = b,				(2)	subject to	A'y + z = c,	
					x \in K								z \in K*

where the conic constraint `x \in K` are partially decomposable. This means that
`x \in K` can be replaced by `p` smaller conic constraints `x_1 \in K_1`, ..., 
`x_p \in K_p`, where `x_1`, ..., `x_p` are (possibly not-disjoint) subsets of the
original optimization variable `x`.

CDCS supports cartesian products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

Currently, CDCS only decomposes semidefinite cones characterized by a chordal 
sparsity pattern. The other supported cone types are not decomposed. 
This means that CDCS is most suitable for large sparse semidefinite programs 
(SDPs), although it can be used for any conic program over the supported cones.


================================================================================
                                 QUICK START
================================================================================

To install CDCS, simply run the installer script in MATLAB:

	>> cdcsInstall;

To test your installation, run 

	>> cdcsTest;
	
CDCS is called with the syntax

	>> [x,y,z,info] = cdcs(At,b,c,K,options);
	
where `At` is the transpose of the matrix `A` in problems (1)-(2) above. 
Note that the inputs and outputs are in the same format used by SeDuMi. Type

	>> help cdcsOpts
	
for a complete list of solver options.
	
NOTE:this is a research code, and is under active development. You may find 
some undocumented inputs and options that are being used for development 
purposes, in the hope that they will become part of the "official" release. If 
you have any suggestions for improvement, or find any bugs, feel free to contact
us (see the Contact Us section below).


================================================================================
                                HOW TO CITE
================================================================================

If you find CDCS useful, please cite

@article{{ZFPGWhsde2016,
    archivePrefix = {arXiv},
	eprint = {1611.01828},
	primaryClass = "math-OC",
	author = {Zheng, Yang and Fantuzzi, Giovanni and Papachristodoulou, Antonis and Goulart, Paul and Wynn, Andrew},
	title = {{Fast ADMM for homogeneous self-dual embeddings of sparse SDPs}}
	}

@article{ZFPGWpd2016,
	archivePrefix = {arXiv},
	eprint = {1609.06068v2},
	primaryClass = "math-OC",
	author = {Zheng, Yang and Fantuzzi, Giovanni and Papachristodoulou, Antonis and Goulart, Paul and Wynn, Andrew},
	title = {{Fast ADMM for Semidefinite Programs with Chordal Sparsity}}
	}
	
@misc{CDCS,
    author       = {Zheng, Yang and Fantuzzi, Giovanni and Papachristodoulou, Antonis and Goulart, Paul and Wynn, Andrew},
    title        = {{CDCS}: Cone Decomposition Conic Solver, version 1.0},
    howpublished = {\url{https://github.com/giofantuzzi/CDCS}},
    month        = Sep,
    year         = 2016
    }


A selection of BibTex styles that support arXiv preprints can be found at
http://arxiv.org/hypertex/bibstyles/


================================================================================
                                 CONTACT US
================================================================================

To contact us about CDCS, suggest improvements and report bugs, email:
Giovanni Fantuzzi: gf910@ic.ac.uk
Yang Zheng	     : yang.zheng@eng.ox.ac.uk


================================================================================
                                   LICENCE
================================================================================

CDCS is free software; you can redistribute it and/or modify it under the terms 
of the GNU Lesser General Public Licence (LGPL) as published by the Free Software
Foundation; either version 3 of the Licence, or (at your option) any later 
version.

CDCS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with CDCS; if not, write to the Free Software Foundation, Inc., 51 Franklin St, 
Fifth Floor, Boston, MA 02110-1301 USA.


================================================================================
                               END OF README
================================================================================