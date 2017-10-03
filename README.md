# CDCS

CDCS (Cone Decomposition Conic Solver) is an open-source MATLAB solver for sparse conic programs with partially decomposable conic constraints. CDCS implements the alternating direction method of multipliers (ADMM)
described in our paper [_Chordal decomposition in operator-splitting methods for
sparse semidefinite programs_](https://arxiv.org/pdf/1707.05058.pdf). 

* Previous conference papers are [_Fast ADMM for Semidefinite Programs with Chordal Sparsity_](https://arxiv.org/pdf/1609.06068v2.pdf) and [_Fast ADMM for homogeneous self-dual embeddings of sparse SDPs_](https://arxiv.org/pdf/1611.01828.pdf) (included in the `doc/` folder)

**Current version:** 1.1.0

**Release notes:** 

* Homogeneous self-dual embedding is the new default method 
* The termination codes have changed. This means that if you use CDCS from YALMIP, the termination code returned by YALMIP will be incorrect. This should be fixed in the next YALMIP release!
* CDCS is based on a temporary research code called ADMM-PDCP, which is no longer maintained. If you downloaded ADMM-PDCP, please replace it with CDCS.


## Contents
* [Description](#Description)
* [Quick start](#QuickStart)
* [How to cite](#References)
* [Contact us](#Contacts)
* [Licence](#Licence)


## Description<a name="Description"></a>

CDCS solves in the standard primal and dual vectorized forms

		minimize 	c'x					maximize 	b'y
	(1)	subject to	Ax = b,				(2)	subject to	A'y + z = c,	
				x \in K							z \in K*

where the conic constraint `x \in K` are partially decomposable. This means that `x \in K` can be replaced by `p` smaller conic constraints `x_1 \in K_1`, ...,  `x_p \in K_p`, where `x_1`, ..., `x_p` are (possibly not-disjoint) subsets of the original optimization variable `x`.

CDCS supports cartesian products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

Currently, CDCS only decomposes semidefinite cones characterized by a chordal sparsity pattern. The other supported cone types are not decomposed.  This means that CDCS is most suitable for large sparse semidefinite programs (SDPs), although it can be used for any conic program over the supported cones.

CDCS offers a choice to solve the primal problem (1) only, the dual problem (2) only, or the homogeneous self-dual embedding of the two problems. _From version 1.1.0, the homogeneous self-dual embedding is the default method._


## Quick start<a name="QuickStart"></a>

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
	
**NOTE:** _this is a research code, and is under active development. You may find 
some undocumented inputs and options that are being used for development 
purposes, in the hope that they will become part of the "official" release. If 
you have any suggestions for improvement, or find any bugs, feel free to [contact us](#Contacts)!_


## How to cite<a name="References"></a>

If you find CDCS useful, please cite at least one of the following papers as appropriate:

```
@article{{ZFPGWchordal2017,
    archivePrefix= {arXiv},
    eprint       = {1707.05058},
    primaryClass = "math-OC",
    author       = {Zheng, Yang and Fantuzzi, Giovanni and Papachristodoulou, Antonis and Goulart, Paul and Wynn, Andrew},
    title        = {{Chordal decomposition in operator-splitting methods for sparse semidefinite programs}}
    }
	
@misc{CDCS,
    author       = {Zheng, Yang and Fantuzzi, Giovanni and Papachristodoulou, Antonis and Goulart, Paul and Wynn, Andrew},
    title        = {{CDCS}: Cone Decomposition Conic Solver, version 1.1},
    howpublished = {\url{https://github.com/giofantuzzi/CDCS}},
    month        = Sep,
    year         = 2016
    }
```
A selection of BibTex styles that support arXiv preprints can be found [here](http://arxiv.org/hypertex/bibstyles/).


## Contact us<a name="Contacts"></a>
To contact us about CDCS, suggest improvements and report bugs, email either [Giovanni Fantuzzi](mailto:gf910@ic.ac.uk?Subject=CDCS) or [Yang Zheng](mailto:yang.zheng@eng.ox.ac.uk?Subject=CDCS).


## Licence<a name="Licence"></a>

CDCS is free software; you can redistribute it and/or modify it under the terms 
of the [GNU Lesser General Public Licence (LGPL)](https://www.gnu.org/licenses/lgpl-3.0.en.html) as published by the Free Software
Foundation; either version 3 of the Licence, or (at your option) any later version.

CDCS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. [See the GNU General Public License for more details](https://www.gnu.org/licenses/gpl-3.0.en.html).

You should have received a copy of the GNU Lesser General Public License along 
with CDCS; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
