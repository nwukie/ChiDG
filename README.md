
<p align="center">
    <img src=/doc/figures/chidg_logo_small.png?raw=true />
</p>

# ChiDG

[![][license img]][license]
[![Build Status](https://travis-ci.org/nwukie/ChiDG.svg?branch=master)](https://travis-ci.org/nwukie/ChiDG)

A Chimera-based, discontinuous Galerkin solver


Overset airfoil grid                              |  NACA 2412 pressure
:------------------------------------------------:|:------------------------------------------------------:
![](doc/figures/naca2412_A4p0_straight_grid.png)  |     ![](doc/figures/naca2412_M0p2_A4p0_cpcontour_P3.png)



4th-order duct grid                               |  Acoustic duct mode radiation. 7th-order accuracy.
:------------------------------------------------:|:------------------------------------------------------:
![](doc/figures/munt_duct_grid.png)               |     ![](doc/figures/mode91_3dview.png)  



Constant pressure outlet boundary condition       |  Fully-implicit nonreflecting outlet boundary condition
:------------------------------------------------:|:-------------------------------------------------------:
<img src="doc/figures/aachen_turbine_reflectingbc.png" hspace="118pt"/> |   <img src="doc/figures/aachen_turbine_nonreflectingbc.png" hspace="118pt"/>





## Documentation

Documentation can be found on the following github page:

[ChiDG Documentation](https://nwukie.github.io/ChiDG/ )











## Installation

### Dependencies

CMake: Build system  
HDF5: File IO
BLAS/LAPACK: Optimized linear algebra
METIS: Domain-decomposition
MPI: Parallelization

[Instructions for building ChiDG](http://nwukie.github.io/ChiDG/getting_started/getting_started.html#build-from-source )





## License
ChiDG is released under the BSD 3-clause license. See LICENSE file.



## Author Acknowledgement:
Nathan A. Wukie   <nwukie@gmail.com>






## Sponsorship Acknowledgement:
This material is based upon work supported by the National Science Foundation Graduate 
Research Fellowship Program under Grant No. 1610397. Any opinions, findings, and 
conclusions or recommendations expressed in this material are those of the author(s) 
and do not necessarily reflect the views of the National Science Foundation.











[license]:LICENSE
[license img]:https://img.shields.io/badge/license-BSD%203--clause-blue.svg























































