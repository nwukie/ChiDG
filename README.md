
# ChiDG

[![][license img]][license]

A Chimera-based, discontinuous Galerkin framework


## Features
 - Discontinuous Galerkin discretization
 - Chimera/Overset grid capability(no hole-cutting as of yet)
 - Intrinsic automatic differentiation
 - Composition-based modeling infrastructure
 - Time-integrators
    - Steady
    - Backward Euler
    - 3-stage diagonally-implicit Runge-Kutta
    - Harmonic Balance
 - Nonlinear solvers
    - Newton
    - Quasi-Newton
 - Linear solvers
    - FGMRES
 - Preconditioners
    - Identity
    - Block Jacobi
    - Block ILU0
    - Restricted Additive Schwarz + ILU0


## Example applications:
 - External aerodynamics 
 - Computational aeroacoustics
 - Turbomachinery
 - Moving grid problems

<img src="doc/figures/naca2412_M0p2_A4p0_cp_P3_thickzoom.png" height="200"> <img src="doc/figures/mode91_3dview_cropped.png" height="200"> <img src="doc/figures/turbine_HB_P3.png" height="200"> <img src="doc/figures/naca_ale.png" height="200">






## Version
Pre-release. Working towards releasing v0.1 soon.

## Documentation

Documentation can be found on the following github page:

[ChiDG Documentation](https://nwukie.github.io/ChiDG/ )


## Installation

[Instructions for building ChiDG](http://nwukie.github.io/ChiDG/getting_started/getting_started.html#build-from-source )


## License
ChiDG is released under the BSD 3-clause license. See LICENSE file.


## Author Acknowledgement:
Nathan A. Wukie   <nathan.wukie@gmail.com>


## Developers:
 - Mayank Sharma
 - Matteo Ugolotti
 - Dr. Eric Wolf


## Sponsorship Acknowledgement:
This material is based upon work supported by the National Science Foundation Graduate 
Research Fellowship Program under Grant No. 1610397. Any opinions, findings, and 
conclusions or recommendations expressed in this material are those of the author(s) 
and do not necessarily reflect the views of the National Science Foundation.


Some contributions to ChiDG were made by the Air Force Research Laboratory at 
Wright-Patterson Airforce Base in Dayton, Ohio. These contributions are marked 
by @author notes in the source appended with (AFRL). These contributions were 
cleared for resease to the public domain under 
Case Numbers:88ABW-2016-5744(5/1/2016-9/30/2016), 88ABW-2017-5059(5/1/2017-8/30/2017).





[license]:LICENSE
[license img]:https://img.shields.io/badge/license-BSD%203--clause-blue.svg























































