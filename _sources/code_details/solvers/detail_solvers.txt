=======
Solvers
=======

The algorithm choices for the simulation environment are set using the 
``chidg%set`` method with a selection parameter and a string option:

::

    call chidg%set('Time Integrator', String)
    call chidg%set('Nonlinear Solver',String)
    call chidg%set('Linear Solver',   String)
    call chidg%set('Preconditioner',  String)


----------------
Time Integrators 
----------------

================================================    ==================
Algorithm                                           String Selector
================================================    ==================
`Steady`_                                           ``Steady``
`Backward Euler`_                                   ``Backward Euler``
`Diagonally-Implicit Runge-Kutta`_                  ``DIRK``
`Harmonic Balance`_                                 ``Harmonic Balance``
`Forward Euler`_                                    ``Forward Euler``
`Explicit Runge Kutta`_                             ``Runge-Kutta Method``
================================================    ==================



ChiDG Representation
====================


Integration of the temporal derivative term


.. math::

    \int \psi \frac{\partial Q}{\partial t} d\Omega

in

.. math::

    \int_{\Omega_e} \boldsymbol{\psi} \frac{\partial Q}{\partial t} d\Omega +
    \int_{\Gamma_e} \boldsymbol{\psi} \vec{F} \cdot \vec{n} d\Gamma - 
    \int_{\Omega_e} \nabla \boldsymbol{\psi} \cdot \vec{F} d\Omega + 
    \int_{\Omega_e} \boldsymbol{\psi} S d\Omega = 0


is handled by the ``time_integrator_t`` class.




.. code-block:: Fortran

    type, abstract, public :: time_integrator_t

    contains

        procedure(data_interface), deferred :: iterate   ! Must define this procedure in the extended type

    end type time_integrator_t


Note that in the terminology used in this section, :math:`M` is the Mass Matrix

.. math::

    \boldsymbol{M} = \int_{\Omega} \boldsymbol{\psi} \boldsymbol{\psi} d\Omega

and :math:`R(Q)` is the right-hand-side vector of spatial operators

.. math::

    R(\hat{Q}) = \int_{\Gamma_e} \boldsymbol{\psi} \vec{F} \cdot \vec{n} d\Gamma - 
                 \int_{\Omega_e} \nabla \boldsymbol{\psi} \cdot \vec{F} d\Omega + 
                 \int_{\Omega_e} \boldsymbol{\psi} S d\Omega

Implicit integrators
====================


Steady
------

.. math::

    \int \psi \frac{\partial Q}{\partial t} d\Omega = 0
    
Nonlinear problem:

.. math::

    R(Q) = 0

Newton Linearization:

.. math::

    \frac{\partial R(\hat{Q}^{n})}{\partial Q} \Delta Q = -R(\hat{Q}^{n})


Backward Euler
--------------

Solution advancement via first-order backward difference of the 
time derivative:

.. math::

    \boldsymbol{M} \frac{\hat{Q}^{n+1} - \hat{Q}^{n}}{\Delta t} = R(\hat{Q}^{n+1})


Nonlinear problem:

.. math:: 

    \frac{\Delta \hat{Q}}{\Delta t}\boldsymbol{M} + R(\hat{Q}^{n+1}) = 0

Newton Linearization:

.. math::

    \bigg(\frac{\boldsymbol{M}}{\Delta t} + \frac{\partial R(\hat{Q}^{m})}{\partial Q}\bigg) \delta Q^{m} & = -\frac{\boldsymbol{M}}{\Delta t}\Delta Q^{m} -R(\hat{Q}^{m})\\
    \hat{Q}^{m} & = \hat{Q}^{n} + \Delta Q^{m}

where, :math:`\delta Q^{m} = \Delta Q^{m + 1} -\Delta Q^{m}` for the :math:`m^{th}` Newton iteration.

.. note::
    
    - The ChiDG nonlinear solver requires the assembly of the linearized system outlined above
    - The nonlinear solver computes :math:`\hat{Q}^{m}` instead of :math:`\Delta \hat{Q}^{m}`
    - For the assembly of the rhs of the linearized system, :math:`\Delta \hat{Q}^{m} = \hat{Q}^{m} - \hat{Q}^{n}` 

Diagonally-Implicit Runge-Kutta
-------------------------------

.. math::

    \boldsymbol{M}\frac{\partial Q}{\partial t} + R(Q) = 0

With the coefficient arrays associated with the diagonally implicit Runge-Kutta method:

.. math::

    \boldsymbol{A} = \left[\begin{array}{ccc}
                           \gamma & 0 & 0 \\
                           \tau_{2} - \gamma & \gamma & 0 \\
                           b_{1} & b_{2} & \gamma \end{array} \right]

.. math::

    \boldsymbol{b} = \left[\begin{array}{ccc}
                             b_{1} & b_{2} & \gamma \end{array} \right]

where :math:`\gamma` is the root of :math:`x^{3} - 3x^{2} + \frac{3}{2}x - \frac{1}{6} = 0 \in \left(\frac{1}{6},\frac{1}{2}\right)` and

.. math::

    \tau_{2} & = (1 + \gamma)/2\\
    b_{1} & = -(6\gamma^{2} - 16\gamma + 1)/4\\
    b_{2} & = (6\gamma^{2} - 20\gamma + 5)/4

The solution is advanced in time as:

.. math::

    \hat{Q}^{n + 1} = \hat{Q}^{n} + b_{1}\Delta \hat{Q}_{1} + b_{2}\Delta \hat{Q}_{2} + b_{3} \Delta \hat{Q}_{3}

Implicit system:

.. math::

    \frac{\Delta \hat{Q}_{i}}{\Delta t}\boldsymbol{M} = -R\left(\hat{Q}^{n} + \sum_{j = 1}^{i}A_{ij}\Delta \hat{Q}_{i}\right)\;\;\;\text{for}\;i = 1,3

Newton linearization:

.. math::

    \left(\boldsymbol{M} + \gamma \Delta t \frac{\partial R(\hat{Q}^{m}_{i})}{\partial Q}\right)\delta \hat{Q}^{m}_{i} = -\boldsymbol{M}\Delta \hat{Q}^{m}_{i} - \Delta t R\left(\hat{Q}^{m}_{i}\right)\;\;\;\text{for}\;i = 1,3

with

.. math::

    \hat{Q}^{m}_{i} & = \hat{Q}^{n} + \sum_{j = 1}^{i - 1}A_{ij}\Delta \hat{Q}_{i} + \gamma \Delta \hat{Q}^{m}_{i}\\
    \delta \hat{Q}^{m}_{i} & = \Delta \hat{Q}^{m + 1}_{i} - \Delta \hat{Q}^{m}_{i}

.. note::

    - The ChiDG nonlinear solver requires the assembly of the stagewise linearized systems
    - The nonlinear solver computes :math:`\hat{Q}^{m}_{i}` instead of :math:`\Delta \hat{Q}^{m}_{i}`
    - For the assembly of the rhs of the linearized system, :math:`\Delta \hat{Q}^{m}_{i} = (\hat{Q}^{m}_{i} - \hat{Q}^{n} - \sum_{j = 1}^{i - 1}A_{ij}\Delta \hat{Q}_{i})/\gamma`

Harmonic Balance
----------------

Consider a set of :math:`N` independent equations:

.. math::

    \frac{\partial \hat{\boldsymbol{Q}}^{*}}{\partial t} + \nabla \cdot \boldsymbol{F}^{*} + \boldsymbol{S}^{*} = 0

where

.. math::

    \hat{\boldsymbol{Q}^{*}} & = \left[\hat{Q}_{1}, \hat{Q}_{2}, \cdots, \hat{Q}^{N}\right]^{T}\\
    \boldsymbol{F}^{*} & = \left[\vec{F}(\hat{Q}_{1}, \nabla \hat{Q}_{1}), \vec{F}(\hat{Q}_{2}, \nabla \hat{Q}_{2}), \cdots, \vec{F}(\hat{Q}_{N}, \nabla \hat{Q}_{N})\right]^{T}\\
    \boldsymbol{S}^{*} & = \left[S(\hat{Q}_{1}, \nabla \hat{Q}_{1}), S(\hat{Q}_{2}, \nabla \hat{Q}_{2}), \cdots, S(\hat{Q}_{N}, \nabla \hat{Q}_{N})\right]^{T}

In the harmonic balance method, a conservative solution vector at any instant of time is represented as a Fourier series in time as:

.. math::

    \hat{Q}_{n} = A_{0} + \sum_{k = 1}^{K}\left[A_{k}\text{sin}(\omega_{k}t_{n}) + B_{k}\text{cos}(\omega_{k}t_{n})\right]

with :math:`K` frequencies, :math:`\boldsymbol{\omega} = [\omega_{1}, \omega_{2}, \cdots, \omega_{K}]` and the instant of time :math:`t_{n}` belongs to the set of time levels,
:math:`\boldsymbol{t} = [t_{1}, t_{2}, \cdots, t_{N}]` with :math:`N = 2K + 1`. Thus, the series of conservative solution vectors can be related to the Fourier coefficients vectors,
:math:`\hat{\boldsymbol{Q}}_{F}` as:

.. math::

    \hat{\boldsymbol{Q}}^{*} = E^{-1}\hat{\boldsymbol{Q}}_{F}

Defining the pseudo spectral operator as,

.. math::

    D = \frac{\partial E^{-1}}{\partial t}E

which couples :math:`\hat{\boldsymbol{Q}}^{*}` such that the conservative solutions satisfy time-varying sinusoidal functions according to their Fourier representation, the governing 
equation can be rewritten as the Harmonic Balance equation:

.. math::

    D\hat{\boldsymbol{Q}}^{*} + \nabla \boldsymbol{F}^{*} + S^{*} = 0

Multiplying with a column of test functions, :math:`\psi` and applying Gauss' divergence theorem provides the working form of the Harmonic Balance equation:

.. math::

    \int_{\Omega_{e}}\psi D\hat{\boldsymbol{Q}}^{*}d\Omega + \int_{\Gamma_{e}}\boldsymbol{F}^{*} \cdot \vec{n}d\Gamma - 
    \int_{\Omega_{e}}\nabla \psi \cdot \boldsymbol{F}^{*}d\Omega + \int_{\Omega_{e}}\psi \boldsymbol{S}^{*}d\Omega = 0

Newton Linearization:

Consider,

    .. math::

        \mathscr{R}^{*} & = \int_{\Gamma_{e}}\boldsymbol{F}^{*} \cdot \vec{n}d\Gamma - \int_{\Omega_{e}}\nabla \psi \cdot \boldsymbol{F}^{*}d\Omega + 
        \int_{\Omega_{e}}\psi \boldsymbol{S}^{*}d\Omega\\
        \mathscr{D}^{*} & = \int_{\Omega_{e}}\psi D\hat{\boldsymbol{Q}}^{*}d\Omega

Then, Newton linearization of the Harmonic Balance system of equations is:

.. math::

    \left(\frac{\partial \mathscr{D}^{*}}{\partial \hat{\boldsymbol{Q}}^{*}} + \frac{\partial \mathscr{R}^{*}}{\partial \hat{\boldsymbol{Q}}^{*}}\right)\Delta \hat{\boldsymbol{Q}^{*}} = 
    -(\mathscr{D}^{*} + \mathscr{R}^{*})

Explicit integrators
====================

Forward Euler
-------------

Solution advancement via a first-order forward difference of the
time derivative:

.. math::

    \boldsymbol{M} \frac{\hat{Q}^{n+1} - \hat{Q}^{n}}{\Delta t} = R(\hat{Q}^{n})




Algebraic problem:

.. math:: 

    \frac{\Delta \hat{Q}}{\Delta t}\boldsymbol{M} + R(\hat{Q}^{n}) = 0


Solution iterated in time as:

.. math::

    \hat{Q}^{n+1} = \hat{Q}^n - {\Delta t} \boldsymbol{M}^{-1}R(\hat{Q}^{n})

Explicit Runge Kutta
--------------------

For a general explicit runge Kutta method with :math:`s` stages:

.. math::

    \boldsymbol{M}\frac{\hat{Q}^{n + 1} - \hat{Q}^{n}}{\Delta t} = \sum_{i = 1}^{s}b_{i}\Delta \hat{Q}_{i}

where

.. math::
    
    \Delta \hat{Q}_{i} = -R\left(\hat{Q}^{n} + \sum_{j = 1}^{i - 1}a_{ij}\Delta \hat{Q}_{j}\right)

Algebraic problem:

.. math::

    \frac{\Delta \hat{Q}}{\Delta t}\boldsymbol{M} - \sum_{i = 1}^{s}b_{i}\Delta \hat{Q}_{i} = 0

Solution iterated in time as:

.. math::

    \hat{Q}^{n + 1} = \hat{Q}^{n} + \Delta t \boldsymbol{M}^{-1}\sum_{i = 1}^{s}b_{i}\Delta \hat{Q}_{i}


|
|
|
|
|
|
|
|


-----------------
Nonlinear solvers
-----------------



================================================    ==================
Algorithm                                           String Selector
================================================    ==================
`Newton`_                                           ``Newton``
`Quasi-Newton`_                                     ``Quasi-Newton``
================================================    ==================




ChiDG Representation
====================
ChiDG includes nonlinear solvers for solving the nonlinear sets of partial 
differential equations. In general, the implicit problem statement here is:

    - Find :math:`Q`, such that :math:`\mathcal{R}(Q) = 0`



.. code-block:: Fortran

    type, abstract, public :: nonlinear_solver_t

    contains

        procedure(data_interface), deferred :: solve   ! Must define this procedure in the extended type

    end type nonlinear_solver_t

    



Newton
======



The Full-Newton solver solves the equation

.. math::

    \mathcal{R}(Q) = 0


where :math:`\mathcal{R}(Q)` is some potentially nonlinear function of the solution. This depends on 
the discretization, the equation set, the solution order, and the time-integration scheme. 
The Newton solver linearizes the problem and computes an update 
of the solution by solving 

.. math::

    \frac{\partial \mathcal{R}}{\partial Q} \Delta Q = -\mathcal{R}

So, at each Newton step, a linear system of equations is being solved for :math:`\Delta Q`.
Once the update is solved for, the solution vector is updated as

.. math::

    Q^{n+1} = Q^{n} + \Delta Q

Considerations:
---------------
One item to consider when using the Full-Newton solver is that the Newton
linearization(direction) is dependent on the current solution. Without a 
reasonable initial guess, Newton's method can diverge by sending the 
solution too far in the wrong direction.


|
|
|
|



Quasi-Newton
============

The Quasi-Newton solver solves a modified set of equations

.. math::

    \int_{\Omega_e} \psi \frac{\partial Q}{\partial \tau} d\Omega + \mathcal{R}(Q) = 0


Note the addition of a pseudo-time term to the nonlinear system of equations. This is
an effort increase robustness of the nonlinear solver by limiting the size of the solution
update in a single Newton step. This is accomplished by adding the time-scaling to the 
diagonal of the Jacobian matrix, increasing the diagonal dominance of the matrix, and 
limiting the size of the soltion update. As the solution progresses, the pseudo-timestep is
increased the pseudo-time derivative goes to zero and the original system of equations
is recovered.

The Quasi-Newton solver linearizes the problem including the pseudo-time scaling of the
solution update to the system of equations as

.. math::

    \int_{\Omega_e} \psi \frac{\Delta Q}{\Delta \tau} d\Omega + \frac{\partial \mathcal{R}}{\partial Q} \Delta Q = -\mathcal{R}


So, at each Quasi-Newton step, a linear system of equations is being solved for 
:math:`\Delta Q`. Once the update is solved for, the solution vector is updated as

.. math::

    Q^{n+1} = Q^{n} + \Delta Q


At each Quasi-Newton step, the pseudo-time step is updated as

.. math::

    d\tau = \frac{CFL^n h_e}{\bar{\lambda_e}}

where :math:`h_e = \sqrt[3]{\Omega_e}` and :math:`\bar{\lambda_e} = |\bar{V_e}| + \bar{c}`
is a mean characteristic speed. The CFL term is computed from the ratio of the initial
and current residual norms as

.. math::

    CFL^n = CFL^0 \frac{||\mathcal{R}(Q^0)||_2}{||\mathcal{R}(Q^n)||_2}


Options:
---------

    - CFL0: The initial CFL factor


|
|
|
|
|
|
|

--------------
Linear Solvers
--------------




================================================    ===============
Algorithm                                           String Selector
================================================    ===============
`Flexible Generalized Minimum Residual`_            ``FGMRES``
================================================    ===============

ChiDG Representation
====================



Flexible Generalized Minimum Residual
=====================================

A flexible version of the Generalized Minimum Residual(FGMRES) algorithm,
which is an iterative method for solving linear systems of equations. The FGMRES
algorithm allows the GMRES algorithm to be preconditioned in a flexible way such that 
the solution can be easily reconstructed.

Options:
--------

    - m:  Number of iterations before the algorithm is restarted 






|
|
|
|
|
|
|





---------------
Preconditioners
---------------



================================================    ===============
Algorithm                                           String Selector
================================================    ===============
`block-Jacobi`_                                     ``Jacobi``
`block-ILU0`_                                       ``ILU0``
`Restricted Additive Schwarz + block-ILU0`_         ``RAS+ILU0``
================================================    ===============


ChiDG Representation
====================

block-Jacobi
============






block-ILU0
==========





Restricted Additive Schwarz + block-ILU0
========================================








