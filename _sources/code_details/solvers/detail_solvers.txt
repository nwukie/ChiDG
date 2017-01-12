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
`Forward Euler`_                                    ``Forward Euler``
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

    \bigg(\frac{\boldsymbol{M}}{\Delta t} + \frac{\partial R(\hat{Q}^{n})}{\partial Q}\bigg) \Delta Q = -R(\hat{Q}^{n})



.. Diagonally-Implicit Runge-Kutta
.. -------------------------------


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








