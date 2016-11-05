.. ChiDG documentation master file, created by
   sphinx-quickstart on Tue Aug  9 21:10:34 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================
Welcome to ChiDG
================


.. .. image:: figures/naca2412.png
    :width: 18 %
    :align: left

.. .. image:: figures/mode91_3dview.png
    :width: 16 %
    :align: center

.. .. image:: figures/aachen_turbine_nonreflectingbc.png
    :width: 10 %
    :align: right


.. +----------+----------+-----------+
.. |  |test1| + |test2|  + |test3|   +
.. +----------+----------+-----------+


|

ChiDG is:
---------

a framework for solving sets of partial differential 
equations using the discontinuous Galerkin method on Chimera-overset grids.





.. image:: figures/icon_dg.png
    :width: 19 %
    :target: ./code_details/equations/detail_equations.html
.. image:: figures/icon_curved.png
    :width: 19%
    :target: ./code_details/mesh/detail_mesh.html
.. image:: figures/icon_chimera.png
    :width: 19%
.. image:: figures/icon_newton.png
    :width: 19%
.. image:: figures/icon_autodiff.png
    :width: 19%
    :target: ./code_details/autodiff/detail_autodiff.html


|
|

ChiDG exists as:
----------------

a Fortran Class that provides an interface to the simulation environment.

::

    type :: chidg_t
        type(chidg_data_t)          data
        type(time_integrator_t)     time_integrator
        type(nonlinear_solver_t)    nonlinear_solver
        type(linear_solver_t)       linear_solver
        type(preconditioner_t)      preconditioner
        ...
    end type chidg_t





|
|



Example applications:
--------------------------------

.. image:: figures/naca2412.png
    :width: 36 %
    :align: left
.. image:: figures/mode91_3dview.png
    :width: 32 %
    :align: right
.. image:: figures/aachen_turbine_nonreflectingbc.png
    :width: 20 %
    :align: center




|
|



.. toctree::
    :maxdepth: 4
    :hidden:
 
    getting_started/getting_started
    examples/example_main
    code_details/details_main




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

