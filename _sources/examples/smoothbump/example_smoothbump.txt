.. _example_smoothbump:

=================
Euler smooth bump
=================


----------
Background
----------
A canonical problem for verifying high-order solvers is the Euler smooth bump test case. 
In this example, the geometry is defined as a channel with a small smooth bump profile 
on the lower wall, as shown here:

.. image:: smoothbump_geometry.png
    :width: 500px
    :align: center


The Euler equations are solved on this geometry, which for subsonic flow are isentropic.
However, in the numerical solution, some entropy is generated due to the discretization.
This gives us a quantity that we can measure in the numerical solution by which we can
gauge the accuracy of the numerical scheme. The measure of the entropy error is integrated 
over the volume of the domain and computed as:

.. math:: 

    Entropy Error = \sqrt{  \frac{\int\Bigg( \frac{ \frac{p}{\rho^\gamma} - \frac{p_\infty}{\rho_\infty^\gamma}
                                                  }
                                                  { \frac{p_\infty}{\rho_\infty^\gamma}     } \Bigg)^2 d\Omega}
                                                  {     \int d\Omega }}




------
Setup
------


Setup from scratch:
-------------------
A Fortran code is attached that can be used to generate a Plot3D grid file for this case:

| :download:`smoothbump_generator.f90 <grids/smoothbump_generator.f90>`
| :download:`chidg.nml <grids/chidg.nml>`
|


| **Step 1:** compile the provided grid generator file using a Fortran compiler
| **Step 2:** run the executable to produce an unformatted double precision Plot3D grid file
| **Step 3:** convert the Plot3D grid file to a ChiDG-format grid file using ``chidg convert gridfile.x``
| **Step 4:** set boundary conditions on the ChiDG grid using the ``chidg edit gridfile.h5`` action
| **Step 5:** download a ``chidg.nml`` configuration file and place it in a run directory with the grid file
| **Step 6:** change to the run directory and execute ``chidg`` from the command line
| **Step 7:** once the simulation is converged, execute ``chidg post gridfile.h5`` to write a Tecplot file for visualization
|


Use setup provided:
-------------------
Alternatively, a ChiDG-formatted grid file is also provided with the boundary 
conditions already set:

| :download:`smoothbump.h5 <grids/smoothbump.h5>`
| :download:`chidg.nml <grids/chidg.nml>`
|

| **Step 1:** download the ``smoothbump.h5`` grid file and the ``chidg.nml`` settings file
| **Step 2:** place the grid and settings file in a new directory
| **Step 3:** change to the run directory and execute ``chidg`` from the command line
| **Step 4:** once the simulation is converged, execute ``chidg post smoothbump.h5`` to write a Tecplot file for visualization
|
|




-------
Results
-------

.. image:: smoothbump_cp_contours.png
    :width: 80 %
    :align: center

.. image:: smoothbump_verification.png
    :width: 80 %
    :align: center


