=========================
Fluid boundary conditions
=========================


-------------
Inlet - Total
-------------

Setting ``Inlet - Total`` imposes total conditions and flow normal vector at an inlet boundary.


=====================   =======================================================================
Parameter               Description
=====================   =======================================================================
``Total Pressure``      Set ``Total Presure`` :math:`p_0` at the inlet boundary
``Total Temperature``   Set ``Total Temperature`` :math:`T_0` at the inlet boundary
``Normal-1``            Set first component of flow normal vector :math:`\vec{n}`
``Normal-2``            Set second component of flow normal vector :math:`\vec{n}`
``Normal-3``            Set third component of flow normal vector :math:`\vec{n}`
=====================   =======================================================================



Momentum is extrapolated from the interior and used to compute velocity, :math:`\vec{V}`. 
With the user-specified normal and Total quantities, static quantities to compute the 
boundary condition solution state :math:`Q_{bc}` are computed as:

.. math::

    T_{bc} &= T_0 -  \frac{|\vec{V}_-|}{2 c_p}

    p_{bc} &= p_0\bigg( \frac{T_{bc}}{T_0} \bigg)^{\frac{\gamma}{\gamma - 1}}

    \vec{V}_{bc} &= \vec{V}_- \cdot \vec{n}


With these static quantities, the boundary condition state is computed as:

.. math::

        \rho_{bc} &= \frac{p_{bc}}{R T_{bc}} 

        \vec{\rho V}_{bc} &= \rho_{bc} \vec{V}_{bc} 

        \rho E_{bc} &= \frac{p_{bc}}{\gamma - 1} + \frac{\rho}{2}( \vec{V}_{bc} \cdot \vec{V}_{bc} )

The gradient of the boundary state is extrapolated:

.. math::

    \nabla Q_{bc} = \nabla Q_{-}




--------------------------
Outlet - Constant Pressure
--------------------------

Setting ``Outlet - Constant Pressure`` imposes constant static pressure on an outlet boundary.

=====================   =======================================================================
Parameter               Description
=====================   =======================================================================
``Static Pressure``     Set ``Static Presure`` :math:`p` at the outlet boundary
=====================   =======================================================================


Density and Momentum are extrapolated from the interior state. The boundary condition state
is computed as:

.. math::

    \rho_{bc} &= \rho_{-}

    \vec{\rho V}_{bc} &= \vec{\rho V}_{-}

    \rho E_{bc} &= \frac{p}{\gamma - 1} + \frac{\rho_{bc}}{2}( \vec{V}_{-} \cdot \vec{V}_{-} )


The gradient of the boundary state is extrapolated:

.. math::

    \nabla Q_{bc} = \nabla Q_{-}



----
Wall
----

A ``Wall`` state function sets the normal component of momentum to zero and 
subtracts its kinetic energy contribution from the extrapolated energy from the interior.
For weakly imposed boundary conditions for the ``Euler`` equations, this defines a Slip Wall. 
For weakly imposed boundary conditions For the ``Navier-Stokes`` equations, 
this defines a No-Slip Wall. The normal flux for either case 

Heat transfer is defined to be adiabatic, :math:`\nabla T = 0`.


=====================   =======================================================================
Parameter               Description
=====================   =======================================================================
no parameters
=====================   =======================================================================


The normal component of the advective flux is given by:

.. math::

    \vec{F}^a(Q) \cdot \vec{n} = 
    \begin{pmatrix}
        \vec{\rho V} \cdot \vec{n} \\
        \rho u \vec{V} \cdot \vec{n} + p n_1 \\
        \rho v \vec{V} \cdot \vec{n} + p n_2 \\
        \rho w \vec{V} \cdot \vec{n} + p n_3 \\
        \rho H \vec{V} \cdot \vec{n}
    \end{pmatrix}


A slip wall is defined as :math:`\vec{V}_{bc} \cdot \vec{n} = 0`. A no-slip wall is defined as
:math:`\vec{V}_{bc} = 0`.

Imposing either of these conditions on the normal flux yields the same result:

.. math::

    
    \vec{F}^a(Q) \cdot \vec{n} = 
    \begin{pmatrix}
        0 \\
        p_- n_1 \\
        p_- n_2 \\
        p_- n_3 \\
        0
    \end{pmatrix}


So, for the Euler equations imposing :math:`\vec{V}_{bc} = 0` yields a slip wall. 
For the Navier-Stokes equations, imposing :math:`\vec{V}_{bc} = 0` yields a no-slip wall.

The boundary solution state is computed by extrapolating density, setting momentum to 
zero, and subtracting the lost kinetic energy from the interior energy:

.. math::

    Q_{bc} = 
    \begin{pmatrix}
        \rho_{-} \\
        0 \\ 
        0 \\ 
        0 \\ 
        \rho E_{-} - \frac{\rho_{-}}{2}(\vec{V}_{-} \cdot \vec{V}_{-})
    \end{pmatrix}


The gradient of the boundary state is computed using the no-slip condition and
also the adiabatic condition :math:`\nabla T = 0`.

If the temperature is a function of the primary fields :math:`T = T(\rho,\vec{\rho V},\rho E)`, 
then the gradient of temperature can be computed using the Chain Rule as:

.. math::

    \nabla T = \frac{\partial T}{\partial \rho} \nabla (\rho)  +  \frac{\partial T}{\partial \rho u} \nabla (\rho u)  +  \frac{\partial T}{\partial \rho v} \nabla(\rho v)  +  \frac{\partial T}{\partial \rho w}\nabla(\rho w)  +  \frac{\partial T}{\partial \rho E} \nabla(\rho E)

One can determine that the jacobian of temperature with respect to momentum goes to zero 
with velocity, :math:`\frac{\partial T}{\vec{\rho V}} \rightarrow 0` 
as :math:`V \rightarrow 0`. The adiabatic condition can then be imposed by setting gradients
of density and energy to zero:

.. math::

    \nabla Q_{bc} = 
    \begin{pmatrix}
        0 \\
        \nabla \vec{\rho V}_{-} \\
        0 \\
    \end{pmatrix}














--------
Symmetry
--------

A ``Symmetry`` condition defines the velocity gradient normal to the boundary as zero.
However, in contrast to the ``Wall`` boundary condition implementation where the 
condition is imposed on the normal component of the flux, here the normal component
of momentum is mirrored about the boundary face such that the normal component
of momentum is effectively zero.

The consequence of this implementation is that is effectively a slip-wall condition
for both the ``Euler`` and ``Navier-Stokes`` equations.

=====================   =======================================================================
Parameter               Description
=====================   =======================================================================
no parameters
=====================   =======================================================================


First, the normal component of momentum is computed. It is then mirrored about the
boundary face such that the normal momentum is zero. The boundary condition state is 
computed as:

.. math::

    Q_{bc} = 
    \begin{pmatrix}
        \rho_{-} \\
        \vec{\rho V}_{bc} - 2(\vec{\rho V}_{-} \cdot \vec{n})\vec{n} \\
        \rho E_{-}
    \end{pmatrix}


The gradient of the boundary state is extrapolated:

.. math::

    \nabla Q_{bc} = \nabla Q_{-}






--------
Farfield
--------

Specify static inflow/outflow conditions at a farfield boundary.


=====================   =======================================================================
Parameter               Description
=====================   =======================================================================
``Density``             Set static ``Density`` :math:`\rho` at the farfield boundary
``Pressure``            Set static ``Pressure`` :math:`p` at the farfield boundary
``Velocity-1``          Set first component of flow velocity :math:`\vec{V}`
``Velocity-2``          Set second component of flow velocity :math:`\vec{V}`
``Velocity-3``          Set third component of flow velocity :math:`\vec{V}`
=====================   =======================================================================


To-be-continued...












