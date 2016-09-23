=========================
Automatic Differentiation
=========================

In the context of Newton's method, one must provide a Jacobian matrix that contains the
partial derivatives of some function with respect to an iterate


.. math::

    \text{Jacobian Matrix} = 
    \begin{pmatrix}
      \frac{\partial R_1}{\partial Q_1}     \quad   \frac{\partial R_1}{\partial Q_2}  \quad   \cdots \\
      \frac{\partial R_2}{\partial Q_1}     \quad   \frac{\partial R_2}{\partial Q_2}  \quad   \cdots \\
      \vdots
    \end{pmatrix} 
    


This presents a developer with a number of undesirable tasks:
    - One must derive a differentiated expression of the original function(R) with respect to the 
      iterate(Q); A tedious and error-prone process.
    - One must implement the differentiated expressions so they can be used by the solver. A tedious
      and error-prone process.
    - Any time an equation or spatial operator is modified or added, these processes must be 
      repeated and reverified for correctness.


**A solution** - Automatic Differentiation

Automatic Differentiation means many things because there are many approaches to this topic, but 
the end goal is that a developer offloads the tasks of deriving and implementing differentiated 
functions to some other tool, library, or preprocessor. Here, we will discuss the approach taken
in ChiDG and how it fits into the framework



--------------------------------------------
DNAD - Dual Number Automatic Differentiation
--------------------------------------------

ChiDG uses a tool previously developed by Wenbin Yu and others:

**References:**
``Yu, W. and Blair, M, "DNAD, a Simple Tool for Automatic Differentiation of Fortran Codes Using Dual Numbers"``
``Spall, R. and Yu, W., "Imbedded Dual-Number Automatic Differentiation for CFD Sensitivity Analysis"``

.. code-block:: Fortran

    type :: AD_D
        real(rk)              :: x_ad_      ! function value
        real(rk), allocatable :: xp_ad_(:)  ! derivatives
    end type AD_D












**References:**

[1]Yu, W. and Blair, M, "DNAD, a Simple Tool for Automatic Differentiation of Fortran Codes
Using Dual Numbers", Computer Physics Communications, Vol. 184, 2013, pp. 1446-1452

[2]Spall, R. and Yu, W., "Imbedded Dual-Number Automatic Differentiation for CFD
Sensitivity Analysis", Journal of Fluids Engineering, Vol 135, 2013





