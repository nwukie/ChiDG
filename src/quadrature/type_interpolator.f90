module type_interpolator
#include <messenger.h>


    !>  An object defining an interpolation from modal polynomials
    !!  to some node sets.
    !!
    !!  Examples would be an interpolation to quadrature nodes or to 
    !!  nodes for post-processing.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/25/2017
    !!
    !---------------------------------------------------------------------
    type, public :: interpolator_t

        real(rk),   allocatable :: nodes(:,:)   ! (nnodes,ncoords)

        real(rk),   allocatable :: val(:,:)     ! (nnodes,nmodes)
        real(rk),   allocatable :: ddxi(:,:)    ! (nnodes,nmodes)
        real(rk),   allocatable :: ddeta(:,:)   ! (nnodes,nmodes)
        real(rk),   allocatable :: ddzeta(:,:)  ! (nnodes,nmodes)

    contains

         

    end type interpolator_t
    !*********************************************************************




contains





end module type_interpolator
