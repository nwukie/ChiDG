module type_spatial_discretization
#include <messenger.h>
    implicit none



    !>  A class to handle the implementation of different spatial discretization loops.
    !!
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public, abstract :: spatial_discretization_t


    contains

        procedure(compute_interface), deferred   :: compute

    end type spatial_discretization_t
    !***************************************************************************************



    !>  Declaration of what the routine 'compute' must look like.
    !!
    !!
    !--------------------------------------------------------------------------------------
    abstract interface
        subroutine compute_interface(self,data,timing,info)
            use type_chidg_data,    only: chidg_data_t
            use mod_kinds,          only: rk, ik
            import spatial_discretization_t
            class(spatial_discretization_t),    intent(in)              :: self
            type(chidg_data_t),                 intent(inout), target   :: data
            real(rk),                           optional                :: timing
            integer(ik),                        optional                :: info
        end subroutine
    end interface
    !***************************************************************************************




contains



end module type_spatial_discretization
