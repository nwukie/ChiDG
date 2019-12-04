module type_system_assembler
#include <messenger.h>
    implicit none





    !>  Object for implementing an algorithm for assembling the equation system, lhs, rhs,
    !!  including contributions from both the spatial and temporal terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2017
    !!
    !----------------------------------------------------------------------------------------
    type, abstract, public :: system_assembler_t


    contains

        procedure(assemble_interface), deferred :: assemble

    end type system_assembler_t
    !****************************************************************************************




    ! Interface for passing a domain_t type
    abstract interface
        subroutine assemble_interface(self,data,differentiate,timing)
            use mod_kinds,              only: rk, ik
            use type_chidg_data,        only: chidg_data_t
            import system_assembler_t
            class(system_assembler_t),  intent(inout)               :: self
            type(chidg_data_t),         intent(inout)               :: data
            integer(ik),                intent(in)                  :: differentiate
            real(rk),                   intent(inout),  optional    :: timing
        end subroutine assemble_interface
    end interface


end module type_system_assembler
