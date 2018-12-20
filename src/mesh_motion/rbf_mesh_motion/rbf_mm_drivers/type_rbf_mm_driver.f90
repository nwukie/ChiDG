! To perform mesh motion using RBFs, we need to supply the values of displacement and velocity
! at the RBF source nodes contained in an rbf_node_patch_mm_t. This might come from evaluating a
! specified prescribed mesh motion function, or in a practical FSI simulation, from an external
! structural solver. 
!
! To accommodate these different scenarios, we start here with an abstract derived type that
! performs the "update" task - that is, to update the dnodes and vnodes at RBF source nodes.
! Extensions of this abstract derive type will implement procedures to accomplish this in 
! different ways.
!
!
module type_rbf_mm_driver
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    implicit none

    type, public, abstract :: rbf_mm_driver_t


        character(len=:),       allocatable, private    :: name_

    contains

        procedure(init_interface),          deferred    :: init             !< Initialize function and register options
        procedure(compute_disp_interface),  deferred    :: compute_disp          !< Elemental function definition
        procedure(compute_vel_interface),   deferred    :: compute_vel          !< Elemental function definition

        procedure                                       :: set_name         !< Add function name. Only use for initialization
        procedure                                       :: get_name         !< Return the function name

    end type rbf_mm_driver_t

    abstract interface
        subroutine init_interface(self)
            import rbf_mm_driver_t 

            class(rbf_mm_driver_t),  intent(inout)  :: self

        end subroutine
    end interface


    abstract interface
        function compute_disp_interface(self, time, node) result(val)
            use mod_kinds,  only: rk
            import rbf_mm_driver_t

            class(rbf_mm_driver_t),  intent(inout)      :: self
            real(rk),                intent(in)          :: time
            real(rk),                intent(in)          :: node(3)

            real(rk)                                    :: val(3)

        end function 
    end interface

    abstract interface
        function compute_vel_interface(self, time, node) result(val)
            use mod_kinds,  only: rk
            import rbf_mm_driver_t

            class(rbf_mm_driver_t),  intent(inout)      :: self
            real(rk),                intent(in)          :: time
            real(rk),                intent(in)          :: node(3)

            real(rk)                                    :: val(3)

        end function

    end interface



contains

    !>  Add a name for the function
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set_name(self,fname)
        class(rbf_mm_driver_t),  intent(inout)   :: self
        character(*),       intent(in)      :: fname

        self%name_ = fname

    end subroutine set_name
    !********************************************************************************************


    !>  Return the name of the function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_name(self) result(fname)
        class(rbf_mm_driver_t),  intent(inout)   :: self

        character(len=:),   allocatable :: fname

        fname = self%name_

    end function get_name
    !*****************************************************************************************





end module type_rbf_mm_driver
