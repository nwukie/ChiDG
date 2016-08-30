module type_boundary_operator
    use mod_kinds,          only: rk, ik
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), abstract, public :: boundary_operator_t

    contains

        procedure(compute_interface), deferred :: compute

    end type boundary_operator_t
    !**********************************************************************************




    abstract interface
        subroutine compute_interface(self,worker,prop)
            use mod_kinds,  only: ik
            import boundary_operator_t
            import chidg_worker_t
            import properties_t

            class(boundary_operator_t), intent(in)      :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
        end subroutine
    end interface




contains








end module type_boundary_operator
