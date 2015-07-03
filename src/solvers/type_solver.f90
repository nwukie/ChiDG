module atype_solver
    use mod_kinds,      only: rk,ik
!    use type_dict,      only: dict_t
!    use type_mesh,      only: mesh_t
    use type_expansion, only: expansion_t
    implicit none

    type(expansion_t), allocatable :: q(:)


    !> solver abstract type definition
    !!
    !-----------------------------------------------------
    type, abstract, public  :: solver_t
!        type(dict_t)    :: options

    contains
        ! Must define these procedures in the extended type
        procedure(self_interface),   deferred  :: init
!        procedure(domain_interface), deferred  :: solve

    end type solver_t

    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================
!    abstract interface
!        subroutine  domain_interface(self,domain)
!            use type_domain,  only: domain_t
!            import solver_t
!            class(solver_t), intent(inout)    :: self
!            class(domain_t), intent(inout)    :: domain
!        end subroutine
!    end interface

    abstract interface
        subroutine self_interface(self)
            import solver_t
            class(solver_t), intent(inout) :: self
        end subroutine
    end interface


contains





end module atype_solver
