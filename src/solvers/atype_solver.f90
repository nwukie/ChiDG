module atype_solver
    use mod_kinds,          only: rk,ik
    use type_domain,        only: domain_t
    implicit none


    !> solver abstract type definition
    !!
    !-----------------------------------------------------
    type, abstract, public  :: solver_t

        real(rk)                        :: testing
        logical                         :: solverInitialized = .false.


    contains
        ! Must define these procedures in the extended type
        procedure(data_interface),   deferred   :: init
        procedure(data_interface),   deferred   :: solve

    end type solver_t

    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================

    abstract interface
        subroutine self_interface(self)
            import solver_t
            class(solver_t),    intent(inout)   :: self
        end subroutine
    end interface




    abstract interface
        subroutine data_interface(self,domain)
            use type_domain,  only: domain_t
            import solver_t
            class(solver_t),    intent(inout)   :: self
            type(domain_t),     intent(inout)   :: domain
        end subroutine
    end interface




contains






end module atype_solver
