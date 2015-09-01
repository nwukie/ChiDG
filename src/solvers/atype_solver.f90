module atype_solver
    use mod_kinds,          only: rk,ik
    use type_domain,        only: domain_t
    use atype_matrixsolver, only: matrixsolver_t
    implicit none


    !> solver abstract type definition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-----------------------------------------------------
    type, abstract, public  :: solver_t

        real(rk)                        :: testing
        logical                         :: solverInitialized = .false.


    contains
        ! Must define these procedures in the extended type
        procedure(init_interface),   deferred   :: init
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
        subroutine init_interface(self,domain,options)
            use type_domain,    only: domain_t
            use type_dict,      only: dict_t
            import solver_t
            class(solver_t),        intent(inout)   :: self
            type(domain_t),         intent(inout)   :: domain
            type(dict_t), optional, intent(inout)   :: options
        end subroutine
    end interface






    ! Interface for passing a domain_t type
    abstract interface
        subroutine data_interface(self,domain,matrixsolver)
            use type_domain,        only: domain_t
            use atype_matrixsolver, only: matrixsolver_t
            import solver_t
            class(solver_t),                 intent(inout)   :: self
            type(domain_t),                  intent(inout)   :: domain
            class(matrixsolver_t), optional, intent(inout)   :: matrixsolver
        end subroutine
    end interface




contains






end module atype_solver
