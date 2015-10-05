module atype_time_scheme
    use mod_kinds,          only: rk,ik
    use type_domain,        only: domain_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_dict,          only: dict_t
    use type_timer,         only: timer_t
    use type_rvector,       only: rvector_t
    use type_ivector,       only: ivector_t
    implicit none


    !> solver abstract type definition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------
    type, abstract, public  :: time_scheme_t

        real(rk)                        :: testing
        logical                         :: solverInitialized = .false.



        ! OPTIONS
        real(rk)        :: dt = 0.001_rk        !< Time-step increment
        real(rk)        :: tol = 1.e-13_rk      !< Convergence tolerance
        integer(ik)     :: nsteps = 100         !< Number of time steps to compute
        integer(ik)     :: nwrite = 10          !< Write data every 'nwrite' steps


        type(timer_t)   :: timer                !< Timer data-type


        ! Data logs
        type(rvector_t)                 :: residual_L2norm
        type(ivector_t)                 :: nnewton_iterations
        type(rvector_t)                 :: iteration_time


    contains
        procedure   :: init
        procedure   :: init_spec

        ! Must define these procedures in the extended type
        procedure(data_interface),   deferred   :: solve

        procedure   :: set
    end type time_scheme_t
    !------------------------------------------------------------------------------









    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================
    
    abstract interface
        subroutine self_interface(self)
            import time_scheme_t
            class(time_scheme_t),   intent(inout)   :: self
        end subroutine
    end interface



    abstract interface
        subroutine init_interface(self,domain,options)
            use type_domain,    only: domain_t
            use type_dict,      only: dict_t
            import time_scheme_t
            class(time_scheme_t),   intent(inout)   :: self
            type(domain_t),         intent(inout)   :: domain
            type(dict_t), optional, intent(inout)   :: options
        end subroutine
    end interface





    ! Interface for passing a domain_t type
    abstract interface
        subroutine data_interface(self,domain,matrixsolver)
            use type_domain,        only: domain_t
            use atype_matrixsolver, only: matrixsolver_t
            import time_scheme_t
            class(time_scheme_t),                 intent(inout)   :: self
            type(domain_t),                  intent(inout)   :: domain
            class(matrixsolver_t), optional, intent(inout)   :: matrixsolver
        end subroutine
    end interface




contains

    !> Common time_scheme initialization interface.
    !!      - Call initialization for options if present
    !!      - Call user-specified initialization routine for concrete type
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   domain      Domain type
    !!  @param[inout]   options     Dictionary containing options
    !----------------------------------------------------------
    subroutine init(self,domain)
        class(time_scheme_t),   intent(inout)   :: self
        type(domain_t),         intent(inout)   :: domain


        !
        ! Call user-specified initialization
        !
        call self%init_spec(domain)



        self%solverInitialized = .true.
    end subroutine






    !> Procedure for setting base time_scheme options
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  options     Dictionary containing base solver options
    !!
    !-----------------------------------------------------------
    subroutine set(self,options)
        class(time_scheme_t),   intent(inout)   :: self
        type(dict_t),           intent(inout)   :: options


        call options%get('dt',self%dt)
        call options%get('tol',self%tol)
        call options%get('nsteps',self%nsteps)
        call options%get('nwrite',self%nwrite)
    end subroutine







    !> Default blank initialization-specialization routine.
    !! This can be overwritten with specific instructions for a conrete
    !! time_scheme.
    !!
    !!  @author Nathan A. Wukie
    !!
    !--------------------------------------------------------------------
    subroutine init_spec(self,domain,options)
        class(time_scheme_t),   intent(inout)   :: self
        type(domain_t),         intent(inout)   :: domain
        type(dict_t), optional, intent(inout)   :: options



    end subroutine














end module atype_time_scheme
