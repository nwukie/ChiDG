module type_domain
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_equations,      only: AssignEquationSet
    use mod_solver,         only: CreateSolver

    use type_mesh,          only: mesh_t
    use atype_solver,       only: solver_t
    use type_point,         only: point_t
    use atype_equationset,  only: equationset_t

    implicit none

    private


    !> Domain data type
    !!      - contains mesh, solution, and equation set information
    !!
    !!   @author Nathan A. Wukie
    !-------------------------------------------------------------------------------------------------
    type, public :: domain_t
        character(100)                      :: name                       !> Domain name -- not currently used
        type(mesh_t)                        :: mesh                       !> Mesh storage
        class(solver_t),      allocatable   :: solver                     !> Solver storage
        class(equationset_t), allocatable   :: eqnset                     !> Equation set solved on this domain

        logical                             :: geomInitialized = .false.
        logical                             :: numInitialized  = .false.

    contains
        procedure       :: init_geom
        procedure       :: init_sol
        final           :: destructor

    end type domain_t
    !-------------------------------------------------------------------------------------------------


contains
    
    !>  Initialize domain geometry
    !!      - call geometry initialization for mesh component
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  nterms_c    Number of terms in the modal representation of the cartesian coordinates
    !!  @param[in]  points      Array of cartesian points defining the element
    !------------------------------------------------------------------
    subroutine init_geom(self,nterms_c,points)
        class(domain_t),    intent(inout)   :: self
        integer(ik),        intent(in)      :: nterms_c
        type(point_t),      intent(in)      :: points(:,:,:)

        if (self%geomInitialized) stop "Error: domain%init_geom -- domain geometry already initialized"

        ! Initialize mesh geometry
        call self%mesh%init_geom(nterms_c,points)

        self%geomInitialized = .true.
    end subroutine





    !>  Initialize domain numerics
    !!      -   call routine to initialize and assign equation set
    !!      -   call numerics initialization for mesh component
    !!      -   allocate and initialize solution storage
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  eqnstring   Character string specifying the equation set being solved
    !!  @param[in]  nterms_s    Number of terms in the modal representation of the solution
    !------------------------------------------------------------------
    subroutine init_sol(self,eqnstring,nterms_s)
        class(domain_t),    intent(inout), target :: self
        character(*),       intent(in)            :: eqnstring
        integer(ik),        intent(in)            :: nterms_s


        associate ( mesh => self%mesh, solver => self%solver, eqnset => self%eqnset)

            if (self%numInitialized) call signal(FATAL,'domain%init_sol -- domain numerics already initialized')

            ! Call factory methods for equationset and solver
            call AssignEquationSet(eqnstring,self%eqnset)          !> Factory method for allocating a equation set
            call CreateSolver('default',self%solver)               !> Factory method for allocating a solver

            ! Initialize mesh solution data
            call self%mesh%init_sol(eqnset%neqns,nterms_s)

            !> Call initialization for solver and solver data
            call solver%init(mesh)

            self%numInitialized = .true.

        end associate
    end subroutine





    !>  Destructor
    !!      -   if allocatable components are allocated, call deallocation routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------
    subroutine destructor(self)
        type(domain_t), intent(inout) :: self


    end subroutine

end module type_domain
