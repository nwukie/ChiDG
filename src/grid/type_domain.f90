module type_domain
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_equations,      only: create_equationset
    use mod_solverdata,     only: create_solverdata
    use mod_bc,             only: create_bc

    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use atype_solverdata,   only: solverdata_t
    use atype_equationset,  only: equationset_t
    use atype_bc,           only: bc_t
    use type_bcset,         only: bcset_t
    use type_dict,          only: dict_t
    implicit none

    private


    !> Domain data type
    !!      - contains mesh, solution, equation set, and boundary condition information
    !!
    !!   @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------------
    type, public :: domain_t
        character(100)                      :: name                       !< Domain name -- not currently used
        type(mesh_t)                        :: mesh                       !< Mesh storage
        type(bcset_t)                       :: bcset                      !< Boundary condition set
        class(equationset_t), allocatable   :: eqnset                     !< Equation set solved on this domain
        class(solverdata_t),  allocatable   :: sdata                      !< Solver data storage

        logical                             :: geomInitialized = .false.
        logical                             :: numInitialized  = .false.

    contains
        procedure       :: init_geom
        procedure       :: init_sol
        procedure       :: init_bcs
        final           :: destructor

    end type domain_t
    !-------------------------------------------------------------------------------------------------


contains
    
    !>  Initialize domain geometry
    !!      - call geometry initialization for mesh component
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  nterms_c    Number of terms in the modal representation of the cartesian coordinates
    !!  @param[in]  points      Array of cartesian points defining the element
    !------------------------------------------------------------------
    subroutine init_geom(self,nterms_c,points)
        class(domain_t),    intent(inout)   :: self
        integer(ik),        intent(in)      :: nterms_c
        type(point_t),      intent(in)      :: points(:,:,:)

        if (self%geomInitialized) stop "Error: domain%init_geom -- domain geometry already initialized"


        call self%mesh%init_geom(nterms_c,points)   ! Initialize mesh geometry
        call self%bcset%init()                      ! Initialize boundary condition storage


        self%geomInitialized = .true.
        
    end subroutine init_geom







    !> Initialize domain boundary conditions
    !!      - Allocate correct boundary condition instances
    !!      - Call boundary condition initialization routines
    !!      - Add boundary condition to self%bcset
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------
    subroutine init_bcs(self,bcstr,iface,options)
        class(domain_t),            intent(inout)   :: self
        character(*),               intent(in)      :: bcstr
        integer(ik),                intent(in)      :: iface
        type(dict_t), optional,     intent(in)      :: options

        class(bc_t), allocatable    :: bc
        integer(ik) :: ibc, nbc


        !
        ! Boundary condition factory for dynamic bc_t allocation
        !
        call create_bc(bcstr,bc)

        !
        ! Call initialization for boundary condition
        ! 
        call bc%init(self%mesh,iface,options)   ! It should be okay to pass the optional argument, as long as it is optional in bc%init

        !
        ! Add initialized boundary condition to bc-set
        !
        call self%bcset%add(bc)


    end subroutine init_bcs







    !>  Initialize domain numerics
    !!      -   call routine to initialize and assign equation set
    !!      -   call numerics initialization for mesh component
    !!      -   allocate and initialize solution storage
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  eqnstring   Character string specifying the equation set being solved
    !!  @param[in]  nterms_s    Number of terms in the modal representation of the solution
    !------------------------------------------------------------------
    subroutine init_sol(self,eqnstring,nterms_s)
        class(domain_t),    intent(inout)   :: self
        character(*),       intent(in)      :: eqnstring
        integer(ik),        intent(in)      :: nterms_s


        if (self%numInitialized) call signal(FATAL,'domain%init_sol -- domain numerics already initialized')


        !
        ! Call factory methods for equationset and solver       
        !
        call create_equationset(eqnstring,self%eqnset)      ! Factory method for allocating a equation set
        call create_solverdata('base',self%sdata)           ! Factory method for allocating solverdata. Allocate base data

        !
        ! Initialize mesh numerics, solution data, and boundary conditions
        !
        call self%mesh%init_sol(self%eqnset%neqns,nterms_s) ! Call initialization for mesh required in solution procedure
        call self%sdata%init(self%mesh)                     ! Call initialization for solver and solver data
        !call self%init_bcs()
        
        
        self%numInitialized = .true.                        ! Confirm initialization
    end subroutine init_sol















    !>  Destructor
    !!      -   if data is allocated via pointers, call deallocation 
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------
    subroutine destructor(self)
        type(domain_t), intent(inout) :: self


    end subroutine

end module type_domain
