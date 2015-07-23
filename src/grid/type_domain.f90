module type_domain
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: HALF, TWO, ONE, DIAG, NFACES, &
                                  XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use mod_equations,      only: AssignEquationSet

    use type_mesh,          only: mesh_t
    use type_point,         only: point_t
    use type_expansion,     only: expansion_t
    use atype_equationset,  only: equationset_t

    implicit none

    private


    !> Domain data type
    !   - contains mesh, solution, and equation set information
    !
    !   @author Nathan A. Wukie
    !---------------------------------------------------------------------------------------
    type, public :: domain_t
        character(100)                    :: name                   !> Domain name -- not currently used
        type(mesh_t)                      :: mesh                   !> Mesh storage
        class(equationset_t), allocatable :: eqnset                 !> Equation set solved on this domain
        type(expansion_t),    allocatable :: q(:)                   !> Array of solution expansions. One for each element.
        type(expansion_t),    pointer     :: q_m(:,:,:) => null()   !> Matrix view of solution expansions

        logical                             :: geomInitialized = .false.
        logical                             :: numInitialized  = .false.

    contains
        procedure       :: init_geom
        procedure       :: init_sol
        final           :: destructor

    end type domain_t
    !---------------------------------------------------------------------------------------


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
        type(expansion_t), pointer                :: temp(:)

        integer(ik) :: ielem

        if (self%numInitialized) stop "Error: domain%init_sol -- domain numerics already initialized"


        ! Set domain equation set
        call AssignEquationSet(eqnstring,self%eqnset)               !> Factory method for allocating an equation set

        ! Initialize mesh solution data
        call self%mesh%init_sol(self%eqnset%neqns,nterms_s)

        ! Initialize solution
        allocate(self%q(self%mesh%nelem))                           !> Allocate an expansion type for each element
        do ielem = 1,self%mesh%nelem
            call self%q(ielem)%init(nterms_s,self%eqnset%neqns)     !> Initialize expansion for each element
        end do


        ! Initialize solution matrix view
        temp => self%q
        self%q_m(1:self%mesh%nelem_xi, 1:self%mesh%nelem_eta, 1:self%mesh%nelem_zeta) => temp(1:self%mesh%nelem)


        self%numInitialized = .true.
    end subroutine





    !>  Destructor
    !!      -   if allocatable components are allocated, call deallocation routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------
    subroutine destructor(self)
        type(domain_t), intent(inout) :: self

!> Shouldn't need to deallocate an 'allocatable'. The compiler is supposed to do that for you
!        if (allocated(self%eqnset)) deallocate(self%eqnset)
!        if (allocated(self%q))      deallocate(self%q)
    end subroutine

end module type_domain
