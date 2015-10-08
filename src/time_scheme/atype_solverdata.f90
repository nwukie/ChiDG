module atype_solverdata
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    !use type_expansion,         only: expansion_t
    use type_blockmatrix,       only: blockmatrix_t
    use type_blockvector,       only: blockvector_t
    use type_mesh,              only: mesh_t
    implicit none


    !> solver abstract type definition
    !!
    !-----------------------------------------------------
    type, abstract, public  :: solverdata_t

        !>  Base solver data
        type(blockvector_t)             :: q                        !< Solution vector
        type(blockvector_t)             :: dq                       !< Change in solution vector
        type(blockvector_t)             :: rhs                      !< Spatial residual vector
        type(blockmatrix_t)             :: lin                      !< Linearization of the spatial scheme

        real(rk), allocatable           :: dt(:)                    !< Element-local time step



        ! Matrix views of expansion storage
        !type(expansion_t),  pointer     :: q_m(:,:,:)   => null()   !< Matrix view of solution expansions
        !type(expansion_t),  pointer     :: dq_m(:,:,:)  => null()   !< Matrix view of change in solution expansions
        !type(expansion_t),  pointer     :: rhs_m(:,:,:) => null()   !< Matrix view of rhs expansion



        logical                         :: solverInitialized = .false.


    contains
        ! Must define these procedures in the extended type
        procedure                               :: init_base    !> Initialize base required data structures. Should be called by specialized init procedures
        procedure(data_interface),   deferred   :: init

    end type solverdata_t

    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================

    abstract interface
        subroutine self_interface(self)
            import solverdata_t
            class(solverdata_t), intent(inout) :: self
        end subroutine
    end interface




    abstract interface
        subroutine data_interface(self,mesh)
            use type_mesh,  only: mesh_t
            import solverdata_t
            class(solverdata_t), intent(inout)    :: self
            type(mesh_t),        intent(in)       :: mesh
        end subroutine
    end interface




contains


    !>  Initialize solver base data structures
    !!      - allocate and initialize q, dq, rhs, and linearization.
    !!      - Should be called by specialized 'init' procedure for derived solvers.
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  mesh    Mesh definition which defines storage requirements
    !----------------------------------------------------------------------------------------------------------
    subroutine init_base(self,mesh)
        class(solverdata_t),    intent(inout), target   :: self
        type(mesh_t),           intent(in)              :: mesh

        integer(ik)                                 :: nterms_s, ielem, nelem, neqns, ierr
        !type(expansion_t), pointer                  :: temp(:)


        nterms_s = mesh%nterms_s
        nelem    = mesh%nelem
        neqns    = mesh%neqns


        !
        ! Initialize and allocate storage
        !
        call self%q%init(  mesh)
        call self%dq%init( mesh)
        call self%rhs%init(mesh)
        call self%lin%init(mesh,'full')

        allocate(self%dt(nelem),stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Initialize solution matrix view
        !
!        temp => self%q
!        self%q_m(1:mesh%nelem_xi, 1:mesh%nelem_eta, 1:mesh%nelem_zeta) => temp(1:mesh%nelem)
!
!        temp => self%dq
!        self%dq_m(1:mesh%nelem_xi, 1:mesh%nelem_eta, 1:mesh%nelem_zeta) => temp(1:mesh%nelem)
!
!        temp => self%rhs
!        self%rhs_m(1:mesh%nelem_xi, 1:mesh%nelem_eta, 1:mesh%nelem_zeta) => temp(1:mesh%nelem)


        !
        ! Confirm solver initialization
        !
        self%solverInitialized = .true.

    end subroutine







end module atype_solverdata
