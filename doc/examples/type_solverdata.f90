module type_solverdata
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_chidgVector,       only: chidgVector_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_mesh,              only: mesh_t
    implicit none


    !> solver type definition
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    !> [solverdata_t]
    type, public  :: solverdata_t

        type(chidgVector_t)             :: q                        ! Solution vector
        type(chidgVector_t)             :: dq                       ! Change in solution vector
        type(chidgVector_t)             :: rhs                      ! Residual of the spatial scheme
        type(chidgMatrix_t)             :: lhs                      ! Linearization of the spatial scheme

    contains
        generic, public       :: init => init_base
        procedure, private    :: init_base

    end type solverdata_t
    !> [solverdata_t]
    !-------------------------------------------------------------------------------------------------------


contains


    !>  Initialize solver base data structures
    !!      - allocate and initialize q, dq, rhs, and linearization.
    !!      - Should be called by specialized 'init' procedure for derived solvers.
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  mesh    Mesh definition which defines storage requirements
    !----------------------------------------------------------------------------------------------------------
    subroutine init_base(self,mesh)
        class(solverdata_t),     intent(inout), target   :: self
        type(mesh_t),            intent(in)              :: mesh(:)

        integer(ik) :: nterms_s, ielem, nelem, neqns, ierr, ndom, maxelems, idom
        logical     :: increase_maxelems = .false.


        !
        ! Initialize and allocate storage
        !
        call self%q%init(  mesh)
        call self%dq%init( mesh)
        call self%rhs%init(mesh)
        call self%lhs%init(mesh,'full')


    
        !
        ! Find maximum number of elements in any domain
        !
        ndom = size(mesh)
        maxelems = 0
        do idom = 1,ndom

            increase_maxelems = ( mesh(idom)%nelem > maxelems )

            if (increase_maxelems) then
                maxelems = mesh(idom)%nelem
            end if

        end do


        !
        ! Allocate timestep storage
        !
        allocate(self%dt(ndom,maxelems),stat=ierr)
        if (ierr /= 0) call AllocationError

        
        !
        ! Confirm solver initialization
        !
        self%solverInitialized = .true.

    end subroutine







end module type_solverdata
