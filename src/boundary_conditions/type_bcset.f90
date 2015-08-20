module type_bcset
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_mesh,          only: mesh_t
    use atype_equationset,  only: equationset_t
    use atype_solverdata,   only: solverdata_t
    use atype_bc,           only: bc_t
    use type_bcwrapper,     only: bcwrapper_t
    implicit none
    private


    !> Abstract base-type for boundary condition set
    !!
    !!
    !!
    !-------------------------------------------------
    type, public :: bcset_t
!        private

        class(bcwrapper_t), allocatable    :: bcs(:)       !> Array of boundary conditions. Using a wrapper here
                                                           !! because we can't allocate an array of polymorphic variables
                                                           !! without SOURCE or MOLD, for which we need a concrete type

    contains
        procedure :: init
        procedure :: apply              !> Spatial application of the boundary condition
        procedure :: add                !> Call for adding a boundary condition

    end type bcset_t




contains

    !> Initialize boundary condition routine
    !!      - Allocate default storage length for boundary condition slots
    !!
    !!  @author Nathan A. Wukie
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(bcset_t),     intent(inout)       :: self

        integer(ik) :: nbc, ierr  !> Default number of boundary conditions

        nbc = 6 !> Six sides to a block

        !> Allocate default number of boundary conditions
        allocate(self%bcs(nbc), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine
    



    !>  Call bc_t%apply for each boundary condition in the set
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   eqnset  Equation set that applies to the boundary condition
    !!  @param[inout]   mesh    Mesh structure containint elements and faces
    !!  @param[inout]   sdata   Solver data structure containing solution, rhs, and linearization
    !!  @param[inout]   iblk    Block direction with respect to which we are computing the linearization
    !-------------------------------------------------------------
    subroutine apply(self,eqnset,mesh,sdata,iblk)
        class(bcset_t),             intent(inout)   :: self
        class(equationset_t),       intent(inout)   :: eqnset
        type(mesh_t),               intent(inout)   :: mesh
        class(solverdata_t),        intent(inout)   :: sdata
        integer(ik),                intent(inout)   :: iblk

        integer(ik) :: ibc


        !> Loop through boundary condition array and call apply for each
        do ibc = 1,size(self%bcs)
            call self%bcs(ibc)%bc%apply(eqnset,mesh,sdata,iblk)
        end do


    end subroutine



    !> Add boundary condition to the boundary condition array
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  bc  Boundary condition that is being added to the list
    !--------------------------------------------------------------
    subroutine add(self,bc)
        class(bcset_t),     intent(inout)   :: self
        class(bc_t),        intent(in)      :: bc

        integer(ik) :: ibc


        !> Look for open boundary condition slot
        ibc = 1
        do while (self%bcs(ibc)%bc%isInitialized)
            ibc = ibc + 1

            !> Make sure we don't go over the bound of allocated slots
            if (ibc > size(self%bcs)) then
                call signal(FATAL,"bcset%add: Number of boundary conditions exceeds allocated slots")
            end if
        end do


        !> Allocate concrete type and assign bc to bcs(ibc)
        allocate(self%bcs(ibc)%bc, source=bc)   !> I think this should source all of the data as well, like an assign.
!        self%bcs(ibc) = bc

    end subroutine



    !------------------------------------------------------------------------------------------------

end module type_bcset
