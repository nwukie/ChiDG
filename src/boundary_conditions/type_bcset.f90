module type_bcset
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_mesh,          only: mesh_t
    use atype_equationset,  only: equationset_t
    use type_solverdata,    only: solverdata_t
    use atype_bc,           only: bc_t
    use type_bcwrapper,     only: bcwrapper_t
    use type_properties,    only: properties_t
    implicit none
    private


    !> Type for a set of boundary conditions that get applied to a block
    !!      - Contains an array of wrapped boundary conditions types that
    !!        can be added.
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------
    type, public :: bcset_t

        class(bcwrapper_t), allocatable    :: bcs(:)       !< Array of boundary conditions. Using a wrapper here
                                                           !! because we can't allocate an array of polymorphic variables
                                                           !! without SOURCE or MOLD, for which we need a concrete type

    contains
        procedure :: init    !< Call for initializing storage for boundary conditions
        procedure :: add     !< Call for adding a boundary condition
        procedure :: apply   !< Spatial application of the boundary condition

    end type bcset_t




contains

    !> Initialize boundary condition routine
    !!      - Allocate default storage length for boundary condition slots
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(bcset_t),     intent(inout)       :: self

        integer(ik) :: nbc, ierr

        nbc = 6 ! Default number of bc's. Six sides to a block

        ! Allocate default number of boundary conditions
        allocate(self%bcs(nbc), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine
    






    !> Add boundary condition to the boundary condition array
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  bc  Boundary condition that is being added to the list
    !------------------------------------------------------------------------------------------
    subroutine add(self,bc)
        class(bcset_t),     intent(inout)   :: self
        class(bc_t),        intent(in)      :: bc

        integer(ik) :: ibc

        !
        ! Look for open boundary condition slot
        ! If the current bc is allocated then we go to the next one to look for an open slot
        !
        ibc = 1
        do while (allocated(self%bcs(ibc)%bc))
            ibc = ibc + 1

            !> Make sure we don't go over the bound of allocated slots
            if (ibc > size(self%bcs)) then
                call signal(FATAL,"bcset%add: Number of boundary conditions exceeds allocated slots")
            end if
        end do


        !
        ! Allocate concrete type, sourced from incoming bc
        !
        allocate(self%bcs(ibc)%bc, source=bc)   ! I think this should source all of the data as well, like an assign.

    end subroutine






    !>  Call bc_t%apply for each boundary condition in the set
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   mesh    mesh_t containing elements and faces
    !!  @param[inout]   sdata   solverdata_t object containing solution, rhs, and linearization
    !!  @param[inout]   iblk    integer block direction with respect to which we are computing the linearization
    !!  @param[inout]   prop    properties_t object with equationset properties, and material_t objects
    !----------------------------------------------------------------------------------------
    subroutine apply(self,mesh,sdata,idom,iblk,prop)
        class(bcset_t),         intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh(:)
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: idom
        integer(ik),            intent(in)      :: iblk
        class(properties_t),    intent(inout)   :: prop

        integer(ik) :: ibc

        !
        ! Loop through boundary condition array and call apply for each
        !
        do ibc = 1,size(self%bcs)

            !
            ! Only apply if there is an allocated boundary condition in the current slot
            !
            if (allocated(self%bcs(ibc)%bc)) then
                call self%bcs(ibc)%bc%apply(mesh,sdata,idom,iblk,prop)
            end if

        end do

    end subroutine
    !----------------------------------------------------------------------------------------







end module type_bcset
