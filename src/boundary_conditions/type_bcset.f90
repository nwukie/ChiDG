module type_bcset
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_mesh,              only: mesh_t
    use type_equationset,       only: equationset_t
    use type_solverdata,        only: solverdata_t
    use type_bc,                only: bc_t
    use type_bcwrapper,         only: bcwrapper_t
    use type_properties,        only: properties_t
    use type_bcset_coupling,    only: bcset_coupling_t
    implicit none
    private


    !>  Type for a set of boundary conditions that get applied to a block
    !!      - Contains an array of wrapped boundary conditions types that
    !!        can be added.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public :: bcset_t

        integer(ik)                     :: nbcs = 0

        class(bcwrapper_t), allocatable :: bcs(:)       !< Array of boundary conditions. Using a wrapper here
                                                        !! because we can't allocate an array of polymorphic variables
                                                        !! without SOURCE or MOLD, for which we need a concrete type

    contains

        procedure   :: add     !< Call for adding a boundary condition
        procedure   :: apply   !< Spatial application of the boundary condition
        final       :: destructor

        procedure   :: get_bcset_coupling

    end type bcset_t
    !******************************************************************************************




contains


    !>  Add boundary condition to the boundary condition array
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[in]  bc  Boundary condition that is being added to the list
    !!  
    !------------------------------------------------------------------------------------------
    subroutine add(self,bc)
        class(bcset_t),     intent(inout)   :: self
        class(bc_t),        intent(in)      :: bc

        integer(ik) :: ibc, ierr

        type(bcwrapper_t), allocatable :: temp_bcs(:)


        !
        ! Increment number of boundary conditioners
        !
        self%nbcs = self%nbcs + 1


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_bcs(self%nbcs), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( self%nbcs > 1) then
            temp_bcs(1:size(self%bcs)) = self%bcs(1:size(self%bcs))
        end if


        !
        ! Allocate new boundary condition
        !
        allocate(temp_bcs(self%nbcs)%bc, source=bc)     ! I think this should source all of the data as well, like an assign


        !
        ! Move allocation to bcset storage
        !
        call move_alloc(temp_bcs,self%bcs)


    end subroutine add
    !******************************************************************************************








    !>  Call bc_t%apply for each boundary condition in the set
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[inout]   mesh    mesh_t containing elements and faces
    !!  @param[inout]   sdata   solverdata_t object containing solution, rhs, and linearization
    !!  @param[inout]   iblk    integer block direction with respect to which we are computing the linearization
    !!  @param[inout]   prop    properties_t object with equationset properties, and material_t objects
    !!
    !----------------------------------------------------------------------------------------
    subroutine apply(self,mesh,sdata,prop,idom)
        class(bcset_t),         intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh(:)
        class(solverdata_t),    intent(inout)   :: sdata
        class(properties_t),    intent(inout)   :: prop
        integer(ik),            intent(in)      :: idom

        integer(ik) :: ibc

        !
        ! Loop through boundary condition array and call apply for each
        !
        do ibc = 1,size(self%bcs)

            !
            ! Only apply if there is an allocated boundary condition in the current slot
            !
            if (allocated(self%bcs(ibc)%bc)) then
                call self%bcs(ibc)%bc%apply(mesh,sdata,prop,idom)
            end if

        end do

    end subroutine apply
    !******************************************************************************************












    !>  Return coupling information for the bcset_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_bcset_coupling(self) result(bcset_coupling)
        class(bcset_t), intent(in)   :: self

        type(bcset_coupling_t)  :: bcset_coupling
        integer(ik)             :: ibc, ierr


        !
        ! Allocate coupling storage for each bc in the set.
        ! 
        allocate(bcset_coupling%bc(self%nbcs), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Loop through bcs and assemble coupling information
        !
        do ibc = 1,self%nbcs
            if ( allocated(self%bcs(ibc)%bc) ) then

                bcset_coupling%bc(ibc)%elems         = self%bcs(ibc)%bc%elems
                bcset_coupling%bc(ibc)%coupled_elems = self%bcs(ibc)%bc%coupled_elems

            end if
        end do ! ibc


    end function get_bcset_coupling
    !******************************************************************************************


















    subroutine destructor(self)
        type(bcset_t),  intent(inout) :: self

    end subroutine


end module type_bcset
