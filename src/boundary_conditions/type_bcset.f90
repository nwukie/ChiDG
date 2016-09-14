module type_bcset
#include <messenger.h>
    use mod_kinds,              only: rk, ik

    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_bc,                only: bc_t
    use type_properties,        only: properties_t
    use type_bcset_coupling,    only: bcset_coupling_t
    implicit none


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

        type(bc_t),     allocatable     :: bcs(:)

    contains

        procedure   :: add     !< Call for adding a boundary condition
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
    function add(self,bc) result(BC_ID)
        class(bcset_t),     intent(inout)   :: self
        type(bc_t),         intent(inout)   :: bc

        integer(ik) :: ibc, ierr, BC_ID

        type(bc_t), allocatable :: temp_bcs(:)


        !
        ! Increment number of boundary conditions
        !
        self%nbcs = self%nbcs + 1


        !
        ! Set BC_ID
        !
        bc%BC_ID = self%nbcs
        BC_ID    = self%nbcs


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
        temp_bcs(self%nbcs) = bc


        !
        ! Move allocation to bcset storage
        !
        call move_alloc(temp_bcs,self%bcs)


    end function add
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
                !bcset_coupling%bc(ibc)%elems         = self%bcs(ibc)%bc%elems
                !bcset_coupling%bc(ibc)%coupled_elems = self%bcs(ibc)%bc%coupled_elems

                ! Only copy if there is data to copy. Or else, coupled_elements will not be allocated
                if (self%bcs(ibc)%bc_patch%ielement_l_%size() > 0) then
                    bcset_coupling%bc(ibc)%elems         = self%bcs(ibc)%bc_patch%ielement_l_%data()
                    bcset_coupling%bc(ibc)%coupled_elems = self%bcs(ibc)%bc_patch%coupled_elements
                end if
        end do ! ibc


    end function get_bcset_coupling
    !******************************************************************************************


















    subroutine destructor(self)
        type(bcset_t),  intent(inout) :: self

    end subroutine


end module type_bcset
