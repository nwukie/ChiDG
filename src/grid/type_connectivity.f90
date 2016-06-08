module type_connectivity
#include <messenger.h>
    use mod_kinds,  only: ik
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: connectivity_t

        integer(ik)                 :: nnodes
        integer(ik),    allocatable :: data(:,:)

    contains
        
        procedure   :: init

        procedure   :: get_element_index        !< Return the index of a given element
        procedure   :: get_element_mapping      !< Return the mapping of a given element
        procedure   :: get_nelements            !< Return the number of elements in the domain
        procedure   :: get_nnodes               !< Return the number of nodes in the domain

    end type connectivity_t
    !*******************************************************************************


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self,nelements)
        class(connectivity_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: nelements
        
        integer(ik) :: ierr


        allocate(self%data(nelements,10), stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !********************************************************************************












    !>  Return the element index from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_index(self,index) result(ielem)
        class(connectivity_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: index

        integer(ik) :: ielem

        ielem = self%data(index,1)

    end function get_element_index
    !********************************************************************************









    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_mapping(self,index) result(mapping)
        class(connectivity_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: index

        integer(ik) :: mapping

        mapping = self%data(index,2)

    end function get_element_mapping
    !********************************************************************************






    !>  Return the number of elements in the connectivity structure.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_nelements(self) result(nelements)
        class(connectivity_t),  intent(in)  :: self

        integer(ik) :: nelements

        nelements = size(self%data,1)

    end function get_nelements
    !********************************************************************************






    !>  Return the number of nodes in the connectivity structure.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_nnodes(self) result(nnodes)
        class(connectivity_t),  intent(in)  :: self

        integer(ik) :: nnodes

        nnodes = self%nnodes

    end function get_nnodes
    !********************************************************************************





end module type_connectivity
