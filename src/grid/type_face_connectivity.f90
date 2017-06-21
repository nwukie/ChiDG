module type_face_connectivity
#include <messenger.h>
    use mod_kinds,                  only: ik
    use mod_constants,              only: NO_PARTITION
    implicit none



    !>  This data type contains the connectivity information for a single face
    !!
    !!  Connectivity data format is as follows:
    !!  data = [ ipt1, ipt2, ipt3, ... ]
    !!
    !!  NOTE: no header/auxiliary data
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: face_connectivity_t

        integer(ik)                 :: header_size = 0  !< No leading auxiliary data
        integer(ik)                 :: mapping
        integer(ik),    allocatable :: data(:)

    contains
        
        procedure   :: init


        procedure   :: set_face_mapping     !< Set the mapping for the face
        procedure   :: set_face_nodes        !< Set the nodes indices defining the face connectivity

        procedure   :: get_face_mapping      !< Return the face mapping 
        procedure   :: get_face_node         !< Return an individual node index
        procedure   :: get_face_nodes        !< Return array of node indices

    end type face_connectivity_t
    !*******************************************************************************


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!  @param[in]  mapping Mapping order in the connectivity list
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self,mapping)
        class(face_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: mapping
        
        integer(ik) :: ierr, header_size, points_size, connectivity_size

        
        header_size = self%header_size
        points_size = (mapping+1)*(mapping+1)
        connectivity_size = header_size + points_size

        allocate(self%data(connectivity_size), stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !********************************************************************************








    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_face_mapping(self) result(mapping)
        class(face_connectivity_t),  intent(in)  :: self

        integer(ik) :: mapping

        mapping = self%mapping

    end function get_face_mapping
    !********************************************************************************




    !>  Set the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_face_mapping(self,mapping)
        class(face_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: mapping


        self%mapping = mapping

    end subroutine set_face_mapping
    !********************************************************************************














!    !>  Return the element mapping from a connectivity
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   5/23/2016
!    !!
!    !!
!    !!
!    !--------------------------------------------------------------------------------
!    function get_element_partition(self) result(ipartition)
!        class(element_connectivity_t),  intent(in)  :: self
!
!        integer(ik) :: ipartition
!
!        ipartition = self%partition
!
!    end function get_element_partition
!    !********************************************************************************
!
!
!
!    !>  Return the element mapping from a connectivity
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   5/23/2016
!    !!
!    !!
!    !!
!    !--------------------------------------------------------------------------------
!    subroutine set_element_partition(self,ipartition)
!        class(element_connectivity_t),  intent(inout)   :: self
!        integer(ik),                    intent(in)      :: ipartition
!
!        self%partition = ipartition
!
!    end subroutine set_element_partition
!    !********************************************************************************
!













    !>  Return an individual node index from the element connectivity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_face_node(self,idx) result(node)
        class(face_connectivity_t), intent(in)  :: self
        integer(ik),                intent(in)  :: idx

        integer(ik) :: node
        integer(ik) :: header

        header = self%header_size 

        node = self%data(header+idx)

    end function get_face_node
    !********************************************************************************






    !>  Return the array of element node indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_face_nodes(self) result(nodes)
        class(face_connectivity_t),  intent(in)  :: self

        integer(ik), allocatable    :: nodes(:)
        integer(ik)                 :: header

        header = self%header_size 

        nodes = self%data((header+1):)

    end function get_face_nodes
    !********************************************************************************


    !>  Set the array of element node indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_face_nodes(self,nodes)
        class(face_connectivity_t), intent(inout)   :: self
        integer(ik),                intent(in)      :: nodes(:)

        integer(ik)                 :: header

        header = self%header_size 

        self%data((header+1):) = nodes

    end subroutine set_face_nodes
    !********************************************************************************





end module type_face_connectivity
