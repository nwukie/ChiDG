module type_element_connectivity
#include <messenger.h>
    use mod_kinds,                  only: ik
    use mod_constants,              only: NO_PARTITION
    implicit none



    !>  This data type contains the connectivity information for a single element
    !!
    !!  Connectivity data format is as follows:
    !!  data = [ idomain_g, ielement_g, mapping, ipt1, ipt2, ipt3, ... ]
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: element_connectivity_t

        integer(ik)                 :: header_size = 3  !< idomain_g, ielement_g, mapping
        integer(ik),    allocatable :: data(:)
        integer(ik)                 :: partition

    contains
        
        procedure   :: init


        procedure   :: set_domain_index         !< Set global domain index
        procedure   :: set_element_index        !< Set domain-global index
        procedure   :: set_element_mapping      !< Set the mapping
        procedure   :: set_element_partition    !< Set partition index that owns the element
        procedure   :: set_element_node        !< Set the nodes indices defining the element connectivity
        procedure   :: set_element_nodes        !< Set the nodes indices defining the element connectivity

        procedure   :: get_domain_index         !< Return the global domain index                   (index within the unpartitioned group of domains)
        procedure   :: get_element_index        !< Return the domain-global index of the element    (index within the unpartitioned domain)
        procedure   :: get_element_mapping      !< Return the element mapping 
        procedure   :: get_element_partition    !< Return the partition index that owns the element
        procedure   :: get_element_node         !< Return an individual node index
        procedure   :: get_element_nodes        !< Return array of node indices

    end type element_connectivity_t
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
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: mapping
        
        integer(ik) :: ierr, header_size, points_size, connectivity_size

        
        header_size = self%header_size
        points_size = (mapping+1)*(mapping+1)*(mapping+1)
        connectivity_size = header_size + points_size

        allocate(self%data(connectivity_size), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Initialize to zero, set mapping.
        self%data    = 0
        self%data(3) = mapping

        ! Initialize to no owner partition
        self%partition = NO_PARTITION

    end subroutine init
    !********************************************************************************








    !>  Return the global domain index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_domain_index(self) result(idomain)
        class(element_connectivity_t),  intent(in)  :: self

        integer(ik) :: idomain

        idomain = self%data(1)

    end function get_domain_index
    !********************************************************************************



    !>  Set the global domain index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_domain_index(self,idomain_g)
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik)                                     :: idomain_g

        self%data(1) = idomain_g

    end subroutine set_domain_index
    !********************************************************************************











    !>  Return the global element index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_index(self) result(ielem)
        class(element_connectivity_t),  intent(in)  :: self

        integer(ik) :: ielem

        ielem = self%data(2)

    end function get_element_index
    !********************************************************************************



    !>  Set the global element index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_element_index(self,ielement_g)
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: ielement_g

        self%data(2) = ielement_g

    end subroutine set_element_index
    !********************************************************************************
























    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_mapping(self) result(mapping)
        class(element_connectivity_t),  intent(in)  :: self

        integer(ik) :: mapping

        mapping = self%data(3)

        if ( mapping == 0 ) then
            call chidg_signal(FATAL,"element_connectivity%get_element_mapping: No mapping is set")
        end if

    end function get_element_mapping
    !********************************************************************************




    !>  Set the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_element_mapping(self,mapping)
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: mapping


        self%data(3) = mapping

    end subroutine set_element_mapping
    !********************************************************************************














    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_partition(self) result(ipartition)
        class(element_connectivity_t),  intent(in)  :: self

        integer(ik) :: ipartition

        ipartition = self%partition

    end function get_element_partition
    !********************************************************************************



    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_element_partition(self,ipartition)
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: ipartition

        self%partition = ipartition

    end subroutine set_element_partition
    !********************************************************************************














    !>  Return an individual node index from the element connectivity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_node(self,idx) result(node)
        class(element_connectivity_t),  intent(in)  :: self
        integer(ik),                    intent(in)  :: idx

        integer(ik) :: node, header, nnodes

        header = self%header_size 

        ! Check bounds
        nnodes = size(self%data) - header
        if ( idx > nnodes ) then
            call chidg_signal(FATAL,"element_connectivity%get_element_node: Node index out of bounds")
        end if

        ! Get node
        node = self%data(header+idx)

    end function get_element_node
    !********************************************************************************



    !>  Return an individual node index from the element connectivity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_element_node(self,idx,node)
        class(element_connectivity_t),  intent(inout)   :: self
        integer(ik),                    intent(in)      :: idx
        integer(ik),                    intent(in)      :: node

        integer(ik) :: header

        header = self%header_size 

        self%data(header+idx) = node

    end subroutine set_element_node
    !********************************************************************************





    !>  Return the array of element node indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_element_nodes(self) result(nodes)
        class(element_connectivity_t),  intent(in)  :: self

        integer(ik), allocatable    :: nodes(:)
        integer(ik)                 :: header

        header = self%header_size 

        nodes = self%data(header+1:)

    end function get_element_nodes
    !********************************************************************************


    !>  Set the array of element node indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine set_element_nodes(self,nodes)
        class(element_connectivity_t),  intent(inout)  :: self
        integer(ik),                    intent(in)  :: nodes(:)

        integer(ik)                 :: header

        header = self%header_size 

        self%data((header+1):) = nodes

    end subroutine set_element_nodes
    !********************************************************************************





end module type_element_connectivity
