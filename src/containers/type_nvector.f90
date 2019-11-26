module type_nvector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_node,          only: node_t
    implicit none




    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, public :: nvector_t

        ! List attributes
        integer(ik)             :: size_        = 0
        integer(ik)             :: capacity_    = 0
        integer(ik)             :: buffer_      = 50

        ! Data storage
        type(node_t), allocatable :: data(:)


    contains

        procedure, public   :: size
        procedure, public   :: capacity


        !< Data modifiers
        procedure, public   :: push_back
        procedure, public   :: clear
        procedure, private  :: increase_capacity


        !< Data accessors
        procedure, public   :: at
        procedure, public   :: search_by_coords
        procedure, public   :: reorder_by_index
        procedure, public   :: get_nodes_coords
        procedure, public   :: get_nodes_sensitivities
        procedure, public   :: get_nodes_domain_g
        procedure, public   :: get_nodes_ID

    end type nvector_t
    !********************************************************************************



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !--------------------------------------------------------------------------------------
    function size(self) result(res)
        class(nvector_t),   intent(inout)   :: self

        integer(ik) :: res

        res = self%size_

    end function size
    !**************************************************************************************



    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !--------------------------------------------------------------------------------------
    function capacity(self) result(res)
        class(nvector_t),   intent(inout)   :: self

        integer(ik) :: res

        res = self%capacity_

    end function capacity
    !**************************************************************************************









    !> Store element at end of vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine push_back(self,new_node)
        class(nvector_t),   intent(inout)   :: self
        type(node_t),       intent(in)      :: new_node

        logical     :: capacity_reached
        integer(ik) :: size


        !
        ! Test if container has storage available. If not, then increase capacity
        !
        capacity_reached = (self%size() == self%capacity())
        if (capacity_reached) then
            call self%increase_capacity()
        end if


        !
        ! Add element to end of vector
        !
        size = self%size()
        self%data(size + 1) = new_node


        !
        ! Increment number of stored elements
        !
        self%size_ = self%size_ + 1


    end subroutine push_back
    !****************************************************************************************






    !> Clear container contents
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine clear(self)
        class(nvector_t),   intent(inout)   :: self

        self%size_      = 0
        self%capacity_  = 0

        if (allocated(self%data)) deallocate(self%data)

    end subroutine clear
    !****************************************************************************************








    !> Access element at index location
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(nvector_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: index

        type(node_t)    :: res
        logical         :: out_of_bounds

        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,'nvector_t%at: out of bounds access')
        end if


        !
        ! Allocate result
        !
        res = self%data(index)

    end function at
    !****************************************************************************************









    !> Increase the storage capacity of the vector by a buffer size predefined in the container
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(nvector_t),   intent(inout)   :: self

        type(node_t), allocatable   :: temp(:)
        integer(ik)                 :: newsize, ierr


        !
        ! Allocate temporary vector of current size plus a buffer
        !
        if ( allocated(self%data) ) then
            newsize = ubound(self%data,1) + self%buffer_
        else
            newsize = self%buffer_
        end if

        allocate(temp(newsize),stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any current data to temporary vector
        !
        if (allocated(self%data)) then
            temp(lbound(self%data,1):ubound(self%data,1))  =  self%data
        end if


        !
        ! Move alloc to move data back to self%data and deallocate temp
        !
        call move_alloc(FROM=temp,TO=self%data)


        !
        ! Reset capacity info
        !
        self%capacity_ = newsize


    end subroutine increase_capacity
    !*****************************************************************************************











    !>  Search node by coordinates, return the index of the node having the given coordiantes
    !!  and returns the node index in the list if found. If not, returns 0. 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function search_by_coords(self,coords,inode_ID) result(index_)
        class(nvector_t),      intent(inout)   :: self
        real(rk),              intent(in)      :: coords(3)
        integer(ik),           intent(in)      :: inode_ID

        integer(ik)     :: index_
        integer(ik)     :: inode
        logical         :: coords_match, id_match

        !
        ! Default, node not found
        !
        index_ = 0
        coords_match = .false.

        do inode = 1,self%size()

            coords_match = self%data(inode)%query_coords(coords(1),coords(2),coords(3))
            id_match     = ( self%data(inode)%node_ID_l == inode_ID )
            if (coords_match .and. id_match) then
                index_ = inode
                exit
            end if

        end do

        !
        ! Sanity check, if the coordinates match, also the local ID should match.
        !
        !if (present(inode_ID) .and. (coords_match)) then
        !    mismatch = (self%data(inode)%node_ID_l /= inode_ID)
        !    if (mismatch) call chidg_signal(FATAL,'nvector_t%search_by_coords: the given node-coordinates match the coordinates of a nodes already processed by their IDs do not match. Check implementation.')
        !end if

    end function search_by_coords
    !****************************************************************************************













    !>  Search node by coordinates, return the index of the node having the given coordiantes
    !!  and returns the node index in the list if found. If not, returns 0. 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine reorder_by_index(self)
        class(nvector_t),       intent(inout)  :: self

        integer(ik)                     :: inode, current_index, ierr
        type(node_t),   allocatable     :: temp_nodes(:)


        !
        ! Allocate vectors 
        !
        allocate(temp_nodes(self%size()),stat=ierr)
        if (ierr/=0) call AllocationError
        
        !
        ! Create a list of all the nodes ordered by node_ID
        !
        do inode = 1,self%size()
            temp_nodes(self%data(inode)%node_ID_l) = self%at(inode)
        end do

        
        !
        ! Move alloc to move data back to self%data and deallocate temp
        !
        call move_alloc(FROM=temp_nodes,TO=self%data)
        
        
    end subroutine reorder_by_index
    !****************************************************************************************









    !>  Return a array with each line containing the coordinates in i, j and k directions 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nodes_coords(self) result(coords)
        class(nvector_t),       intent(inout)  :: self

        real(rk),   allocatable :: coords(:,:)
        integer(ik)             :: inode, ierr


        !
        ! Allocate output array
        !
        allocate(coords(self%size(),3),stat=ierr)
        if (ierr/=0) call AllocationError

        do inode = 1,self%size()

            coords(inode,:) = self%data(inode)%get_coords()

        end do

    end function get_nodes_coords
    !****************************************************************************************









    !>  Return a array with each line containing the sensitivities in i, j and k directions 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nodes_sensitivities(self) result(sens)
        class(nvector_t),       intent(inout)  :: self

        real(rk),   allocatable :: sens(:,:)
        integer(ik)             :: inode, ierr


        !
        ! Allocate output array
        !
        allocate(sens(self%size(),3),stat=ierr)
        if (ierr/=0) call AllocationError

        do inode = 1,self%size()

            sens(inode,:) = self%data(inode)%get_sensitivities()

        end do

    end function get_nodes_sensitivities
    !****************************************************************************************










    !>  Return a array with each line containing the global domain index that the ith node 
    !!  belongs to.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nodes_domain_g(self) result(dom_indeces)
        class(nvector_t),       intent(inout)  :: self

        integer(ik),    allocatable :: dom_indeces(:)
        integer(ik)                 :: inode, ierr


        !
        ! Allocate output array
        !
        allocate(dom_indeces(self%size()),stat=ierr)
        if (ierr/=0) call AllocationError

        do inode = 1,self%size()

            dom_indeces(inode) = self%data(inode)%domain_g

        end do

    end function get_nodes_domain_g
    !****************************************************************************************










    !>  Return a array with each line containing node index of the ith node 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nodes_ID(self) result(node_indeces)
        class(nvector_t),       intent(inout)  :: self

        integer(ik),    allocatable :: node_indeces(:)
        integer(ik)                 :: inode, ierr


        !
        ! Allocate output array
        !
        allocate(node_indeces(self%size()),stat=ierr)
        if (ierr/=0) call AllocationError

        do inode = 1,self%size()

            node_indeces(inode) = self%data(inode)%node_ID_l

        end do

    end function get_nodes_ID
    !****************************************************************************************





end module type_nvector
