module type_node
#include <messenger.h>
    use mod_kinds,      only: ik,rk
    implicit none


    !>  Node data type, containing: 
    !!      - node_ID_l
    !!      - node_ID_g
    !!      - three spatial coordinates
    !!      - list of element indeces that contributed to sensitivitites
    !!      - node sensitivitites
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2018
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: node_t

        integer(ik)             :: node_ID_l         ! node ID local (from connectivity)
        integer(ik)             :: node_ID_g         ! node ID global !NOT USED
        integer(ik)             :: domain_g          ! Global domain index of the domain
        real(rk)                :: coords(3)         ! x, y, z
        real(rk)                :: sensitivities(3)  ! Grid-node sensitivities
        integer(ik)             :: coordinate_system ! 1. CARTESIAN, 2. CYLINDRICAL

    contains

        procedure   :: init_node
        procedure   :: x
        procedure   :: y
        procedure   :: z
        procedure   :: query_coords
        procedure   :: add_sensitivities
        procedure   :: get_coords
        procedure   :: get_sensitivities

    end type node_t
    !******************************************************************************************

contains

    
    
    
    !>  initialize node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_node(self,ID_l,idomain_g,coords_system,coords,sensitivities)
        class(node_t),      intent(inout)   :: self
        integer(ik),        intent(in)      :: ID_l
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: coords_system
        real(rk),           intent(in)      :: coords(3)
        real(rk),           intent(in)      :: sensitivities(3)


        self%node_ID_l          = ID_l
        self%domain_g           = idomain_g
        self%coordinate_system  = coords_system
        self%coords             = coords
        self%sensitivities      = sensitivities 




    end subroutine init_node
    !*****************************************************************************************








    !>  get x coordiante of the node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function x(self) result(x_coord)
        class(node_t),      intent(inout)   :: self

        real(rk)        :: x_coord

        x_coord = self%coords(1)


    end function x
    !*****************************************************************************************







    !>  get y coordiante of the node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function y(self) result(y_coord)
        class(node_t),      intent(inout)   :: self

        real(rk)        :: y_coord

        y_coord = self%coords(2)


    end function y
    !*****************************************************************************************







    !>  get z coordiante of the node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function z(self) result(z_coord)
        class(node_t),      intent(inout)   :: self

        real(rk)        :: z_coord

        z_coord = self%coords(3)


    end function z
    !*****************************************************************************************







    !>  Check if a set of coordinates corresponed to those of the current node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function query_coords(self,x,y,z) result(match)
        class(node_t),      intent(inout)   :: self
        real(rk),           intent(in)      :: x
        real(rk),           intent(in)      :: y
        real(rk),           intent(in)      :: z

        logical     :: match, x_match, y_match, z_match

        match = .false.

        x_match = (self%coords(1) == x)
        y_match = (self%coords(2) == y)
        z_match = (self%coords(3) == z)

        if (x_match .and. y_match .and. z_match) then
            match = .true.
        end if

    end function query_coords
    !*****************************************************************************************









    !>  Add contribution to node sensitivities
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine add_sensitivities(self,node_sens) 
        class(node_t),      intent(inout)   :: self
        real(rk),           intent(in)      :: node_sens(3)


        self%sensitivities(1) = self%sensitivities(1) + node_sens(1)
        self%sensitivities(2) = self%sensitivities(2) + node_sens(2)
        self%sensitivities(3) = self%sensitivities(3) + node_sens(3)


    end subroutine add_sensitivities
    !*****************************************************************************************







    !>  Get coordinates of the node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_coords(self) result(coords)
        class(node_t),      intent(inout)   :: self

        real(rk)        :: coords(3)

        
        coords(1) = self%coords(1)
        coords(2) = self%coords(2)
        coords(3) = self%coords(3)
        

    end function get_coords
    !*****************************************************************************************








    !>  Get sensitivities of the node 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_sensitivities(self) result(sens)
        class(node_t),      intent(inout)   :: self

        real(rk)        :: sens(3)

        
        sens(1) = self%sensitivities(1)
        sens(2) = self%sensitivities(2)
        sens(3) = self%sensitivities(3)
        

    end function get_sensitivities
    !*****************************************************************************************


end module type_node
