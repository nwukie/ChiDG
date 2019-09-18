module type_rbf_node
    use mod_kinds,                          only: ik, rk
    implicit none

    type, public :: rbf_node_t

        real(rk)                        :: node_center(3)
        real(rk)                        :: node_radius(3)
    
    contains

        procedure :: set_node

    end type rbf_node_t

contains

    subroutine set_node(self, node_center_in, node_radius_in) 
        class(rbf_node_t), intent(inout) :: self
        real(rk), intent(in) :: node_center_in(3)
        real(rk), intent(in) :: node_radius_in(3)


        self%node_center = node_center_in
        self%node_radius = node_radius_in

    end subroutine set_node

end module type_rbf_node
