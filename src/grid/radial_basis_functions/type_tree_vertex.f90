module type_tree_vertex
    use mod_kinds, only: ik, rk
    implicit none

    type, public :: tree_vertex_t

        integer(ik)                 :: surface_grid_index   !Index of the corresponding surface grid point

        integer(ik)                 :: tree_level

        integer(ik)                 :: vertex_ID            !Index of the vertex in the tree%vertices(:) array
        integer(ik), allocatable    :: containedInVertex(:) !Array of the vertex indices whose support contains the present vertex
        integer(ik), allocatable    :: containsVertex(:)    !Array fo the v-indices contained in the support of the present vertex

        !
        ! Radial Basis Function Parameters
        !
        ! This should be moved into the RBF type!
        real(rk)                    :: support_radius
        real(rk)                    :: support_node(3)
        

        real(rk)                    :: rbf_values

        !
        ! Note: we could have arrays of
        !
        class(radial_basis_function_f)  :: rbf
    contains

        procedure :: init_tree_vertex
        procedure :: get_rbf_term

    end type tree_vertex_t

contains

    function get_rbf_term(self,node) result(val)
        class(tree_vertex_t)    :: self
        real(rk)                :: node(3)

        real(rk) :: val

        val = self%rbf%compute(node)

    end function getRBFContribution
end module type_tree_vertex
