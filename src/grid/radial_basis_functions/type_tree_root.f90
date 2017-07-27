module type_tree_root
    use mod_kinds, only: ik, rk
    implicit none

    type, public :: tree_root_t


        real(rk), dimension(:,:), allocatable   ::  rbf_matrix_base, inv_rbf_matrix_base

        integer(ik)                 :: surface_grid_index   !Index of the corresponding surface grid point

        integer(ik)                 :: root_ID            !Index of the root in the tree%vertices(:) array
        integer(ik), allocatable    :: containedInroot(:) !Array of the root indices whose support contains the present root
        integer(ik), allocatable    :: containsroot(:)    !Array fo the v-indices contained in the support of the present root

        !
        ! Radial Basis Function Parameters
        !
        ! This should be moved into the RBF type!
        real(rk)                    :: supportRadius
        real(rk)                    :: pointCoords(3)
        

        real(rk)                    :: rbf_values

        !
        ! Note: we could have arrays of
        !
        class(radial_basis_function_f)  :: rbf
    contains

        procedure :: init_tree_root

    end type tree_root_t


end module type_tree_root
