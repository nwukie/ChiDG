module type_rbf_node_patch_mm
    use type_rbf_node_patch
    use mod_kinds,                      only: ik, rk
    implicit none

    type, extends(rbf_node_patch_t) :: rbf_node_patch_mm_t


        real(rk), allocatable   :: dnodes(:,:)
        real(rk), allocatable   :: vnodes(:,:)

    contains

        procedure :: init_mm

    end type rbf_node_patch_mm_t


contains


end module type_rbf_node_patch_mm
