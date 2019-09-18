module type_rbf_node_patch_mm_pmm
    use type_rbf_node_patch_mm
    use mod_kinds,                      only: ik, rk
    implicit none

    type, extends(rbf_node_patch_mm_t) :: rbf_node_patch_mm_pmm_t


        real(rk), allocatable   :: pmmf_ID(:)

    contains

        procedure :: init_mm

    end type rbf_node_patch_mm_pmm_t


contains


end module type_rbf_node_patch_mm
