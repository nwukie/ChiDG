module type_prescribed_mesh_motion_group_wrapper
    use mod_kinds, only: ik
    use type_prescribed_mesh_motion_group,          only: prescribed_mesh_motion_group_t
    implicit none

    type    :: prescribed_mesh_motion_group_wrapper_t

        integer(ik) :: ngroups = 0 
        type(prescribed_mesh_motion_group_t),  allocatable :: pmm_groups(:)

    end type prescribed_mesh_motion_group_wrapper_t



end module type_prescribed_mesh_motion_group_wrapper
