module type_mesh_motion_group_wrapper
    use mod_kinds, only: ik
    use type_mesh_motion_group,          only: mesh_motion_group_t
    implicit none

    type    :: mesh_motion_group_wrapper_t

        integer(ik) :: ngroups = 0 
        class(mesh_motion_group_t),  allocatable :: mm_groups(:)

    end type mesh_motion_group_wrapper_t



end module type_mesh_motion_group_wrapper
