module type_mesh_motion_wrapper
    use mod_kinds, only: ik
    use type_mesh_motion,          only: mesh_motion_t
    implicit none

    type    :: mesh_motion_wrapper_t

        class(mesh_motion_t),  allocatable :: mm

    end type mesh_motion_wrapper_t



end module type_mesh_motion_wrapper
