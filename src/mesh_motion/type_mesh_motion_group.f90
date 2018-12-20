module type_mesh_motion_group
    use type_mesh_motion,        only: mesh_motion_t
    implicit none
    
    type :: mesh_motion_group_t

        character(:), allocatable            :: name
        class(mesh_motion_t), allocatable     :: mm

    end type mesh_motion_group_t
    

end module type_mesh_motion_group
