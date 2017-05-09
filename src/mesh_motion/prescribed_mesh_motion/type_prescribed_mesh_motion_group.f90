module type_prescribed_mesh_motion_group
    use type_prescribed_mesh_motion,        only: prescribed_mesh_motion_t
    implicit none
    
    type :: prescribed_mesh_motion_group_t

        character(:), allocatable            :: name
        class(prescribed_mesh_motion_t), allocatable     :: pmm

    end type prescribed_mesh_motion_group_t
    

end module type_prescribed_mesh_motion_group
