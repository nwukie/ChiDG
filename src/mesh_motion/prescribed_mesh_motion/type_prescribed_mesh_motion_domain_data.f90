module type_prescribed_mesh_motion_domain_data
    
    type :: prescribed_mesh_motion_domain_data_t

        character(:), allocatable        :: domain_name
        character(:), allocatable        :: pmm_group_name


    end type prescribed_mesh_motion_domain_data_t

end module type_prescribed_mesh_motion_domain_data
