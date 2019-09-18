module type_mesh_motion_domain_data
    
    type :: mesh_motion_domain_data_t

        character(:), allocatable        :: domain_name
        character(:), allocatable        :: mm_group_name


    end type mesh_motion_domain_data_t

end module type_mesh_motion_domain_data
