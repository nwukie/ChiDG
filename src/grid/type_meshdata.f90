module type_meshdata
    use mod_kinds,                  only: ik
    use type_point,                 only: point_t
    use type_domain_connectivity,   only: domain_connectivity_t



    !>  Data type for returning mesh-data from a file-read routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: meshdata_t

        character(:),   allocatable :: name             !< Name of the current domain
        type(point_t),  allocatable :: points(:)        !< Array containing mesh points
        type(domain_connectivity_t) :: connectivity     !< Connectivity data for each element
        character(:),   allocatable :: eqnset           !< Equation set to allocate for the domain
        character(:),   allocatable :: coord_system     !< 'Cartesian' or 'Cylindrical'
        integer(ik)                 :: nelements_g      !< Number of elements in the unpartitioned domain
        integer(ik)                 :: spacedim         !< Number of spatial dimensions
        integer(ik)                 :: nterms_c         !< number of terms in the coordinate expansion
        integer(ik)                 :: proc             !< Processor assignment

    end type meshdata_t
    !*****************************************************************************************







end module type_meshdata
