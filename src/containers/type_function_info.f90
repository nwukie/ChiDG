module type_function_info
    use mod_kinds,  only: ik
    use type_seed,  only: seed_t


    !> Container for flux information.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  Added dtype for defining the type of differentiation required ('dQ','dX','none')
    !!  @author Matteo Ugolotti
    !!  @date   9/3/2018
    !!
    !---------------------------------------------------------------------
    type, public :: function_info_t

        integer(ik)     :: type     !< Function type. ex: 'boundary_advective_flux', 'boundary_diffusive_flux', etc.
                                    !< Also used to distinguish 'primary' and 'auxiliary' modal coefficients in update_cache 
        integer(ik)     :: ifcn     ! Index of the function of the give type being computed.
        integer(ik)     :: idepend  ! Dependency index of a given element to the function. Chimera could be > 1.
        integer(ik)     :: idiff    ! Index of the direction being linearized. XI_MIN, ETA_MAX, DIAG, etc.
        integer(ik)     :: dtype    !< Differentiation type: 'dQ_DIFF=1, 'dX_DIFF=2', 'dBC_DIFF=3', 'NO_DIFF=0'
        type(seed_t)    :: seed     ! Indices of the element being linearized

        ! BC parameters for BC differentiation
        character(:),   allocatable :: bc_param       !< BC property for AdjointBC
        logical                     :: bc_group_match !< BC group for AdjointBC

    end type function_info_t
    !*********************************************************************

end module type_function_info
