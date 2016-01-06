

    !> The Euler Equations, governing inviscid fluid flows.
    !!
    !!   @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    !> [euler_e]
    type, extends(equationset_t), public :: euler_e


    contains
        procedure   :: init

    end type euler_e
    !> [euler_e]
    !----------------------------------------------------------------------------



contains


    !> Initialize the Euler equations
    !!
    !!   - initialize properties
    !!   - add equations
    !!   - add flux components
    !!
    !----------------------------------------------------------------------------
    !> [euler_e init]
    subroutine init(self)
        class(euler_e), intent(inout) :: self


        type(EULER_volume_advective_flux_t)             :: volume_flux
        type(EULER_boundary_average_advective_flux_t)   :: average_flux
        type(EULER_Roe_flux_t)                          :: roe
        type(EULER_properties_t)                        :: prop

        type(perfect_gas_t)                             :: perfect_gas


        !
        ! Set equation set name
        !
        self%name    = 'Euler'


        !
        ! Equation set properties
        !
        call prop%add_fluid(perfect_gas)    ! set up properties
        call self%add_properties(prop)      ! add to equation set


        !
        ! Add equations
        !
        call self%add_equation("rho",1)
        call self%add_equation("rhou",2)
        call self%add_equation("rhov",3)
        call self%add_equation("rhow",4)
        call self%add_equation("rhoE",5)



        !
        ! Add flux components
        !
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(roe)
        call self%add_volume_advective_flux(volume_flux)


    end subroutine
    !> [euler_e init]
    !----------------------------------------------------------------------------







end module eqn_euler
