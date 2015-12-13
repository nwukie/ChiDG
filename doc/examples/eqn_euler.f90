

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


        self%name    = 'Euler'

        !
        ! Allocate equation set properties.
        !
        if (allocated(self%prop)) deallocate(self%prop)
        allocate(EULER_properties_t::self%prop)

        select type (prop => self%prop)
            type is (EULER_properties_t) 
                allocate(perfect_gas_t::prop%fluid)
        end select



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
