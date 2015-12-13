module eqn_euler
    use mod_kinds,                              only: rk,ik

    use type_equationset,                       only: equationset_t
    use perfect_gas,                            only: perfect_gas_t

    use EULER_boundary_average_advective_flux,  only: EULER_boundary_average_advective_flux_t
    use EULER_volume_advective_flux,            only: EULER_volume_advective_flux_t
    use EULER_LaxFriedrichs_flux,               only: EULER_LaxFriedrichs_flux_t
    use EULER_Roe_flux,                         only: EULER_Roe_flux_t
    use EULER_properties,                       only: EULER_properties_t
    implicit none
    private




    !> The Euler Equations, governing inviscid fluid flows.
    !!
    !!   @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    type, extends(equationset_t), public :: euler_e


    contains
        procedure   :: init

    end type euler_e
    !----------------------------------------------------------------------------









contains


    !> Initialize the Euler equations
    !!
    !!   - initialize properties
    !!   - add equations
    !!   - assign flux components
    !!
    !!   @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    subroutine init(self)
        class(euler_e), intent(inout) :: self


        type(EULER_volume_advective_flux_t)             :: volume_flux
        type(EULER_boundary_average_advective_flux_t)   :: average_flux
        type(EULER_LaxFriedrichs_flux_t)                :: lax
        type(EULER_Roe_flux_t)                          :: roe


        self%name    = 'Euler'

        !
        ! Allocate equation set properties. TODO: Delegate this to type-bound procedure of equationset_t. 'add_properties(prop)'
        !
        if (allocated(self%prop)) deallocate(self%prop)
        allocate(EULER_properties_t::self%prop)

        select type (prop => self%prop)
            type is (EULER_properties_t) 
                allocate(perfect_gas_t::prop%fluid)
        end select



        !
        ! Allocate and initialize equations
        !
        call self%add_equation("rho",1)
        call self%add_equation("rhou",2)
        call self%add_equation("rhov",3)
        call self%add_equation("rhow",4)
        call self%add_equation("rhoE",5)



        !
        ! Allocate flux components to specific types for the equation set
        !
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(roe)
        !call self%add_boundary_advective_flux(lax)
        call self%add_volume_advective_flux(volume_flux)


    end subroutine
    !----------------------------------------------------------------------------







end module eqn_euler
