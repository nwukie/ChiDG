module eqn_linearized_euler
    use mod_kinds,                              only: rk,ik

    use type_equationset,                       only: equationset_t
    use perfect_gas,                            only: perfect_gas_t

    use LINEULER_boundary_average_advective_flux_real,  only: LINEULER_boundary_average_advective_flux_real_t
    use LINEULER_boundary_average_advective_flux_imag,  only: LINEULER_boundary_average_advective_flux_imag_t
    use LINEULER_volume_advective_flux_real,            only: LINEULER_volume_advective_flux_real_t
    use LINEULER_volume_advective_flux_imag,            only: LINEULER_volume_advective_flux_imag_t
    use LINEULER_LaxFriedrichs_flux_real,               only: LINEULER_LaxFriedrichs_flux_real_t
    use LINEULER_LaxFriedrichs_flux_imag,               only: LINEULER_LaxFriedrichs_flux_imag_t
    use LINEULER_volume_advective_source_real,          only: LINEULER_volume_advective_source_real_t
    use LINEULER_volume_advective_source_imag,          only: LINEULER_volume_advective_source_imag_t
    use LINEULER_volume_advective_sourceterms_real,     only: LINEULER_volume_advective_sourceterms_real_t

    use LINEULER_properties,                            only: LINEULER_properties_t
    implicit none
    private




    !> The Euler Equations, governing inviscid fluid flows.
    !!
    !!   @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    type, extends(equationset_t), public :: linearized_euler_e


    contains
        procedure   :: init

    end type linearized_euler_e
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
        class(linearized_euler_e), intent(inout) :: self


        type(LINEULER_volume_advective_flux_real_t)             :: volume_flux_real
        type(LINEULER_volume_advective_flux_imag_t)             :: volume_flux_imag
        type(LINEULER_boundary_average_advective_flux_real_t)   :: average_flux_real
        type(LINEULER_boundary_average_advective_flux_imag_t)   :: average_flux_imag
        type(LINEULER_LaxFriedrichs_flux_real_t)                :: lax_real
        type(LINEULER_LaxFriedrichs_flux_imag_t)                :: lax_imag
        type(LINEULER_volume_advective_source_real_t)           :: volume_source_real
        type(LINEULER_volume_advective_source_imag_t)           :: volume_source_imag
        type(LINEULER_volume_advective_sourceterms_real_t)      :: volume_sourceterms_real


        self%name    = 'Linearized Euler'

        !
        ! Allocate equation set properties
        !
        if (allocated(self%prop)) deallocate(self%prop)
        allocate(LINEULER_properties_t::self%prop)

        select type (prop => self%prop)
            type is (LINEULER_properties_t) 
                allocate(perfect_gas_t::prop%fluid)
        end select



        !
        ! Allocate and initialize equations
        !
        call self%add_equation("rho_r",1)
        call self%add_equation("rhou_r",2)
        call self%add_equation("rhov_r",3)
        call self%add_equation("rhow_r",4)
        call self%add_equation("rhoE_r",5)

        call self%add_equation("rho_i",6)
        call self%add_equation("rhou_i",7)
        call self%add_equation("rhov_i",8)
        call self%add_equation("rhow_i",9)
        call self%add_equation("rhoE_i",10)


        !
        ! Allocate flux components to specific types for the equation set
        !
        call self%add_boundary_advective_flux(average_flux_real)
        call self%add_boundary_advective_flux(average_flux_imag)
        call self%add_boundary_advective_flux(lax_real)
        call self%add_boundary_advective_flux(lax_imag)
        call self%add_volume_advective_flux(volume_flux_real)
        call self%add_volume_advective_flux(volume_flux_imag)
        call self%add_volume_advective_flux(volume_source_real)
        call self%add_volume_advective_flux(volume_source_imag)
        call self%add_volume_advective_flux(volume_sourceterms_real)


    end subroutine
    !----------------------------------------------------------------------------







end module eqn_linearized_euler
