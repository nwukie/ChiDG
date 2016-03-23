module eqn_primitive_linearized_euler
    use mod_kinds,                                      only: rk,ik

    use type_equationset,                               only: equationset_t
    use perfect_gas,                                    only: perfect_gas_t

    use PRIMLINEULER_boundary_average_advective_flux_real,  only: PRIMLINEULER_boundary_average_advective_flux_real_t
    use PRIMLINEULER_boundary_average_advective_flux_imag,  only: PRIMLINEULER_boundary_average_advective_flux_imag_t
    use PRIMLINEULER_volume_advective_flux_real,            only: PRIMLINEULER_volume_advective_flux_real_t
    use PRIMLINEULER_volume_advective_flux_imag,            only: PRIMLINEULER_volume_advective_flux_imag_t
    use PRIMLINEULER_LaxFriedrichs_flux_real,               only: PRIMLINEULER_LaxFriedrichs_flux_real_t
    use PRIMLINEULER_LaxFriedrichs_flux_imag,               only: PRIMLINEULER_LaxFriedrichs_flux_imag_t
    use PRIMLINEULER_volume_advective_source_real,          only: PRIMLINEULER_volume_advective_source_real_t
    use PRIMLINEULER_volume_advective_source_imag,          only: PRIMLINEULER_volume_advective_source_imag_t
    use PRIMLINEULER_properties,                            only: PRIMLINEULER_properties_t
    implicit none

    private






    !>  The Euler Equations, governing inviscid fluid flows.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !------------------------------------------------------------------------------------------------------
    type, extends(equationset_t), public :: primitive_linearized_euler_e


    contains

        procedure   :: init

    end type primitive_linearized_euler_e
    !******************************************************************************************************








contains






    !>  Initialize the Euler equations
    !!
    !!   - initialize properties
    !!   - add equations
    !!   - assign flux components
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(primitive_linearized_euler_e), intent(inout) :: self


        type(PRIMLINEULER_volume_advective_flux_real_t)             :: volume_flux_real
        type(PRIMLINEULER_volume_advective_flux_imag_t)             :: volume_flux_imag
        type(PRIMLINEULER_boundary_average_advective_flux_real_t)   :: average_flux_real
        type(PRIMLINEULER_boundary_average_advective_flux_imag_t)   :: average_flux_imag
        type(PRIMLINEULER_LaxFriedrichs_flux_real_t)                :: lax_real
        type(PRIMLINEULER_LaxFriedrichs_flux_imag_t)                :: lax_imag
        type(PRIMLINEULER_volume_advective_source_real_t)           :: volume_source_real
        type(PRIMLINEULER_volume_advective_source_imag_t)           :: volume_source_imag
        type(PRIMLINEULER_properties_t)                             :: prop

        type(perfect_gas_t)                                     :: perfect_gas

        !
        ! Set equationset name.
        !
        call self%set_name("PrimitiveLinearizedEuler")


        !
        ! Set up properties
        !
        call prop%add_fluid(perfect_gas)


        !
        ! Allocate equation set properties
        !
        call self%add_properties(prop)


        !
        ! Allocate and initialize equations
        !
        call self%add_equation("rho_r", 1)
        call self%add_equation("u_r",   2)
        call self%add_equation("v_r",   3)
        call self%add_equation("w_r",   4)
        call self%add_equation("p_r",   5)
        call self%add_equation("rho_i", 6)
        call self%add_equation("u_i",   7)
        call self%add_equation("v_i",   8)
        call self%add_equation("w_i",   9)
        call self%add_equation("p_i",   10)



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


    end subroutine init
    !*****************************************************************************************************







end module eqn_primitive_linearized_euler
