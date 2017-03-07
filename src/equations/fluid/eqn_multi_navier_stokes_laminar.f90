module eqn_multi_navier_stokes_laminar
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    use perfect_gas,                only: perfect_gas_t
    implicit none


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: multi_navier_stokes_laminar

    contains

        procedure   :: init
        procedure   :: build

    end type multi_navier_stokes_laminar
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(multi_navier_stokes_laminar),   intent(inout)  :: self

        call self%set_name('Multi Navier Stokes Laminar')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(equation_set)
        class(multi_navier_stokes_laminar),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        type(equation_set_t)            :: equation_set
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call equation_set%set_name(self%get_name())
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call equation_set%add_operator('Euler Volume Flux')
                call equation_set%add_operator('Euler Boundary Average Flux')
                call equation_set%add_operator('Euler Roe Flux')
                call equation_set%add_operator('Euler BC Flux')
                call equation_set%add_operator('Euler Volume Cylindrical Source')

                call equation_set%add_operator('Fluid Viscous Volume Operator')
                call equation_set%add_operator('Fluid Viscous Boundary Average Operator')
                call equation_set%add_operator('Fluid Viscous BC Operator')
                call equation_set%add_operator('Fluid Viscous Volume Cylindrical Source')

                call equation_set%add_model('Ideal Gas')
                call equation_set%add_model('Constant Viscosity Laminar')
                call equation_set%add_model('Stokes Hypothesis')
                call equation_set%add_model('Reynolds Analogy')
                call equation_set%add_model('Zero Turbulent Model Fields')
                call equation_set%add_model('Shear Stress')
                call equation_set%add_model('Temperature Gradient')


                call equation_set%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_laminar_navier_stokes: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************






end module eqn_multi_navier_stokes_laminar
