module eqn_laminar_navier_stokes_ale
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    implicit none


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: laminar_navier_stokes_ale

    contains

        procedure   :: init
        procedure   :: build

    end type laminar_navier_stokes_ale
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(laminar_navier_stokes_ale),   intent(inout)  :: self

        call self%set_name('Laminar Navier Stokes ALE')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(laminar_navier_stokes_ale_eqns)
        class(laminar_navier_stokes_ale),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        type(equation_set_t)            :: laminar_navier_stokes_ale_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call laminar_navier_stokes_ale_eqns%set_name(self%get_name())
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call laminar_navier_stokes_ale_eqns%add_operator('Euler ALE Volume Flux')
                call laminar_navier_stokes_ale_eqns%add_operator('Euler ALE Boundary Average Flux')
                call laminar_navier_stokes_ale_eqns%add_operator('Euler ALE Roe Flux')
                call laminar_navier_stokes_ale_eqns%add_operator('Euler ALE BC Flux')

                call laminar_navier_stokes_ale_eqns%add_operator('Fluid Viscous ALE Volume Operator')
                call laminar_navier_stokes_ale_eqns%add_operator('Fluid Viscous ALE Boundary Average Operator')
                call laminar_navier_stokes_ale_eqns%add_operator('Fluid Viscous ALE BC Operator')

                call laminar_navier_stokes_ale_eqns%add_model('Ideal Gas')
                !call laminar_navier_stokes_ale_eqns%add_model('Fluid Advection Velocity')
                !call laminar_navier_stokes_ale_eqns%add_model('Sutherlands Law')
                call laminar_navier_stokes_ale_eqns%add_model('Constant Viscosity')
                call laminar_navier_stokes_ale_eqns%add_model('Stokes Hypothesis')
                call laminar_navier_stokes_ale_eqns%add_model('Reynolds Analogy')
                call laminar_navier_stokes_ale_eqns%add_model('Zero Turbulent Model Fields')
                call laminar_navier_stokes_ale_eqns%add_model('Shear Stress')
                call laminar_navier_stokes_ale_eqns%add_model('Temperature Gradient')


                call laminar_navier_stokes_ale_eqns%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_laminar_navier_stokes_ale: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************






end module eqn_laminar_navier_stokes_ale
