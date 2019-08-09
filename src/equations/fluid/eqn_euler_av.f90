module eqn_euler_av
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
    type, public, extends(equation_builder_t) :: euler_av

    contains

        procedure   :: init
        procedure   :: build

    end type euler_av
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_av),   intent(inout)  :: self

        call self%set_name('Euler AV')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(euler_av_eqns)
        class(euler_av),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        type(equation_set_t)            :: euler_av_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call euler_av_eqns%set_name('Euler AV')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call euler_av_eqns%add_operator('Euler Volume Flux')
                call euler_av_eqns%add_operator('Euler Boundary Average Flux')
                call euler_av_eqns%add_operator('Euler Roe Flux')
                call euler_av_eqns%add_operator('Euler BC Flux')
                call euler_av_eqns%add_operator('Euler Volume Cylindrical Source')

                call euler_av_eqns%add_operator('Fluid Laplacian AV Volume Operator')
                call euler_av_eqns%add_operator('Fluid Laplacian AV Boundary Average Operator')
                call euler_av_eqns%add_operator('Fluid Laplacian AV BC Operator')

                call euler_av_eqns%add_model('Regularized Fluid Primary Fields')
                call euler_av_eqns%add_model('Regularized Ideal Gas')
                call euler_av_eqns%add_model('Pressure Gradient')
                call euler_av_eqns%add_model('Velocity Gradient')
                call euler_av_eqns%add_model('Velocity Divergence and Curl')
                call euler_av_eqns%add_model('Critical Sound Speed')
                call euler_av_eqns%add_model('MNPH Shock Sensor')
                call euler_av_eqns%add_model('MNPH Artificial Viscosity')

                call euler_av_eqns%add_pseudo_timestep(fluid_pseudo_time)

            case default
                call chidg_signal_one(FATAL, "build_euler_av: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************






end module eqn_euler_av
