module eqn_euler_ale
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_solverdata,            only: solverdata_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    implicit none


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: euler_ale

    contains

        procedure   :: init
        procedure   :: build

    end type euler_ale
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale),   intent(inout)  :: self

        call self%set_name('Euler ALE')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(euler_ale_eqns)
        class(euler_ale),   intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: euler_ale_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call euler_ale_eqns%set_name('Euler ALE')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call euler_ale_eqns%add_operator('Euler ALE Volume Flux')
                call euler_ale_eqns%add_operator('Euler ALE Boundary Average Flux')
                call euler_ale_eqns%add_operator('Euler ALE Roe Flux')
                call euler_ale_eqns%add_operator('Euler ALE BC Flux')


                call euler_ale_eqns%add_model('Ideal Gas')
            !    call euler_ale_eqns%add_model('Fluid Advection Velocity')

                call euler_ale_eqns%add_pseudo_timestep(fluid_pseudo_time)


                call euler_ale_eqns%add_operator('Geometric Conservation Volume Operator')
                call euler_ale_eqns%add_operator('Geometric Conservation Boundary Average Operator')
!                call euler_ale_eqns%add_operator('Geometric Conservation LaxFriedrichs Operator')
                call euler_ale_eqns%add_operator('Geometric Conservation BC Operator')

            case default
                call chidg_signal_one(FATAL, "build_euler_ale: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************

















end module eqn_euler_ale
