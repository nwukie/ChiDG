module eqn_graddemo
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
    !-------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: graddemo

    contains

        procedure   :: init
        procedure   :: build

    end type graddemo
    !*************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo),   intent(inout)  :: self

        call self%set_name('graddemo')

    end subroutine init
    !*************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    function build(self,blueprint) result(graddemo_eqns)
        class(graddemo),    intent(in)  :: self
        character(*),       intent(in)  :: blueprint

        type(equation_set_t)            :: graddemo_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call graddemo_eqns%set_name('graddemo')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call graddemo_eqns%add_operator('Graddemo Volume Flux')

                call graddemo_eqns%add_model('Ideal Gas')
                call graddemo_eqns%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_graddemo: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !***********************************************************************************



end module eqn_graddemo
