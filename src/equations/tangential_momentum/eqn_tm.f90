module eqn_tm
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
    !------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: tm

    contains

        procedure   :: init
        procedure   :: build

    end type tm
    !******************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(tm),   intent(inout)  :: self

        call self%set_name('TM')

    end subroutine init
    !******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function build(self,blueprint) result(tm_eqns)
        class(tm),      intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: tm_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call tm_eqns%set_name('TM')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call tm_eqns%add_operator('TM Volume Flux')
                call tm_eqns%add_operator('TM Boundary Average Flux')
                call tm_eqns%add_operator('TM BC Flux')
                call tm_eqns%add_operator('TM LaxFriedrichs Flux')
                call tm_eqns%add_operator('TM Volume Cylindrical Source')

                call tm_eqns%add_model('RAC')

            case default
                call chidg_signal_one(FATAL, "build_tm: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !*******************************************************************************************



end module eqn_tm
