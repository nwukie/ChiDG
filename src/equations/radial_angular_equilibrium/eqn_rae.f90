module eqn_rae
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
    type, public, extends(equation_builder_t) :: rae

    contains

        procedure   :: init
        procedure   :: build

    end type rae
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(rae),   intent(inout)  :: self

        call self%set_name('RAE')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(rae_eqns)
        class(rae),     intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rae_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call rae_eqns%set_name('RAE')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call rae_eqns%add_operator('RAE Volume Flux')
                call rae_eqns%add_operator('RAE Boundary Average Flux')
                call rae_eqns%add_operator('RAE BC Flux')
                call rae_eqns%add_operator('RAE Volume Cylindrical Source')
                call rae_eqns%add_operator('RAE LaxFriedrichs Flux')

                call rae_eqns%add_model('RAE')

            case default
                call chidg_signal_one(FATAL, "build_rae: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************



end module eqn_rae
