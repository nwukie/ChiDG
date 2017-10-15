module eqn_rac
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
    type, public, extends(equation_builder_t) :: rac

    contains

        procedure   :: init
        procedure   :: build

    end type rac
    !*************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self)
        class(rac),   intent(inout)  :: self

        call self%set_name('RAC')

    end subroutine init
    !*************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    function build(self,blueprint) result(rac_eqns)
        class(rac),     intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rac_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call rac_eqns%set_name('RAC')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call rac_eqns%add_operator('RAC Volume Flux')
                call rac_eqns%add_operator('RAC Boundary Average Flux')
                call rac_eqns%add_operator('RAC BC Flux')
                call rac_eqns%add_operator('RAC Volume Cylindrical Source')
                call rac_eqns%add_operator('RAC LaxFriedrichs Flux')

                call rac_eqns%add_model('RAC')

            case default
                call chidg_signal_one(FATAL, "build_rac: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !***********************************************************************************



end module eqn_rac
