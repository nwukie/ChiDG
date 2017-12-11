module eqn_pgradtest
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
    type, public, extends(equation_builder_t) :: pgradtest

    contains

        procedure   :: init
        procedure   :: build

    end type pgradtest
    !*************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self)
        class(pgradtest),   intent(inout)  :: self

        call self%set_name('pgradtest')

    end subroutine init
    !*************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    function build(self,blueprint) result(pgradtest_eqns)
        class(pgradtest),    intent(in)  :: self
        character(*),       intent(in)  :: blueprint

        type(equation_set_t)            :: pgradtest_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call pgradtest_eqns%set_name('pgradtest')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')

                call pgradtest_eqns%add_operator('pgradtest BC Operator')
                call pgradtest_eqns%add_operator('pgradtest Boundary Average Operator')
                call pgradtest_eqns%add_operator('pgradtest Volume Operator')

                !call pgradtest_eqns%add_model('Ideal Gas')
                !call pgradtest_eqns%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_pgradtest: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !***********************************************************************************



end module eqn_pgradtest
