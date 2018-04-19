module eqn_empty
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_solverdata,            only: solverdata_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use mod_operators,              only: operator_t, operator_factory
    implicit none


    !>
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: empty

    contains

        procedure   :: init => init_eqn
        procedure   :: build

    end type empty
    !************************************************************************************

    !>
    !!  @author Nathan A. Wukie
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: empty_volume_advective_operator_t

    contains
        procedure   :: init => init_op
        procedure   :: compute
    end type empty_volume_advective_operator_t
    !*************************************************************************



contains

    !>
    !!  @author Nathan A. Wukie
    !--------------------------------------------------------------------------------
    subroutine init_op(self)
        class(empty_volume_advective_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Empty Advection Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')

    end subroutine init_op
    !********************************************************************************

    !>
    !!  @author Nathan A. Wukie
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(empty_volume_advective_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

    end subroutine compute
    !****************************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/19/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine init_eqn(self)
        class(empty),   intent(inout)  :: self

        call self%set_name('Empty')

    end subroutine init_eqn
    !*************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/19/2018
    !!
    !-------------------------------------------------------------------------------------
    function build(self,blueprint) result(empty_eqn)
        class(empty),   intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: empty_eqn
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time
        type(empty_volume_advective_operator_t)                  :: empty_op

        call operator_factory%register(empty_op)

        ! Set equation set name
        call empty_eqn%set_name('Empty')

        ! Add spatial operators
        select case (trim(blueprint))

            case('default')
                call empty_eqn%add_operator('Empty Advection Volume Operator')

            case default
                call chidg_signal_one(FATAL, "build_empty: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************

















end module eqn_empty
