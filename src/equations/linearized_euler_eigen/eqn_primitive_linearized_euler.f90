module eqn_primitive_linearized_euler
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_solverdata,            only: solverdata_t
!    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    implicit none


    !>
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: primitive_linearized_euler

    contains

        procedure   :: init
        procedure   :: build

    end type primitive_linearized_euler
    !************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self)
        class(primitive_linearized_euler),   intent(inout)  :: self

        call self%set_name('Linearized Euler Eigen')

    end subroutine init
    !*************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    function build(self,blueprint) result(primitive_linearized_euler_eqns)
        class(primitive_linearized_euler),   intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: primitive_linearized_euler_eqns

        !
        ! Set equation set name
        !
        call primitive_linearized_euler_eqns%set_name('Linearized Euler Eigen')
        

        !
        ! Add operators
        !
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Boundary Average')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER LaxFriedrichs')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Volume Advection')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER BC')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Temporal Source')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Axial Source')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Circumferential Source')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Equation Source')
        call primitive_linearized_euler_eqns%add_operator('PRIMLINEULER Divergence Source')


    end function build
    !**********************************************************************************










end module eqn_primitive_linearized_euler
