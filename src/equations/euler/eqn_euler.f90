module eqn_euler
#include <messenger.h>
    use type_equation_set,      only: equation_set_t
    use type_equation_builder,  only: equation_builder_t
    use perfect_gas,            only: perfect_gas_t
    implicit none


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: euler

    contains

        procedure   :: init
        procedure   :: build

    end type euler
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(euler),   intent(inout)  :: self

        call self%set_name('Euler')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(euler_eqns)
        class(euler),       intent(in)  :: self
        character(len=*),   intent(in)  :: blueprint

        type(equation_set_t)    :: euler_eqns
        type(perfect_gas_t)     :: perfect_gas

        !
        ! Set equation set name
        !
        call euler_eqns%set_name("Euler")
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call euler_eqns%add_operator('Euler Volume Flux')
                call euler_eqns%add_operator('Euler Boundary Average Flux')
                call euler_eqns%add_operator('Euler Roe Flux')

                call euler_eqns%prop%add_fluid(perfect_gas)



            case default
                call chidg_signal_one(FATAL, "build_euler: I didn't recognize the construction parameter &
                                              that was passed to build the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************






end module eqn_euler
