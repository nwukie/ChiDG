module eqn_laminar_navier_stokes_io
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
    type, public, extends(equation_builder_t) :: laminar_navier_stokes_io

    contains

        procedure   :: init
        procedure   :: build

    end type laminar_navier_stokes_io
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(laminar_navier_stokes_io),   intent(inout)  :: self

        call self%set_name('Laminar Navier Stokes IO')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(eqn)
        class(laminar_navier_stokes_io),    intent(in)  :: self
        character(*),                       intent(in)  :: blueprint

        type(equation_set_t)    :: eqn

        !
        ! Set equation set name
        !
        call eqn%set_name(self%get_name())
        


        call eqn%add_operator('Euler Volume Flux')
        call eqn%add_model('Ideal Gas')


        call eqn%prop%clear_io_fields()
        call eqn%prop%add_io_field('Pressure')



        !call eqn%add_model('Fluid Advection Velocity')
        !call eqn%add_model('Constant Viscosity')
        !call eqn%add_model('Stokes Hypothesis')
        !call eqn%add_model('Reynolds Analogy')
        !call eqn%add_model('Zero Turbulent Model Fields')
        !call eqn%add_model('Shear Stress')
        !call eqn%add_model('Temperature Gradient')




    end function build
    !**********************************************************************************************






end module eqn_laminar_navier_stokes_io
