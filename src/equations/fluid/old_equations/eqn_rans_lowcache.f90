module eqn_rans_lowcache
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
    type, public, extends(equation_builder_t) :: rans_lowcache

    contains

        procedure   :: init
        procedure   :: build

    end type rans_lowcache
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_lowcache),   intent(inout)  :: self

        call self%set_name('RANS Low-Cache')

    end subroutine init
    !*********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function build(self,blueprint) result(equation_set)
        class(rans_lowcache),   intent(in)  :: self
        character(*),           intent(in)  :: blueprint

        type(equation_set_t)            :: equation_set
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call equation_set%set_name('RANS Low-Cache')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call equation_set%add_operator('RANS Volume Advection')
                call equation_set%add_operator('RANS Volume Diffusion')
                call equation_set%add_operator('RANS Boundary Advection')
                call equation_set%add_operator('RANS Boundary Diffusion')
                call equation_set%add_operator('RANS BC Advection')
                call equation_set%add_operator('RANS BC Diffusion')
                call equation_set%add_operator('RANS Source')

                call equation_set%add_model('Ideal Gas')
                !call equation_set%add_model('Fluid Advection Velocity')
                call equation_set%add_model('Constant Viscosity')

                call equation_set%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_rans_lowcache: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !**********************************************************************************************






end module eqn_rans_lowcache
