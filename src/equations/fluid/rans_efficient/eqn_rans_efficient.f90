module eqn_rans_efficient
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    use mod_rans_efficient,         only: turbulence_model, viscosity_model, state_equation, fluid
    implicit none


    !>
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: rans_efficient

    contains

        procedure   :: init
        procedure   :: build

    end type rans_efficient
    !*******************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_efficient),   intent(inout)  :: self

        call self%set_name('RANS Efficient')

    end subroutine init
    !*******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function build(self,blueprint) result(rans_efficient_eqns)
        class(rans_efficient),    intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rans_efficient_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time


        integer         :: unit, msg
        logical         :: file_exists



        ! Set equation set name
        call rans_efficient_eqns%set_name('RANS Efficient')
        
        ! Add spatial operators
        select case (trim(blueprint))

            case('default')


                ! Check if input from 'models.nml' is available.
                !   1: if available, read 
                !   2: if not available, do nothing and turbulence_model retains default value
                inquire(file='models.nml', exist=file_exists)
                if (file_exists) then
                    open(newunit=unit,form='formatted',file='models.nml')
                    read(unit,nml=fluid,iostat=msg)
                    close(unit)
                end if



                ! Add models dependent upon namelist input in models.nml
                call rans_efficient_eqns%add_model(trim(state_equation))
                call rans_efficient_eqns%add_model(trim(viscosity_model))


                call rans_efficient_eqns%add_operator('RANS Boundary Average Advection')
                call rans_efficient_eqns%add_operator('RANS Boundary Average Diffusion')
                call rans_efficient_eqns%add_operator('RANS BC Advection')
                call rans_efficient_eqns%add_operator('RANS BC Diffusion')
                call rans_efficient_eqns%add_operator('RANS Volume Advection')
                call rans_efficient_eqns%add_operator('RANS Volume Diffusion')
                call rans_efficient_eqns%add_operator('RANS Upwind Operator')
                call rans_efficient_eqns%add_operator('RANS Volume Source')

                call rans_efficient_eqns%add_pseudo_timestep(fluid_pseudo_time)

            case default
                call chidg_signal_one(FATAL, "build_rans_efficient: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !******************************************************************************************






end module eqn_rans_efficient
