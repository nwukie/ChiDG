module eqn_rans
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    implicit none


    !>
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: rans

    contains

        procedure   :: init
        procedure   :: build

    end type rans
    !*******************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self)
        class(rans),   intent(inout)  :: self

        call self%set_name('RANS')

    end subroutine init
    !*******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function build(self,blueprint) result(rans_eqns)
        class(rans),    intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rans_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call rans_eqns%set_name('RANS')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call rans_eqns%add_operator('Euler Volume Flux')
                call rans_eqns%add_operator('Euler Boundary Average Flux')
                call rans_eqns%add_operator('Euler Roe Flux')
                call rans_eqns%add_operator('Euler BC Flux')
                call rans_eqns%add_operator('Euler Volume Cylindrical Source')

                call rans_eqns%add_operator('Fluid Viscous Volume Operator')
                call rans_eqns%add_operator('Fluid Viscous Boundary Average Operator')
                call rans_eqns%add_operator('Fluid Viscous BC Operator')
                call rans_eqns%add_operator('Fluid Viscous Volume Cylindrical Source')

                call rans_eqns%add_model('Ideal Gas')
                !call rans_eqns%add_model('Sutherlands Law')
                call rans_eqns%add_model('Constant Viscosity')
                call rans_eqns%add_model('Stokes Hypothesis')
                call rans_eqns%add_model('Reynolds Analogy')


                call rans_eqns%add_operator('Spalart-Allmaras Source Operator')
                call rans_eqns%add_operator('Spalart-Allmaras Advection Boundary Average Operator')
                call rans_eqns%add_operator('Spalart-Allmaras LaxFriedrichs Operator')
                call rans_eqns%add_operator('Spalart-Allmaras Volume Advection Operator')
                call rans_eqns%add_operator('Spalart-Allmaras BC Advection Operator')
                call rans_eqns%add_operator('Spalart-Allmaras Boundary Diffusion Operator')
                call rans_eqns%add_operator('Spalart-Allmaras Volume Diffusion Operator')
                call rans_eqns%add_operator('Spalart-Allmaras BC Diffusion Operator')


                ! Add shear stress after turbulence viscosity models from SA so they are computed first
                call rans_eqns%add_model('Shear Stress')
                call rans_eqns%add_model('Temperature Gradient')


                call rans_eqns%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_rans: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !******************************************************************************************






end module eqn_rans
