module eqn_rans_sst
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
    type, public, extends(equation_builder_t) :: rans_sst

    contains

        procedure   :: init
        procedure   :: build

    end type rans_sst
    !*******************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_sst),   intent(inout)  :: self

        call self%set_name('RANS SST')

    end subroutine init
    !*******************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function build(self,blueprint) result(rans_sst_eqns)
        class(rans_sst),    intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rans_sst_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call rans_sst_eqns%set_name('RANS SST')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call rans_sst_eqns%add_operator('Euler Volume Flux')
                call rans_sst_eqns%add_operator('Euler Boundary Average Flux')
                call rans_sst_eqns%add_operator('Euler Roe Flux')
                call rans_sst_eqns%add_operator('Euler BC Flux')
                call rans_sst_eqns%add_operator('Euler Volume Cylindrical Source')

                call rans_sst_eqns%add_operator('Fluid Viscous Volume Operator')
                call rans_sst_eqns%add_operator('Fluid Viscous Boundary Average Operator')
                call rans_sst_eqns%add_operator('Fluid Viscous BC Operator')
                call rans_sst_eqns%add_operator('Fluid Viscous Volume Cylindrical Source')

                call rans_sst_eqns%add_model('SST Turbulence Kinetic Energy')
                call rans_sst_eqns%add_model('Ideal Gas SST')
                call rans_sst_eqns%add_model('Constant Viscosity')
                !call rans_sst_eqns%add_model('Sutherlands Law')
                call rans_sst_eqns%add_model('Wall Distance : p-Poisson Normalization')
                call rans_sst_eqns%add_model('Velocity Gradient')
                call rans_sst_eqns%add_model('Velocity Divergence and Curl')
                call rans_sst_eqns%add_model('Rotation Rate')
                call rans_sst_eqns%add_model('Strain Rate')
                call rans_sst_eqns%add_model('Stokes Hypothesis')
                call rans_sst_eqns%add_model('Reynolds Analogy')
                call rans_sst_eqns%add_model('Zero Reynolds Stress')
                !call rans_sst_eqns%add_model('Fluid Advection Velocity')
!                call rans_sst_eqns%add_model('Zero Turbulent Model Fields')
                call rans_sst_eqns%add_model('SST Turbulence Quantities')
                call rans_sst_eqns%add_model('SST Blended Coefficients')
                call rans_sst_eqns%add_model('SST Source Terms')





                call rans_sst_eqns%add_operator('SST Source Operator')
                call rans_sst_eqns%add_operator('SST Advection Boundary Average Operator')
                !call rans_sst_eqns%add_operator('SST Roe Flux')
                call rans_sst_eqns%add_operator('SST LaxFriedrichs Operator')
                call rans_sst_eqns%add_operator('SST Volume Advection Operator')
                call rans_sst_eqns%add_operator('SST BC Advection Operator')
                call rans_sst_eqns%add_operator('SST Boundary Diffusion Operator')
                call rans_sst_eqns%add_operator('SST Volume Diffusion Operator')
                call rans_sst_eqns%add_operator('SST BC Diffusion Operator')
                call rans_sst_eqns%add_operator('SST Artificial Viscosity Volume Operator')
                call rans_sst_eqns%add_operator('SST Artificial Viscosity BC Operator')
                call rans_sst_eqns%add_operator('SST Artificial Viscosity Boundary Average Operator')

                ! Add shear stress after turbulence viscosity models from SA so they are computed first
                call rans_sst_eqns%add_model('Shear Stress')
                call rans_sst_eqns%add_model('Temperature Gradient')




                call rans_sst_eqns%add_pseudo_timestep(fluid_pseudo_time)


                !call rans_sst_eqns%add_io_field('SST k Source Term')
                !call rans_sst_eqns%add_io_field('SST Omega Source Term')
                !call rans_sst_eqns%add_io_field('SST sigma_d')
                !call rans_sst_eqns%add_io_field('Laminar Viscosity')
                !call rans_sst_eqns%add_io_field('Velocity 1 - Gradient 2')
                call rans_sst_eqns%add_io_field('Turbulent Viscosity')
                !call rans_sst_eqns%add_io_field('k')
                !call rans_sst_eqns%add_io_field('SST k Production Term')
                !call rans_sst_eqns%add_io_field('Wall Distance')
            case default
                call chidg_signal_one(FATAL, "build_rans_sst: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !******************************************************************************************






end module eqn_rans_sst
