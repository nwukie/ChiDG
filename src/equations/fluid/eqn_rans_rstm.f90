module eqn_rans_rstm
#include <messenger.h>
    use type_equation_set,          only: equation_set_t
    use type_equation_builder,      only: equation_builder_t
    use type_fluid_pseudo_timestep, only: fluid_pseudo_timestep_t
    implicit none

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(equation_builder_t) :: rans_rstm

    contains

        procedure   :: init
        procedure   :: build

    end type rans_rstm
    !*******************************************************************************************




contains


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_rstm),   intent(inout)  :: self

        call self%set_name('RANS_RSTM')

    end subroutine init
    !*******************************************************************************************



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    function build(self,blueprint) result(rans_rstm_eqns)
        class(rans_rstm),    intent(in)  :: self
        character(*),   intent(in)  :: blueprint

        type(equation_set_t)            :: rans_rstm_eqns
        type(fluid_pseudo_timestep_t)   :: fluid_pseudo_time

        !
        ! Set equation set name
        !
        call rans_rstm_eqns%set_name('RANS_RSTM')
        

        !
        ! Add spatial operators
        !
        select case (trim(blueprint))


            case('default')
                call rans_rstm_eqns%add_operator('Euler Volume Flux')
                call rans_rstm_eqns%add_operator('Euler Boundary Average Flux')
                call rans_rstm_eqns%add_operator('Euler Roe Flux')
                call rans_rstm_eqns%add_operator('Euler BC Flux')
                call rans_rstm_eqns%add_operator('Euler Volume Cylindrical Source')

                call rans_rstm_eqns%add_operator('Fluid Viscous Volume Operator')
                call rans_rstm_eqns%add_operator('Fluid Viscous Boundary Average Operator')
                call rans_rstm_eqns%add_operator('Fluid Viscous BC Operator')
                call rans_rstm_eqns%add_operator('Fluid Viscous Volume Cylindrical Source')



! NATHAN COMMENTED THIS OUT DURING MERGE, BECAUSE IT WASN'T INCLUDED IN RSTM DIRECTORY
!                call rans_rstm_eqns%add_model('RSTM Turbulence Kinetic Energy')
!



                !call rans_rstm_eqns%add_model('Ideal Gas RSTM')
                call rans_rstm_eqns%add_model('Ideal Gas')
                call rans_rstm_eqns%add_model('Constant Viscosity')
                !call rans_rstm_eqns%add_model('Sutherlands Law')
                call rans_rstm_eqns%add_model('Stokes Hypothesis')
                call rans_rstm_eqns%add_model('Reynolds Analogy')
                !call rans_rstm_eqns%add_model('Zero Reynolds Stress')
                !call rans_rstm_eqns%add_model('Fluid Advection Velocity')
                !call rans_rstm_eqns%add_model('Constant Viscosity')
                call rans_rstm_eqns%add_model('Velocity Gradient')
                call rans_rstm_eqns%add_model('Strain Rate')
                call rans_rstm_eqns%add_model('Rotation Rate')
                call rans_rstm_eqns%add_model('Wall Distance : p-Poisson Normalization')
                call rans_rstm_eqns%add_model('Zero Turbulent Model Fields')
                
                
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Realizable Reynolds Stress')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Turbulence Quantities')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW LRR Coefficients')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Production')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Isotropic Dissipation')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Pressure-Strain Correlation')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Realizability Source')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Artificial Viscosity')

                ! Add RSTM diffusion model
                ! 'RSTMSSGLRRW Generalized Diffusion' is the standard model
                ! 'RSTMSSGLRRW Simple Diffusion' can be used for greater numerical robustness
                !call rans_rstm_eqns%add_model('RSTMSSGLRRW Generalized Diffusion')
                call rans_rstm_eqns%add_model('RSTMSSGLRRW Simple Diffusion')

                !call rans_rstm_eqns%add_operator('Spalart-Allmaras Source Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras Advection Boundary Average Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras LaxFriedrichs Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras Volume Advection Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras BC Advection Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras Boundary Diffusion Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras Volume Diffusion Operator')
                !call rans_rstm_eqns%add_operator('Spalart-Allmaras BC Diffusion Operator')




                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Source Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Advection Boundary Average Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW LaxFriedrichs Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Volume Advection Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW BC Advection Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Boundary Diffusion Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Volume Diffusion Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW BC Diffusion Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Artificial Viscosity Volume Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Artificial Viscosity BC Operator')
                call rans_rstm_eqns%add_operator('RSTMSSGLRRW Artificial Viscosity Boundary Average Operator')


                
                

                ! Add shear stress after turbulence viscosity models from SA so they are computed first
                call rans_rstm_eqns%add_model('Reynolds Shear Stress')
                call rans_rstm_eqns%add_model('Temperature Gradient')


                call rans_rstm_eqns%add_pseudo_timestep(fluid_pseudo_time)


            case default
                call chidg_signal_one(FATAL, "build_rans_rstm: I didn't recognize the &
                                              construction parameter that was passed to build &
                                              the equation set.", blueprint)

        end select


    end function build
    !******************************************************************************************






end module eqn_rans_rstm
