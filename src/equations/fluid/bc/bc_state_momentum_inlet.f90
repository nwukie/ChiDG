module bc_state_momentum_inlet
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO
    use mod_fluid,              only: gam

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: momentum_inlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type momentum_inlet_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(momentum_inlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Momentum Inlet")
        call self%set_family("Inlet")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("Momentum-1")
!        call self%set_equation("Momentum-2")
!        call self%set_equation("Momentum-3")
!        call self%set_equation("Energy"    )


        !
        ! Add functions
        !
        call self%bcproperties%add('Density',   'Required')
        call self%bcproperties%add('X-Velocity','Required')
        call self%bcproperties%add('Y-Velocity','Required')
        call self%bcproperties%add('Z-Velocity','Required')



    end subroutine init
    !********************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(momentum_inlet_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM



        ! Equation indices


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,  p_m,    &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,   &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            flux_x, flux_y, flux_z, integrand,                  &
            u_m,    v_m,    w_m,                                &
            u_bc,   v_bc,   w_bc


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_field("Density"   , 'value', 'face interior')
        mom1_m    = worker%get_field("Momentum-1", 'value', 'face interior')
        mom2_m    = worker%get_field("Momentum-2", 'value', 'face interior')
        mom3_m    = worker%get_field("Momentum-3", 'value', 'face interior')
        energy_m  = worker%get_field("Energy"    , 'value', 'face interior')

        grad1_density_m = worker%get_field("Density"   , 'grad1', 'face interior')
        grad2_density_m = worker%get_field("Density"   , 'grad2', 'face interior')
        grad3_density_m = worker%get_field("Density"   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field("Momentum-1", 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_field("Momentum-1", 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_field("Momentum-1", 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_field("Momentum-2", 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_field("Momentum-2", 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field("Momentum-2", 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_field("Momentum-3", 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_field("Momentum-3", 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field("Momentum-3", 'grad3', 'face interior')
        
        grad1_energy_m  = worker%get_field("Energy"    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_field("Energy"    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_field("Energy"    , 'grad3', 'face interior')


        !
        ! Initialize variables
        !
        density_bc = density_m
        u_bc = density_m
        v_bc = density_m
        w_bc = density_m


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        density_bc = self%bcproperties%compute("Density",    worker%time(), worker%coords())
        u_bc       = self%bcproperties%compute("X-Velocity", worker%time(), worker%coords())
        v_bc       = self%bcproperties%compute("Y-Velocity", worker%time(), worker%coords())
        w_bc       = self%bcproperties%compute("Z-Velocity", worker%time(), worker%coords())


        !
        ! Compute bc momentum
        !
        mom1_bc = density_bc * u_bc
        mom2_bc = density_bc * v_bc
        mom3_bc = density_bc * w_bc

        !
        ! Compute interior pressure
        !
        p_m = worker%get_field('Pressure', 'value', 'face interior')



        !
        ! Compute bc energy
        !
        energy_bc = p_m/(gam - ONE) + (density_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )



        !
        ! Store computed boundary state
        !
        call worker%store_bc_state("Density"   , density_bc, 'value')
        call worker%store_bc_state("Momentum-1", mom1_bc,    'value')
        call worker%store_bc_state("Momentum-2", mom2_bc,    'value')
        call worker%store_bc_state("Momentum-3", mom3_bc,    'value')
        call worker%store_bc_state("Energy"    , energy_bc,  'value')




        call worker%store_bc_state("Density"   , grad1_density_m, 'grad1')
        call worker%store_bc_state("Density"   , grad2_density_m, 'grad2')
        call worker%store_bc_state("Density"   , grad3_density_m, 'grad3')

        call worker%store_bc_state("Momentum-1", grad1_mom1_m,    'grad1')
        call worker%store_bc_state("Momentum-1", grad2_mom1_m,    'grad2')
        call worker%store_bc_state("Momentum-1", grad3_mom1_m,    'grad3')
                                                
        call worker%store_bc_state("Momentum-2", grad1_mom2_m,    'grad1')
        call worker%store_bc_state("Momentum-2", grad2_mom2_m,    'grad2')
        call worker%store_bc_state("Momentum-2", grad3_mom2_m,    'grad3')
                                                
        call worker%store_bc_state("Momentum-3", grad1_mom3_m,    'grad1')
        call worker%store_bc_state("Momentum-3", grad2_mom3_m,    'grad2')
        call worker%store_bc_state("Momentum-3", grad3_mom3_m,    'grad3')
                                                
        call worker%store_bc_state("Energy"    , grad1_energy_m,  'grad1')
        call worker%store_bc_state("Energy"    , grad2_energy_m,  'grad2')
        call worker%store_bc_state("Energy"    , grad3_energy_m,  'grad3')










    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_momentum_inlet
