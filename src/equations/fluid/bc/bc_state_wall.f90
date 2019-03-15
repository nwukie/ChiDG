module bc_state_wall
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ONE
    use mod_fluid,              only: gam, Rgas
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name('Wall')
        call self%set_family('Wall')


!        !
!        ! Set operator equations
!        !
!        call self%set_equation('Density'   )
!        call self%set_equation('Momentum-1')
!        call self%set_equation('Momentum-2')
!        call self%set_equation('Momentum-3')
!        call self%set_equation('Energy'    )


    end subroutine init
    !********************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(wall_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            grad1_density_bc, grad1_mom1_bc, grad1_mom2_bc, grad1_mom3_bc, grad1_energy_bc,  &
            grad2_density_bc, grad2_mom1_bc, grad2_mom2_bc, grad2_mom3_bc, grad2_energy_bc,  &
            grad3_density_bc, grad3_mom1_bc, grad3_mom2_bc, grad3_mom3_bc, grad3_energy_bc,  &
            grad1_vel1, grad2_vel1, grad3_vel1, &
            grad1_vel2, grad2_vel2, grad3_vel2, &
            grad1_vel3, grad2_vel3, grad3_vel3, &
            grad1_p, grad2_p, grad3_p, &
            dvel1_ddensity, dvel1_dmom1, &
            dvel2_ddensity, dvel2_dmom2, &
            dvel3_ddensity, dvel3_dmom3, &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy, &
            vel1_m, vel2_m, vel3_m, p_m, T_m, invdensity

        real(rk),   allocatable, dimension(:)   :: r
        real(rk),   allocatable, dimension(:,:) :: grid_velocity




        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_field('Density'    , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1' , 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2' , 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3' , 'value', 'face interior')
        energy_m  = worker%get_field('Energy'     , 'value', 'face interior')
        p_m       = worker%get_field('Pressure'   , 'value', 'face interior')
        T_m       = worker%get_field('Temperature', 'value', 'face interior')


    
        grad1_density_m = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m    = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field('Momentum-3', 'grad3', 'face interior')

        grad1_energy_m  = worker%get_field('Energy'    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_field('Energy'    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_field('Energy'    , 'grad3', 'face interior')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)
        end if



        ! Initialize arrays
        vel1_m = mom1_m/density_m
        vel2_m = mom2_m/density_m
        vel3_m = mom3_m/density_m

        density_bc = density_m
        mom1_bc    = density_m
        mom2_bc    = density_m
        mom3_bc    = density_m
        energy_bc  = density_m



        !
        ! Set relative velocity to zero: ie set fluid velocity equal to grid/wall velocity.
        !
        grid_velocity = worker%get_grid_velocity_face('face interior')
        mom1_bc = density_m*grid_velocity(:,1)
        mom2_bc = density_m*grid_velocity(:,2)
        mom3_bc = density_m*grid_velocity(:,3)


        !
        ! Compute energy by extrapolating pressure and adding the fluid velocity on the wall
        !
        energy_bc = p_m/(gam-ONE) + HALF*density_m*(grid_velocity(:,1)**TWO + grid_velocity(:,2)**TWO + grid_velocity(:,3)**TWO)
        !energy_bc = energy_m - HALF*density_m*(vel1_m**TWO + vel2_m**TWO + vel3_m**TWO)


        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')








        !
        ! Impose adiabatic
        !


        !
        ! compute velocity jacobians
        !
        invdensity = ONE/density_m
        dvel1_ddensity = -invdensity*invdensity*mom1_m
        dvel2_ddensity = -invdensity*invdensity*mom2_m
        dvel3_ddensity = -invdensity*invdensity*mom3_m

        dvel1_dmom1 = invdensity
        dvel2_dmom2 = invdensity
        dvel3_dmom3 = invdensity



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_vel1 = dvel1_ddensity*grad1_density_m  +  dvel1_dmom1*grad1_mom1_m
        grad2_vel1 = dvel1_ddensity*grad2_density_m  +  dvel1_dmom1*grad2_mom1_m
        grad3_vel1 = dvel1_ddensity*grad3_density_m  +  dvel1_dmom1*grad3_mom1_m

        grad1_vel2 = dvel2_ddensity*grad1_density_m  +  dvel2_dmom2*grad1_mom2_m
        grad2_vel2 = dvel2_ddensity*grad2_density_m  +  dvel2_dmom2*grad2_mom2_m
        grad3_vel2 = dvel2_ddensity*grad3_density_m  +  dvel2_dmom2*grad3_mom2_m

        grad1_vel3 = dvel3_ddensity*grad1_density_m  +  dvel3_dmom3*grad1_mom3_m
        grad2_vel3 = dvel3_ddensity*grad2_density_m  +  dvel3_dmom3*grad2_mom3_m
        grad3_vel3 = dvel3_ddensity*grad3_density_m  +  dvel3_dmom3*grad3_mom3_m

        !
        ! Pressure jacobians
        !
        dp_ddensity =  (gam-ONE)*HALF*(mom1_m*mom1_m + mom2_m*mom2_m + mom3_m*mom3_m)/(density_m*density_m)
        dp_dmom1    = -(gam-ONE)*mom1_m/density_m
        dp_dmom2    = -(gam-ONE)*mom2_m/density_m
        dp_dmom3    = -(gam-ONE)*mom3_m/density_m
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)



        ! Compute pressure gradient using Chain-rule
        grad1_p = dp_ddensity * grad1_density_m  + &
                  dp_dmom1    * grad1_mom1_m     + &
                  dp_dmom2    * grad1_mom2_m     + &
                  dp_dmom3    * grad1_mom3_m     + &
                  dp_denergy  * grad1_energy_m

        grad2_p = dp_ddensity * grad2_density_m  + &
                  dp_dmom1    * grad2_mom1_m     + &
                  dp_dmom2    * grad2_mom2_m     + &
                  dp_dmom3    * grad2_mom3_m     + &
                  dp_denergy  * grad2_energy_m

        grad3_p = dp_ddensity * grad3_density_m  + &
                  dp_dmom1    * grad3_mom1_m     + &
                  dp_dmom2    * grad3_mom2_m     + &
                  dp_dmom3    * grad3_mom3_m     + &
                  dp_denergy  * grad3_energy_m



        grad1_density_bc = (ONE/(Rgas*T_m))*grad1_p
        grad2_density_bc = (ONE/(Rgas*T_m))*grad2_p
        grad3_density_bc = (ONE/(Rgas*T_m))*grad3_p


        grad1_mom1_bc = (p_m/(Rgas*T_m))*grad1_vel1  +  (grid_velocity(:,1)/(Rgas*T_m))*grad1_p
        grad2_mom1_bc = (p_m/(Rgas*T_m))*grad2_vel1  +  (grid_velocity(:,1)/(Rgas*T_m))*grad2_p
        grad3_mom1_bc = (p_m/(Rgas*T_m))*grad3_vel1  +  (grid_velocity(:,1)/(Rgas*T_m))*grad3_p

        grad1_mom2_bc = (p_m/(Rgas*T_m))*grad1_vel2  +  (grid_velocity(:,2)/(Rgas*T_m))*grad1_p
        grad2_mom2_bc = (p_m/(Rgas*T_m))*grad2_vel2  +  (grid_velocity(:,2)/(Rgas*T_m))*grad2_p
        grad3_mom2_bc = (p_m/(Rgas*T_m))*grad3_vel2  +  (grid_velocity(:,2)/(Rgas*T_m))*grad3_p

        grad1_mom3_bc = (p_m/(Rgas*T_m))*grad1_vel3  +  (grid_velocity(:,3)/(Rgas*T_m))*grad1_p
        grad2_mom3_bc = (p_m/(Rgas*T_m))*grad2_vel3  +  (grid_velocity(:,3)/(Rgas*T_m))*grad2_p
        grad3_mom3_bc = (p_m/(Rgas*T_m))*grad3_vel3  +  (grid_velocity(:,3)/(Rgas*T_m))*grad3_p

        grad1_energy_bc = (p_m*grid_velocity(:,1)/(Rgas*T_m))*grad1_vel1 +  &
                          (p_m*grid_velocity(:,2)/(Rgas*T_m))*grad1_vel2 +  &
                          (p_m*grid_velocity(:,3)/(Rgas*T_m))*grad1_vel3 +  &
                          (ONE/(gam-ONE) + HALF*(ONE/(Rgas*T_m))*(grid_velocity(:,1)**TWO + grid_velocity(:,2)**TWO + grid_velocity(:,3)**TWO))*grad1_p

        grad2_energy_bc = (p_m*grid_velocity(:,1)/(Rgas*T_m))*grad2_vel1 +  &
                          (p_m*grid_velocity(:,2)/(Rgas*T_m))*grad2_vel2 +  &
                          (p_m*grid_velocity(:,3)/(Rgas*T_m))*grad2_vel3 +  &
                          (ONE/(gam-ONE) + HALF*(ONE/(Rgas*T_m))*(grid_velocity(:,1)**TWO + grid_velocity(:,2)**TWO + grid_velocity(:,3)**TWO))*grad2_p

        grad3_energy_bc = (p_m*grid_velocity(:,1)/(Rgas*T_m))*grad3_vel1 +  &
                          (p_m*grid_velocity(:,2)/(Rgas*T_m))*grad3_vel2 +  &
                          (p_m*grid_velocity(:,3)/(Rgas*T_m))*grad3_vel3 +  &
                          (ONE/(gam-ONE) + HALF*(ONE/(Rgas*T_m))*(grid_velocity(:,1)**TWO + grid_velocity(:,2)**TWO + grid_velocity(:,3)**TWO))*grad3_p


        ! Convert gradients of tangential momentum to angular momentum
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            grad1_mom2_bc = r*grad1_mom2_bc + density_bc*grid_velocity(:,2) !mom2_bc
            grad2_mom2_bc = r*grad2_mom2_bc
            grad3_mom2_bc = r*grad3_mom2_bc
        end if


        call worker%store_bc_state('Density'   , grad1_density_bc, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density_bc, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density_bc, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', grad1_mom1_bc, 'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1_bc, 'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1_bc, 'grad3')
                                                
        call worker%store_bc_state('Momentum-2', grad1_mom2_bc, 'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2_bc, 'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2_bc, 'grad3')
                                                
        call worker%store_bc_state('Momentum-3', grad1_mom3_bc, 'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3_bc, 'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3_bc, 'grad3')

        call worker%store_bc_state('Energy'    , grad1_energy_bc, 'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy_bc, 'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy_bc, 'grad3')


















!        grad1_density_m = ZERO
!        grad2_density_m = ZERO
!        grad3_density_m = ZERO
!        call worker%store_bc_state('Density'   , grad1_density_m, 'grad1')
!        call worker%store_bc_state('Density'   , grad2_density_m, 'grad2')
!        call worker%store_bc_state('Density'   , grad3_density_m, 'grad3')
!                                                
!        call worker%store_bc_state('Momentum-1', grad1_mom1_m, 'grad1')
!        call worker%store_bc_state('Momentum-1', grad2_mom1_m, 'grad2')
!        call worker%store_bc_state('Momentum-1', grad3_mom1_m, 'grad3')
!                                                
!        call worker%store_bc_state('Momentum-2', grad1_mom2_m, 'grad1')
!        call worker%store_bc_state('Momentum-2', grad2_mom2_m, 'grad2')
!        call worker%store_bc_state('Momentum-2', grad3_mom2_m, 'grad3')
!                                                
!        call worker%store_bc_state('Momentum-3', grad1_mom3_m, 'grad1')
!        call worker%store_bc_state('Momentum-3', grad2_mom3_m, 'grad2')
!        call worker%store_bc_state('Momentum-3', grad3_mom3_m, 'grad3')
!
!        grad1_energy_m = ZERO
!        grad2_energy_m = ZERO
!        grad3_energy_m = ZERO
!        call worker%store_bc_state('Energy'    , grad1_energy_m, 'grad1')
!        call worker%store_bc_state('Energy'    , grad2_energy_m, 'grad2')
!        call worker%store_bc_state('Energy'    , grad3_energy_m, 'grad3')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_wall
