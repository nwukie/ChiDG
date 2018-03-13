module bc_state_auxiliary_interior
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO
    use mod_fluid,              only: gam
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: auxiliary_interior_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type auxiliary_interior_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(auxiliary_interior_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("Auxiliary Gradient Interior")
        call self%set_family("Wall")

    end subroutine init
    !********************************************************************************





    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !---------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(auxiliary_interior_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                  &
            density_o,  grad1_density_o,    grad2_density_o,    grad3_density_o,    &
            vel1_o,     grad1_vel1_o,       grad2_vel1_o,       grad3_vel1_o,       &
                        grad1_vel2_o,       grad2_vel2_o,       grad3_vel2_o,       &
                        grad1_vel3_o,       grad2_vel3_o,       grad3_vel3_o,       &
            p_o,        grad1_p_o,          grad2_p_o,          grad3_p_o,          &
            density_bc, grad1_density_bc,   grad2_density_bc,   grad3_density_bc,   &
            vel1_bc,    grad1_vel1_bc,      grad2_vel1_bc,      grad3_vel1_bc,      &
            p_bc,       grad1_p_bc,         grad2_p_bc,         grad3_p_bc,         &
            mom1_o
            


        !
        ! Retrieve values from original problem
        !
        density_o = worker%get_field('Density',    'value', 'face interior')
        mom1_o    = worker%get_field('Momentum-1', 'value', 'face interior')
        p_o       = worker%get_field('Pressure',   'value', 'face interior')
        vel1_o    = mom1_o/density_o


        call compute_pressure_gradient(worker,grad1_p_o, grad2_p_o, grad3_p_o)
        call compute_velocity_gradients(worker, &
                                        grad1_vel1_o,grad2_vel1_o,grad3_vel1_o,  &
                                        grad1_vel2_o,grad2_vel2_o,grad3_vel2_o,  &
                                        grad1_vel3_o,grad2_vel3_o,grad3_vel3_o)
        call compute_density_gradient(worker,grad1_density_o, grad2_density_o, grad3_density_o)




        !
        ! Set boundary condition on auxiliary problem from original problem
        !
        density_bc       = density_o
        grad1_density_bc = grad1_density_o
        grad2_density_bc = grad2_density_o
        grad3_density_bc = grad3_density_o

        vel1_bc       = vel1_o
        grad1_vel1_bc = grad1_vel1_o
        grad2_vel1_bc = grad2_vel1_o
        grad3_vel1_bc = grad3_vel1_o

        p_bc       = p_o
        grad1_p_bc = grad1_p_o
        grad2_p_bc = grad2_p_o
        grad3_p_bc = grad3_p_o



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density_TEMP',    density_bc,       'value')
        call worker%store_bc_state('Density_TEMP',    grad1_density_bc, 'grad1')
        call worker%store_bc_state('Density_TEMP',    grad2_density_bc, 'grad2')
        call worker%store_bc_state('Density_TEMP',    grad3_density_bc, 'grad3')

        call worker%store_bc_state('Velocity-1_TEMP', vel1_bc,          'value')
        call worker%store_bc_state('Velocity-1_TEMP', grad1_vel1_bc,    'grad1')
        call worker%store_bc_state('Velocity-1_TEMP', grad2_vel1_bc,    'grad2')
        call worker%store_bc_state('Velocity-1_TEMP', grad3_vel1_bc,    'grad3')

        call worker%store_bc_state('Pressure_TEMP',   p_bc,             'value')
        call worker%store_bc_state('Pressure_TEMP',   grad1_p_bc,       'grad1')
        call worker%store_bc_state('Pressure_TEMP',   grad2_p_bc,       'grad2')
        call worker%store_bc_state('Pressure_TEMP',   grad3_p_bc,       'grad3')




    end subroutine compute_bc_state
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_pressure_gradient(worker,grad1_p, grad2_p, grad3_p) 
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_p 
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_p 
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_p 

        type(AD_D), allocatable, dimension(:)   ::                              &
            density,       mom1,       mom2,       mom3,       energy,          &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,    &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,    &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy,    &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy

        real(rk),   allocatable, dimension(:)   :: r

        !
        ! Interpolate solution to quadrature nodes
        !
        density       = worker%get_field('Density',    'value', 'face interior')
        mom1          = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2          = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3          = worker%get_field('Momentum-3', 'value', 'face interior')
        energy        = worker%get_field('Energy',     'value', 'face interior')

        grad1_density = worker%get_field('Density',    'grad1', 'face interior', override_lift=.true.)
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', 'face interior', override_lift=.true.)
        grad1_energy  = worker%get_field('Energy',     'grad1', 'face interior', override_lift=.true.)


        grad2_density = worker%get_field('Density',    'grad2', 'face interior',override_lift=.true.)
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', 'face interior',override_lift=.true.)
        grad2_energy  = worker%get_field('Energy',     'grad2', 'face interior',override_lift=.true.)


        grad3_density = worker%get_field('Density',    'grad3', 'face interior',override_lift=.true.)
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', 'face interior',override_lift=.true.)
        grad3_energy  = worker%get_field('Energy',     'grad3', 'face interior',override_lift=.true.)


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if



        !
        ! Compute pressure jacobians
        !
        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
        dp_dmom1    = -(gam-ONE)*mom1/density
        dp_dmom2    = -(gam-ONE)*mom2/density
        dp_dmom3    = -(gam-ONE)*mom3/density
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)


        !
        ! Compute pressure gradient
        !
        grad1_p = dp_ddensity * grad1_density  + &
                  dp_dmom1    * grad1_mom1     + &
                  dp_dmom2    * grad1_mom2     + &
                  dp_dmom3    * grad1_mom3     + &
                  dp_denergy  * grad1_energy

        grad2_p = dp_ddensity * grad2_density  + &
                  dp_dmom1    * grad2_mom1     + &
                  dp_dmom2    * grad2_mom2     + &
                  dp_dmom3    * grad2_mom3     + &
                  dp_denergy  * grad2_energy

        grad3_p = dp_ddensity * grad3_density  + &
                  dp_dmom1    * grad3_mom1     + &
                  dp_dmom2    * grad3_mom2     + &
                  dp_dmom3    * grad3_mom3     + &
                  dp_denergy  * grad3_energy


    end subroutine compute_pressure_gradient
    !******************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_velocity_gradients(worker,  &
                                          grad1_vel1, grad2_vel1, grad3_vel1,   &
                                          grad1_vel2, grad2_vel2, grad3_vel2,   &
                                          grad1_vel3, grad2_vel3, grad3_vel3)
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel3
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel3
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel3

        type(AD_D), allocatable, dimension(:)   ::              &
            density,       mom1,       mom2,       mom3,        &
                           vel1,       vel2,       vel3,        &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3,  &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3,  &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3

        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density       = worker%get_field('Density',    'value', 'face interior')
        mom1          = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2          = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3          = worker%get_field('Momentum-3', 'value', 'face interior')

        grad1_density = worker%get_field('Density',    'grad1', 'face interior', override_lift=.true.)
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', 'face interior', override_lift=.true.)


        grad2_density = worker%get_field('Density',    'grad2', 'face interior',override_lift=.true.)
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', 'face interior',override_lift=.true.)


        grad3_density = worker%get_field('Density',    'grad3', 'face interior',override_lift=.true.)
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', 'face interior',override_lift=.true.)

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if

        
        !
        ! Compute velocities
        !
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density


        !
        ! Compute velocity gradient
        !
        grad1_vel1 = (ONE/density)*grad1_mom1  -  (vel1/density)*grad1_density
        grad2_vel1 = (ONE/density)*grad2_mom1  -  (vel1/density)*grad2_density
        grad3_vel1 = (ONE/density)*grad3_mom1  -  (vel1/density)*grad3_density

        grad1_vel2 = (ONE/density)*grad1_mom2  -  (vel2/density)*grad1_density
        grad2_vel2 = (ONE/density)*grad2_mom2  -  (vel2/density)*grad2_density
        grad3_vel2 = (ONE/density)*grad3_mom2  -  (vel2/density)*grad3_density

        grad1_vel3 = (ONE/density)*grad1_mom3  -  (vel3/density)*grad1_density
        grad2_vel3 = (ONE/density)*grad2_mom3  -  (vel3/density)*grad2_density
        grad3_vel3 = (ONE/density)*grad3_mom3  -  (vel3/density)*grad3_density


    end subroutine compute_velocity_gradients
    !*******************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_density_gradient(worker, grad1_density, grad2_density, grad3_density)
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_density
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_density
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_density


        grad1_density = worker%get_field('Density', 'grad1', 'face interior', override_lift=.true.)
        grad2_density = worker%get_field('Density', 'grad2', 'face interior', override_lift=.true.)
        grad3_density = worker%get_field('Density', 'grad3', 'face interior', override_lift=.true.)


    end subroutine compute_density_gradient
    !*******************************************************************************




end module bc_state_auxiliary_interior
