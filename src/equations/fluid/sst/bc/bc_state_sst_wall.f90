module bc_state_sst_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_sst
    implicit none
    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: sst_wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type sst_wall_t
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
        class(sst_wall_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('SST Wall')
        call self%set_family('Wall')

    end subroutine init
    !********************************************************************************







    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(sst_wall_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            distance_elem, density_m, mu_l, tau_wall, shear_11, shear_22, shear_33,&
            shear_12, shear_13, shear_23,&
            r_1, r_2, r_3, density_omega_lb, &
            density_k_m, grad1_density_k_m, grad2_density_k_m, grad3_density_k_m, k_m, mu_t,& 
            density_omega_m, grad1_density_omega_m, grad2_density_omega_m, grad3_density_omega_m, y_plus_c

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3, div_vel,                  &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            du_ddensity,   dv_ddensity,   dw_ddensity,      &
            du_dmom1,      dv_dmom2,      dw_dmom3,         &
            invdensity


        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3, distance_field
        real(rk),   allocatable, dimension(:,:) :: h_sm
        real(rk) :: h(3), distance, alpha_p, y_plus
        integer(ik) :: ii, nnodes, order


        !
        ! Get omega wall bc correction factor, from Schoenwa and Hartmann (2014)
        !

        order = worker%solution_order('interior')
        alpha_p = 1.0_rk/sqrt(10.0_rk)
        if (order==0) then
            alpha_p = 3.6788e-1_rk 
        else if (order==1) then
            alpha_p = 8.2085e-2_rk 
        else if (order==2) then
            alpha_p = 3.5674e-2_rk 
        else if (order==3) then
            alpha_p = 1.9908e-2_rk 
        else if (order==4) then
            alpha_p = 1.2694e-2_rk 
        else if (order==5) then
            alpha_p = 8.7973e-3_rk 
        else if (order==6) then
            alpha_p = 6.4555e-3_rk 
        else if (order==7) then
            alpha_p = 4.9386e-3_rk 
        else if (order==8) then
            alpha_p = 3.9000e-3_rk 
        else if (order==9) then
            alpha_p = 3.1578e-3_rk 
        else if (order==10) then
            alpha_p = 2.6091e-3_rk
        end if


        !
        ! Construct wall shear stress
        !
        ! NOTE: We are not able to use the shear stress model because of the structure of the 
        !       cache handler. It computes BC before computing models with f(Grad(Q)) dependency.
        !       Because lifting operators need the BC state to compute a jump, this is not easily resolved.
        !
        density_m = worker%get_field('Density',     'value', 'face interior')
        density = worker%get_field('Density',    'value', 'face interior')
        mom1    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3    = worker%get_field('Momentum-3', 'value', 'face interior')

        k_m = worker%get_field('k',     'value', 'face interior')
        density_omega_m = worker%get_field('Density * Omega',     'value', 'face interior')
        mu_t = density_m*k_m*exp(-density_omega_m/density_m)

        grad1_density    = worker%get_field('Density'   , 'grad1',  'face interior', override_lift=.true.)
        grad2_density    = worker%get_field('Density'   , 'grad2',  'face interior', override_lift=.true.)
        grad3_density    = worker%get_field('Density'   , 'grad3',  'face interior', override_lift=.true.)

        grad1_mom1       = worker%get_field('Momentum-1', 'grad1',  'face interior', override_lift=.true.)
        grad2_mom1       = worker%get_field('Momentum-1', 'grad2',  'face interior', override_lift=.true.)
        grad3_mom1       = worker%get_field('Momentum-1', 'grad3',  'face interior', override_lift=.true.)

        grad1_mom2       = worker%get_field('Momentum-2', 'grad1',  'face interior', override_lift=.true.)
        grad2_mom2       = worker%get_field('Momentum-2', 'grad2',  'face interior', override_lift=.true.)
        grad3_mom2       = worker%get_field('Momentum-2', 'grad3',  'face interior', override_lift=.true.)

        grad1_mom3       = worker%get_field('Momentum-3', 'grad1',  'face interior', override_lift=.true.)
        grad2_mom3       = worker%get_field('Momentum-3', 'grad2',  'face interior', override_lift=.true.)
        grad3_mom3       = worker%get_field('Momentum-3', 'grad3',  'face interior', override_lift=.true.)

        invdensity = ONE/(density )


        
        !
        ! compute velocity jacobians
        !
        du_ddensity  = -invdensity*invdensity*mom1
        dv_ddensity  = -invdensity*invdensity*mom2
        dw_ddensity  = -invdensity*invdensity*mom3

        du_dmom1 = invdensity
        dv_dmom2 = invdensity
        dw_dmom3 = invdensity



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
        grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1
        grad3_u = du_ddensity*grad3_density  +  du_dmom1*grad3_mom1

        grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
        grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2
        grad3_v = dv_ddensity*grad3_density  +  dv_dmom2*grad3_mom2

        grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
        grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
        grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3

      

        mu_l = worker%get_field('Laminar Viscosity', 'value', 'face interior')
        div_vel = grad1_u + grad2_v + grad3_w
        shear_11 = (mu_l+mu_t)*(grad1_u + grad1_u - (TWO/THREE)*div_vel)
        shear_22 = (mu_l+mu_t)*(grad2_v + grad2_v - (TWO/THREE)*div_vel)
        shear_33 = (mu_l+mu_t)*(grad3_w + grad3_w - (TWO/THREE)*div_vel)
        shear_12 = (mu_l+mu_t)*(grad1_v + grad2_u)
        shear_13 = (mu_l+mu_t)*(grad1_w + grad3_u)
        shear_23 = (mu_l+mu_t)*(grad2_w + grad3_v)

        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)

        r_1 = shear_11*unorm_1 + shear_12*unorm_2 + shear_13*unorm_3
        r_2 = shear_12*unorm_1 + shear_22*unorm_2 + shear_23*unorm_3
        r_3 = shear_13*unorm_1 + shear_23*unorm_2 + shear_33*unorm_3
        
        tau_wall = shear_11

        tau_wall = sqrt(r_1**TWO+ r_2**TWO + r_3**TWO)



        
        !distance_elem       = worker%get_field('Wall Distance', 'value', 'element')

        h = worker%element_size('interior')
        distance = minval(abs(h))


        h_sm = worker%h_smooth('boundary')
        distance_field =  h_sm(:,1)*unorm_1 + h_sm(:,2)*unorm_2 + h_sm(:,3)*unorm_3
        !
        ! Store boundary condition state - Dirichlet Zero
        !


        y_plus_c = mu_l
        
        y_plus_c = distance*sqrt(tau_wall*density_m)/mu_l
        !if (any(1.0_rk<y_plus_c(:)%x_ad_)) print *, 'y_plus computed:', y_plus_c(:)%x_ad_
        ! NOTE: From my reading of Schoenawa and Hartmann (2014), they set y_plus to a specified constant value of 2.5
        !       instead of computing it from the mean flow variables, which instead are used to compute hte wall stress tau_wall.
        y_plus = 2.5_rk
        density_omega_m = density_m


        density_omega_m = density_m*log(6.0_rk*tau_wall/(mu_l*sst_w_beta_k*(alpha_p*y_plus)**TWO+1.0e-15_rk))

        ! 
        ! Set a lower bound
        !
        density_omega_lb = density_m
        density_omega_lb = density_m*log(6.0_rk*mu_l/(density_m*sst_w_beta_k*(alpha_p*distance_field)**TWO+1.0e-15_rk))

        !density_omega_m = density_omega_lb
        nnodes = size(density_omega_m)
        do ii = 1, nnodes
            if (density_omega_m(ii)%x_ad_ < density_omega_lb(ii)%x_ad_) density_omega_m(ii) = density_omega_lb(ii)
            !print *, density_omega_lb(ii)%x_ad_
        end do

        call worker%store_bc_state('Density * Omega', density_omega_m, 'value')

        !
        ! Store boundary condition gradient - Extrapolate
        !
        !
        ! Interpolate interior solution to quadrature nodes
        !
        grad1_density_omega_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_omega_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_omega_m = worker%get_field('Density * Omega', 'grad3', 'face interior')


        call worker%store_bc_state('Density * Omega', grad1_density_omega_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_omega_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_omega_m, 'grad3')
                                                

        ! k 
        density_k_m     = worker%get_field('Density * k',     'value', 'face interior')
        grad1_density_k_m = worker%get_field('Density * k', 'grad1', 'face interior')
        grad2_density_k_m = worker%get_field('Density * k', 'grad2', 'face interior')
        grad3_density_k_m = worker%get_field('Density * k', 'grad3', 'face interior')

        density_k_m = ZERO
        call worker%store_bc_state('Density * k', density_k_m, 'value')
        call worker%store_bc_state('Density * k', grad1_density_k_m, 'grad1')
        call worker%store_bc_state('Density * k', grad2_density_k_m, 'grad2')
        call worker%store_bc_state('Density * k', grad3_density_k_m, 'grad3')


    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_sst_wall
