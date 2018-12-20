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
            distance_elem, density_m, grad2_density, mom1, mu_l, tau_wall, shear_11, shear_22, shear_33,&
            shear_12, shear_13, shear_23,&
            r_1, r_2, r_3, density_omega_lb, &
            density_reynolds_m, grad1_density_reynolds_m, grad2_density_reynolds_m, grad3_density_reynolds_m, & 
            density_nutilde_m, grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m

        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3
        real(rk) :: h(3), distance, alpha_p, y_plus
        integer(ik) :: ii, nnodes, order

        order = worker%solution_order('interior')
        density_m = worker%get_field('Density',     'value', 'face interior')

        !
        ! Get omega wall bc correction factor, from Schoenwa and Hartmann (2014)
        !
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


!        shear_11 = density_m
!        !shear_11 = worker%get_field('Shear-11', 'value', 'face interior')
!        shear_22 = worker%get_field('Shear-22', 'value', 'face interior')
!        shear_33 = worker%get_field('Shear-33', 'value', 'face interior')
!        shear_12 = worker%get_field('Shear-12', 'value', 'face interior')
!        shear_13 = worker%get_field('Shear-13', 'value', 'face interior')
!        shear_23 = worker%get_field('Shear-23', 'value', 'face interior')
!
!        unorm_1 = worker%unit_normal(1)
!        unorm_2 = worker%unit_normal(2)
!        unorm_3 = worker%unit_normal(3)
!
!        r_1 = shear_11*unorm_1 + shear_12*unorm_2 + shear_13*unorm_3
!        r_2 = shear_12*unorm_1 + shear_22*unorm_2 + shear_23*unorm_3
!        r_3 = shear_13*unorm_1 + shear_23*unorm_2 + shear_33*unorm_3
!        
!        tau_wall = shear_11
!
!        tau_wall = sqrt(r_1**TWO+ r_2**TWO + r_3**TWO)
!


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_nutilde_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_nutilde_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_field('Density * Omega', 'grad3', 'face interior')


        !distance_elem       = worker%get_field('Wall Distance', 'value', 'element')

        h = worker%element_size('interior')
        distance = minval(abs(h))

        !
        ! Store boundary condition state - Dirichlet Zero
        !


        mu_l = worker%get_field('Laminar Viscosity', 'value', 'face interior')
        
        ! NOTE: From my reading of Schoenawa and Hartmann (2014), they set y_plus to a specified constant value of 2.5
        !       instead of computing it from the mean flow variables, which instead are used to compute hte wall stress tau_wall.
        y_plus = 2.5_rk
        density_nutilde_m = density_m


        !density_nutilde_m = density_m*log(6.0_rk*tau_wall/(mu_l*sst_w_beta_k*(alpha_p*y_plus)**TWO+1.0e-13_rk))

        ! 
        ! Set a lower bound
        !
        density_omega_lb = density_m
        density_omega_lb = density_m*log(6.0_rk*mu_l/(density_m*sst_w_beta_k*(alpha_p*distance)**TWO+1.0e-12_rk))
        density_nutilde_m = density_omega_lb
        nnodes = size(density_nutilde_m)
!        do ii = 1, nnodes
!            !if (density_nutilde_m(ii)%x_ad_ < density_omega_lb(ii)%x_ad_) density_nutilde_m(ii) = density_omega_lb(ii)
!            !print *, density_omega_lb(ii)%x_ad_
!        end do
!
        call worker%store_bc_state('Density * Omega', density_nutilde_m, 'value')

        !
        ! Store boundary condition gradient - Extrapolate
        !
        call worker%store_bc_state('Density * Omega', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_nutilde_m, 'grad3')
                                                

        ! k 
        density_reynolds_m       = worker%get_field('Density * k', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * k', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * k', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * k', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * k', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * k', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * k', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * k', grad3_density_reynolds_m, 'grad3')


    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_sst_wall
