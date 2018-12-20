module bc_state_rstm_ssglrrw_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_rstm_lrr_lnw
    implicit none
    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: rstm_ssglrrw_wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type rstm_ssglrrw_wall_t
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
        class(rstm_ssglrrw_wall_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Wall')
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
        class(rstm_ssglrrw_wall_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            distance_elem, density_m, mu_l, &
            density_reynolds_m, grad1_density_reynolds_m, grad2_density_reynolds_m, grad3_density_reynolds_m, & 
            density_nutilde_m, grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m

        real(rk) :: h(3), distance
        integer(ik) :: ii, nnodes, order

        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_nutilde_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_nutilde_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_field('Density * Omega', 'grad3', 'face interior')


!        grad1_density_nutilde_m = ZERO
!        grad2_density_nutilde_m = ZERO
!        grad3_density_nutilde_m = ZERO

        distance_elem       = worker%get_field('Wall Distance', 'value', 'element')

        h = worker%element_size('interior')
        order = worker%solution_order('interior')
        !distance = minval(abs(h))/(real(order, rk)+1.0_rk)
        distance = minval(abs(h))
        !distance = (minval(distance_elem%x_ad_))
        if (ieee_is_nan(distance)) print *, 'rstm wall bc distance is nan'

        !
        ! Store boundary condition state - Dirichlet Zero
        !
        density_m = worker%get_field('Density',     'value', 'face interior')
        if (any(ieee_is_nan(density_m(:)%x_ad_))) print *, 'rstm wall bc density is nan, source : ', worker%interpolation_source
        !if (any(ieee_is_nan(density_m(:)%x_ad_))) stop

        mu_l = worker%get_field('Laminar Viscosity', 'value', 'face interior')
        if (any(ieee_is_nan(mu_l(:)%x_ad_))) print *, 'rstm wall bc mu_l is nan'
        density_nutilde_m = density_m
        !density_nutilde_m = density_m*log(10.0_rk*6.0_rk*mu_l/(density_m*SSG_LRRW_beta_w*distance**TWO+1.0e-16_rk))
        !density_nutilde_m = density_m*log(1.0_rk*6.0_rk*mu_l/(density_m*lrr_lnw_beta_0*distance**TWO+1.0e-16_rk))
        !density_nutilde_m = 10.0_rk*6.0_rk*1.634303551710997e-05/(density_m*lrr_lnw_beta_0*distance**TWO+1.0e-16_rk)
        density_nutilde_m = 8716.72515337_rk
        density_nutilde_m = log(density_nutilde_m)
        density_nutilde_m = density_m*density_nutilde_m
        !density_nutilde_m = density_m*(log(10.0_rk*6.0_rk*1.634303551710997e-05)-log((abs(density_m)*lrr_lnw_beta_0*distance**TWO+1.0e-6_rk)))
        nnodes = size(density_nutilde_m)
        do ii = 1, nnodes
            if (distance < 1.0e-16_rk) print *, 'density omega bc tiny distance'
            if (ieee_is_nan(density_nutilde_m(ii)%x_ad_)) print *, 'density omega wall bc is nan, distance ', distance
            if (any(ieee_is_nan(density_nutilde_m(ii)%xp_ad_))) print *, 'density omega wall bc deriv is nan, density = ', density_m(ii)%x_ad_
            if (any(ieee_is_nan(density_nutilde_m(ii)%xp_ad_))) print *, 'density omega wall bc deriv is nan, densityomega_bc = ', density_nutilde_m(ii)%x_ad_
            if (any(ieee_is_nan(density_nutilde_m(ii)%xp_ad_))) print *, 'density omega wall bc deriv is nan, distance = ', distance
            if (any(ieee_is_nan(density_m(ii)%xp_ad_))) print *, 'density omega wall bc - density deriv is nan'
            if (any(ieee_is_nan(density_m(ii)%xp_ad_))) print *, 'omega wall bc - density deriv is nan'
            if (any(ieee_is_nan(density_nutilde_m(ii)%xp_ad_))) density_nutilde_m(ii)%xp_ad_ = ZERO
        end do

        call worker%store_bc_state('Density * Omega', density_nutilde_m, 'value')


        !
        ! Store boundary condition gradient - Extrapolate
        !
        call worker%store_bc_state('Density * Omega', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_nutilde_m, 'grad3')
                                                

        ! R_11
        density_reynolds_m       = worker%get_field('Density * Reynolds-11', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-11', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-11', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-11', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-11', grad3_density_reynolds_m, 'grad3')

        ! R_22
        density_reynolds_m       = worker%get_field('Density * Reynolds-22', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-22', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-22', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-22', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-22', grad3_density_reynolds_m, 'grad3')

        ! R_33
        density_reynolds_m       = worker%get_field('Density * Reynolds-33', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-33', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-33', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-33', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-33', grad3_density_reynolds_m, 'grad3')

        ! R_12
        density_reynolds_m       = worker%get_field('Density * Reynolds-12', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-12', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-12', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-12', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-12', grad3_density_reynolds_m, 'grad3')

        ! R_13
        density_reynolds_m       = worker%get_field('Density * Reynolds-13', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-13', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-13', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-13', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-13', grad3_density_reynolds_m, 'grad3')

        ! R_23
        density_reynolds_m       = worker%get_field('Density * Reynolds-23', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad3', 'face interior')

        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-23', density_reynolds_m, 'value')
        call worker%store_bc_state('Density * Reynolds-23', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-23', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-23', grad3_density_reynolds_m, 'grad3')


    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_rstm_ssglrrw_wall
