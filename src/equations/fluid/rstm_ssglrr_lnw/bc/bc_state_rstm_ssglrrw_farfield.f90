module bc_state_rstm_ssglrrw_farfield
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, THREE, HALF, ZERO, RKTOL
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use ieee_arithmetic,        only: ieee_is_nan
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use mod_rstm_ssglrrw
    implicit none
    


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: rstm_ssglrrw_farfield_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type rstm_ssglrrw_farfield_t
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
        class(rstm_ssglrrw_farfield_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Farfield')
        call self%set_family('Farfield')


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
        class(rstm_ssglrrw_farfield_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_m, mom1_m, mom2_m, mom3_m, normal_momentum, &
            density_bc, mom1_bc, mom2_bc, mom3_bc,  &
            u_infty_sq, k_t_ff, mu_t_ff,                     &
            density_reynolds_ff, density_reynolds_m, density_reynolds_bc, &
            grad1_density_reynolds_m, grad2_density_reynolds_m, grad3_density_reynolds_m, &
            density_omega_ff, density_omega_m, density_omega_bc,              &
            grad1_density_omega_m, grad2_density_omega_m, grad3_density_omega_m, mu_m

        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3
        logical,    allocatable, dimension(:)   :: inflow, outflow

        integer(ik) :: ii


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_field('Density',    'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
        end if


        
        !
        ! Get unit normal vector
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)


        !
        ! Determine which modes are inflow/outflow
        !
        normal_momentum = mom1_m*unorm_1 + mom2_m*unorm_2 + mom3_m*unorm_3

        !allocate(inflow(size(unorm_1)))
        !allocate(outflow(size(unorm_1)))
        inflow  = ( normal_momentum <= RKTOL )
        outflow = ( normal_momentum >  RKTOL )




        !
        ! Set boundary values for density * omega.
        !   1: Extrapolate all values (assuming outflow)
        !   2: Where inflow is detected, set from user specified parameter
        !
        
        density_bc = worker%get_field('Density',    'value', 'boundary')
        mom1_bc    = worker%get_field('Momentum-1', 'value', 'boundary')
        mom2_bc    = worker%get_field('Momentum-2', 'value', 'boundary')
        mom3_bc    = worker%get_field('Momentum-3', 'value', 'boundary')

        u_infty_sq = (mom1_bc**TWO + mom2_bc**TWO + mom3_bc**TWO)/density_bc**TWO

        ! Get farfield turbulent kinetic energy 
        k_t_ff  = (THREE/TWO)*(SSG_LRRW_tu_infty)**TWO*u_infty_sq
        k_t_ff = rstm_ssglrrw_k_infty

        density_reynolds_ff = (TWO/THREE)*density_bc*k_t_ff
        !density_reynolds_ff = 6.0e-9_rk*347.27942639897344_rk**TWO

        ! Get farfield turbulent viscosity
        mu_m    = worker%get_field('Laminar Viscosity', 'value', 'face interior')
        mu_t_ff = SSG_LRRW_f_mu_infty*mu_m


        density_omega_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_omega_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_omega_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_omega_m = worker%get_field('Density * Omega', 'grad3', 'face interior')



        ! Omega
        ! Compute FF values
        density_omega_ff = density_omega_m
        !density_omega_ff = density_bc*log(density_bc*k_t_ff/mu_t_ff)
        density_omega_ff = density_bc*log(rstm_ssglrrw_omega_infty)
        !density_omega_ff = 1.0e-6_rk*1.1765047303964247_rk*347.27942639897344_rk**TWO/1.634303551710997e-05_rk
        density_omega_bc = density_omega_m
        density_omega_bc = density_omega_ff
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_omega_bc(ii) = density_omega_ff(ii)
            end if

        end do

    !    where (inflow)
    !        density_omega_bc = density_omega_ff
    !    end where
    !

        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Omega', density_omega_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_omega_m = ZERO
        grad2_density_omega_m = ZERO
        grad3_density_omega_m = ZERO
        call worker%store_bc_state('Density * Omega', grad1_density_omega_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_omega_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_omega_m, 'grad3')


        ! R_11
        ! Compute FF values
        density_reynolds_m  = worker%get_field('Density * Reynolds-11', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad3', 'face interior')


        density_reynolds_bc = density_reynolds_m
        density_reynolds_bc = density_reynolds_ff
    !    where(inflow)
    !        density_reynolds_bc = density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-11', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-11', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-11', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-11', grad3_density_reynolds_m, 'grad3')

        ! R_22
        ! Compute FF values
        if (allocated(density_reynolds_m)) deallocate(density_reynolds_m)
        density_reynolds_m  = worker%get_field('Density * Reynolds-22', 'value', 'face interior')
        density_reynolds_bc = density_reynolds_m
        density_reynolds_bc = density_reynolds_ff
    !    where(inflow)
    !        density_reynolds_bc = density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-22', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-22', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-22', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-22', grad3_density_reynolds_m, 'grad3')

         ! R_33
        ! Compute FF values
        if (allocated(density_reynolds_m)) deallocate(density_reynolds_m)
        density_reynolds_m  = worker%get_field('Density * Reynolds-33', 'value', 'face interior')
        density_reynolds_bc = density_reynolds_m
        density_reynolds_bc = density_reynolds_ff
    !    where(inflow)
    !        density_reynolds_bc = density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-33', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-33', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-33', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-33', grad3_density_reynolds_m, 'grad3')

        ! R_12
        ! Compute FF values
        if (allocated(density_reynolds_m)) deallocate(density_reynolds_m)
        density_reynolds_m  = worker%get_field('Density * Reynolds-12', 'value', 'face interior')
        density_reynolds_bc = ZERO*density_reynolds_m
    !    where(inflow)
    !        density_reynolds_bc = ZERO*density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = ZERO*density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-12', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-12', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-12', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-12', grad3_density_reynolds_m, 'grad3')

         ! R_13
        ! Compute FF values
        if (allocated(density_reynolds_m)) deallocate(density_reynolds_m)
        density_reynolds_m  = worker%get_field('Density * Reynolds-13', 'value', 'face interior')
        density_reynolds_bc = ZERO*density_reynolds_m
    !    where(inflow)
    !        density_reynolds_bc = ZERO*density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = ZERO*density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-13', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-13', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-13', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-13', grad3_density_reynolds_m, 'grad3')
                                    
        ! R_23
        ! Compute FF values
        if (allocated(density_reynolds_m)) deallocate(density_reynolds_m)
        density_reynolds_m  = worker%get_field('Density * Reynolds-23', 'value', 'face interior')
        density_reynolds_bc = ZERO*density_reynolds_m
    !    where(inflow)
    !        density_reynolds_bc = ZERO*density_reynolds_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_reynolds_bc(ii) = ZERO*density_reynolds_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * Reynolds-23', density_reynolds_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-23', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-23', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-23', grad3_density_reynolds_m, 'grad3')
     
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_rstm_ssglrrw_farfield
