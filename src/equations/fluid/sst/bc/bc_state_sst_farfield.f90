module bc_state_sst_farfield
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, THREE, HALF, ZERO, RKTOL
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use ieee_arithmetic,        only: ieee_is_nan
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use mod_sst
    implicit none
    


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: sst_farfield_t
        real(rk)    :: k_infty      = 1.0e-3_rk
        real(rk)    :: omega_infty  = 345.0_rk

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type sst_farfield_t
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
        class(sst_farfield_t),   intent(inout) :: self
        real(rk)            :: k_infty, omega_infty
        integer             :: unit, msg
        logical             :: file_exists

        namelist /k_omega/   k_infty, omega_infty 


        !
        ! Set operator name
        !
        call self%set_name('SST Farfield')
        call self%set_family('Farfield')
        
        
        call self%bcproperties%add('Freestream Density', 'Required')
        call self%bcproperties%add('Freestream Flow Speed', 'Required')
        call self%set_fcn_option('Freestream Density', 'val', 1.17_rk)
        call self%set_fcn_option('Freestream Flow Speed', 'val', 100._rk)
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%k_infty/omega_infty
        !   2: if not available, do nothing and retain default values
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=k_omega,iostat=msg)
            if (msg == 0) self%k_infty = k_infty 
            if (msg == 0) self%omega_infty = omega_infty 
            close(unit)
        end if


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
        class(sst_farfield_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_m, mom1_m, mom2_m, mom3_m, normal_momentum, &
            density_bc, mom1_bc, mom2_bc, mom3_bc,  &
            u_infty_sq,                      &
            density_k_ff, density_k_m, density_k_bc, &
            grad1_density_k_m, grad2_density_k_m, grad3_density_k_m, &
            density_omega_ff, density_omega_m, density_omega_bc,              &
            grad1_density_omega_m, grad2_density_omega_m, grad3_density_omega_m, mu_m

        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3, u_infty, density, k_t_ff, mu_t_ff
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

        u_infty = self%bcproperties%compute('Freestream Flow Speed', worker%time(), worker%coords() )
        density = self%bcproperties%compute('Freestream Density', worker%time(), worker%coords() )
        
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
        
        
        !density_bc = worker%get_field('Density',    'value', 'boundary')
        
        !mom1_bc    = worker%get_field('Momentum-1', 'value', 'boundary')
        
        !mom2_bc    = worker%get_field('Momentum-2', 'value', 'boundary')
        
        !mom3_bc    = worker%get_field('Momentum-3', 'value', 'boundary')

        !u_infty_sq = (mom1_bc**TWO + mom2_bc**TWO + mom3_bc**TWO)/density_bc**TWO

        ! Get farfield turbulent kinetic energy 
        !k_t_ff  = (THREE/TWO)*(sst_tu_infty)**TWO*u_infty**TWO

        density_k_ff = density_m
        !density_k_ff = density*k_t_ff
        !density_k_ff = density*sst_k_infty
        density_k_ff = density*self%k_infty
        !print *, 'density k ff'
        !print *, density_k_ff(:)%x_ad_
        !density_k_ff = 6.0e-9_rk*347.27942639897344_rk**TWO

        ! Get farfield turbulent viscosity
        mu_m    = worker%get_field('Laminar Viscosity', 'value', 'face interior')
        mu_t_ff = sst_f_mu_infty*mu_m(:)%x_ad_


        density_omega_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_omega_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_omega_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_omega_m = worker%get_field('Density * Omega', 'grad3', 'face interior')



        ! Omega
        ! Compute FF values
        density_omega_ff = density_omega_m
        !density_omega_ff = density*log(density*k_t_ff/mu_t_ff)
        !density_omega_ff = density*log(sst_omega_infty)
        density_omega_ff = density*log(self%omega_infty)
        density_omega_bc = density_omega_m
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
        density_k_m  = worker%get_field('Density * k', 'value', 'face interior')
        grad1_density_k_m = worker%get_field('Density * k', 'grad1', 'face interior')
        grad2_density_k_m = worker%get_field('Density * k', 'grad2', 'face interior')
        grad3_density_k_m = worker%get_field('Density * k', 'grad3', 'face interior')


        density_k_bc = density_k_m
        !print *, 'density k m'
        !print *, density_k_bc(:)%x_ad_
        !density_k_bc = density_k_ff
    !    where(inflow)
    !        density_k_bc = density_k_ff
    !    end where
        do ii = 1, size(inflow) 
            if (inflow(ii)) then
                density_k_bc(ii) = density_k_ff(ii)
            end if

        end do



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * k', density_k_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_k_m = ZERO
        grad2_density_k_m = ZERO
        grad3_density_k_m = ZERO
        call worker%store_bc_state('Density * k', grad1_density_k_m, 'grad1')
        call worker%store_bc_state('Density * k', grad2_density_k_m, 'grad2')
        call worker%store_bc_state('Density * k', grad3_density_k_m, 'grad3')

        
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_sst_farfield
