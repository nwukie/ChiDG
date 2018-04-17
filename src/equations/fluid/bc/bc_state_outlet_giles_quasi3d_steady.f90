module bc_state_outlet_giles_quasi3d_steady
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use bc_state_fluid_averaging,   only: bc_fluid_averaging_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_AllReduce, mpi_comm, &
                                      MPI_INTEGER, MPI_BCast, MPI_MIN, MPI_MAX
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none





    !>  Steady Quasi-3D Giles Nonreflecting Outlet Boundary Condition
    !!
    !!  Options:
    !!      : Average Pressure
    !!      : Pitch
    !!
    !!  References:
    !!              
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, extends(bc_fluid_averaging_t) :: outlet_giles_quasi3d_steady_t

        integer(ik) :: nr = 10
        integer(ik) :: nfourier_space = 8

        real(rk),   allocatable :: r(:)
        real(rk)                :: theta_ref

        type(element_info_t),   allocatable :: donor(:,:)
        real(rk),               allocatable :: donor_node(:,:,:)
        real(rk),               allocatable :: theta(:,:)

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_postcomm     ! Implement specialized initialization
        procedure   :: compute_bc_state     ! boundary condition function implementation

        procedure   :: compute_fourier_decomposition
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

    end type outlet_giles_quasi3d_steady_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_giles_quasi3d_steady_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Giles Quasi3D Steady')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure','Required')
        call self%bcproperties%add('Pitch',           'Required')

    end subroutine init
    !********************************************************************************




    !>  Default specialized initialization procedure. This is called from the base 
    !!  bc%init procedure and can be overwritten by derived types to implement 
    !!  specialized initiailization details.
    !!
    !!  By default, this routine does nothing. However, a particular bc_state_t could 
    !!  reimplement this routine to perform some specialized initialization calculations 
    !!  during initialization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine init_bc_postcomm(self,mesh,group_ID,bc_comm)
        class(outlet_giles_quasi3d_steady_t),   intent(inout)   :: self
        type(mesh_t),                           intent(inout)   :: mesh
        integer(ik),                            intent(in)      :: group_ID
        type(mpi_comm),                         intent(in)      :: bc_comm

        call self%analyze_bc_geometry(mesh,group_ID,bc_comm)

        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm)

    end subroutine init_bc_postcomm
    !*************************************************************************************





    !>  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_giles_quasi3d_steady_t),   intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop
        type(mpi_comm),                                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,                            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,                           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            vel1_bc, vel2_bc, vel3_bc, pressure_bc,                                     &
            vel1_m,  vel2_m,  vel3_m,  pressure_m,                                      &
            c1,    c2,    c3,    c4,    c5,                                             &
            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar,             &
            ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero

        type(AD_D), allocatable, dimension(:,:) ::                                              &
            density_hat_real, vel1_hat_real, vel2_hat_real, vel3_hat_real, pressure_hat_real,   &
            density_hat_imag, vel1_hat_imag, vel2_hat_imag, vel3_hat_imag, pressure_hat_imag,   &
            c1_hat_real,      c2_hat_real,   c3_hat_real,   c4_hat_real,   c5_hat_real,         &
            c1_hat_imag,      c2_hat_imag,   c3_hat_imag,   c4_hat_imag,   c5_hat_imag,         &
            c5_hat_real_gq,   c5_hat_imag_gq


        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg,     &
                       c_avg, ddensity_mean, dvel1_mean, dvel2_mean, dvel3_mean,    &
                       dpressure_mean, density_bar_r, vel1_bar_r, vel2_bar_r,       &
                       vel3_bar_r, pressure_bar_r, c_bar_r, A3_real, A3_imag,       &
                       A4_real, A4_imag, beta

        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset
        type(point_t),  allocatable                 :: coords(:)
        integer                                     :: i, ngq, ivec, imode, iradius, ierr, igq, nmodes


        ! Get back pressure from function.
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',           worker%time(),worker%coords())


        ! Interpolate interior solution to face quadrature nodes
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')


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



        ! Store boundary gradient state. Grad(Q_bc). Do this here, before we
        ! compute any transformations for cylindrical.
        call worker%store_bc_state('Density'   , grad1_density_m, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density_m, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', grad1_mom1_m,    'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1_m,    'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-2', grad1_mom2_m,    'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2_m,    'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-3', grad1_mom3_m,    'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3_m,    'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3_m,    'grad3')
                                                
        call worker%store_bc_state('Energy'    , grad1_energy_m,  'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy_m,  'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy_m,  'grad3')


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)
        end if

        ! Compute velocity and pressure
        vel1_m = mom1_m/density_m
        vel2_m = mom2_m/density_m
        vel3_m = mom3_m/density_m
        pressure_m = worker%get_field('Pressure', 'value', 'face interior')

        ! Update average pressure
        call self%compute_averages(worker,bc_COMM,vel1_avg, vel2_avg, vel3_avg, density_avg, pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)


        ! Compute Fourier decomposition at set of radial stations: 
        !   : U_hat(nmodes,nradius)
        call self%compute_fourier_decomposition(worker,bc_COMM,                        &
                                                density_hat_real,  density_hat_imag,   &
                                                vel1_hat_real,     vel1_hat_imag,      &
                                                vel2_hat_real,     vel2_hat_imag,      &
                                                vel3_hat_real,     vel3_hat_imag,      &
                                                pressure_hat_real, pressure_hat_imag,  &
                                                c1_hat_real,       c1_hat_imag,        &
                                                c2_hat_real,       c2_hat_imag,        &
                                                c3_hat_real,       c3_hat_imag,        &
                                                c4_hat_real,       c4_hat_imag,        &
                                                c5_hat_real,       c5_hat_imag)



        ! Solve for c5 using nonreflecting condition
        nmodes = size(c5_hat_real,1)
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_hat_real( 1,iradius)
            vel1_bar_r     = vel1_hat_real(    1,iradius)
            vel2_bar_r     = vel2_hat_real(    1,iradius)
            vel3_bar_r     = vel3_hat_real(    1,iradius)
            pressure_bar_r = pressure_hat_real(1,iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)

!            ! The imaginary part of beta has already been accounted for in
!            ! the expressions for A2 and A3
!            beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
!            A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!            A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!
!            A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
!            A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)
!
!
!            ! Compute c5 according to nonreflecting condition
!            !
!            !   hat{c5} = A3*hat{c3}  -  A4*hat{c4}
!            !
!            do imode = 2,size(c5_hat_real,1) 
!                c5_hat_real(imode,iradius) = (A3_real*c3_hat_real(imode,iradius) - A3_imag*c3_hat_imag(imode,iradius))  &   ! A3*c3 (real)
!                                           - (A4_real*c4_hat_real(imode,iradius) - A4_imag*c4_hat_imag(imode,iradius))      ! A4*c4 (real)
!                c5_hat_imag(imode,iradius) = (A3_imag*c3_hat_real(imode,iradius) + A3_real*c3_hat_imag(imode,iradius))  &   ! A3*c3 (imag)
!                                           - (A4_imag*c4_hat_real(imode,iradius) + A4_real*c4_hat_imag(imode,iradius))      ! A4*c4 (imag)
!            end do !imode

            do imode = 2,nmodes ! -1 here because the first mode is treated with 1D characteristics
                ! The imaginary part of beta has already been accounted for in
                ! the expressions for A2 and A3
                if (imode <= (nmodes-1)/2 + 1) then
                    beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                else if (imode > (nmodes-1)/2 + 1) then
                    beta = -sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                end if

                A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
                A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)

                A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
                A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)

                c5_hat_real(imode,iradius) = (A3_real*c3_hat_real(imode,iradius) - A3_imag*c3_hat_imag(imode,iradius))  &   ! A3*c3 (real)
                                           - (A4_real*c4_hat_real(imode,iradius) - A4_imag*c4_hat_imag(imode,iradius))      ! A4*c4 (real)
                c5_hat_imag(imode,iradius) = (A3_imag*c3_hat_real(imode,iradius) + A3_real*c3_hat_imag(imode,iradius))  &   ! A3*c3 (imag)
                                           - (A4_imag*c4_hat_real(imode,iradius) + A4_real*c4_hat_imag(imode,iradius))      ! A4*c4 (imag)
            end do !imode


        end do !iradius


        ! Interpolate c5 to the correct radial stations for quadrature nodes
        coords = worker%coords()
        allocate(c5_hat_real_gq(size(c5_hat_real,1),size(coords)), &
                 c5_hat_imag_gq(size(c5_hat_imag,1),size(coords)), stat=ierr)
        if (ierr /= 0) call AllocationError
        c5_hat_real_gq = c5_hat_real(1,1)
        c5_hat_imag_gq = c5_hat_real(1,1)
        c5_hat_real_gq = ZERO   
        c5_hat_imag_gq = ZERO
        do igq = 1,size(coords)
            do imode = 2,size(c5_hat_real,1) ! not interpolating mode1, so it remains zero and isn't present in idft
                c5_hat_real_gq(imode,igq) = interpolate_linear_ad(self%r,c5_hat_real(imode,:),coords(igq)%c1_)
                c5_hat_imag_gq(imode,igq) = interpolate_linear_ad(self%r,c5_hat_imag(imode,:),coords(igq)%c1_)
            end do
        end do


        ! Evaluate c5 at radius to correct theta
        expect_zero = [AD_D(1)]
        c5_3d = c5_hat_real_gq(1,:)
        c5_3d = ZERO
        do igq = 1,size(coords)
            theta_offset = coords(igq)%c2_ - self%theta_ref
            ! We include all modes here for generality, but we already set mode1 to zero
            ! so we are only getting the perturbation part.
            !c5_3d(igq:igq) = idft_eval(c5_hat_real_gq(:,igq),c5_hat_imag_gq(:,igq),[theta_offset]/pitch(1))
            call idft_eval(c5_hat_real_gq(:,igq),   &
                           c5_hat_imag_gq(:,igq),   &
                           [theta_offset]/pitch(1), &
                           c5_3d(igq:igq),          &
                           expect_zero)
        end do


        ! Handle perturbation from local radial mean (m /= 0)
        c1_3d = c5_3d
        c2_3d = c5_3d
        c3_3d = c5_3d
        c4_3d = c5_3d
        c1_3d = ZERO
        c2_3d = ZERO
        c3_3d = ZERO
        c4_3d = ZERO
        density_bar  = c5_3d
        vel1_bar     = c5_3d
        vel2_bar     = c5_3d
        vel3_bar     = c5_3d
        pressure_bar = c5_3d
        density_bar  = ZERO
        vel1_bar     = ZERO
        vel2_bar     = ZERO
        vel3_bar     = ZERO
        pressure_bar = ZERO
        do igq = 1,size(coords)
            density_bar(igq)  = interpolate_linear_ad(self%r, density_hat_real( 1,:), coords(igq)%c1_)
            vel1_bar(igq)     = interpolate_linear_ad(self%r, vel1_hat_real(    1,:), coords(igq)%c1_)
            vel2_bar(igq)     = interpolate_linear_ad(self%r, vel2_hat_real(    1,:), coords(igq)%c1_)
            vel3_bar(igq)     = interpolate_linear_ad(self%r, vel3_hat_real(    1,:), coords(igq)%c1_)
            pressure_bar(igq) = interpolate_linear_ad(self%r, pressure_hat_real(1,:), coords(igq)%c1_)
        end do
        c_bar = sqrt(gam*pressure_bar/density_bar)


        ! Compute perturbation from radius-local mean
        ddensity  = density_m  - density_bar
        dvel1     = vel1_m     - vel1_bar
        dvel2     = vel2_m     - vel2_bar
        dvel3     = vel3_m     - vel3_bar
        dpressure = pressure_m - pressure_bar

        ! Compute characteristics 1-4 associated with local perturbation. 
        ! 5 was already handled from nonreflecting condition on Fourier modes.
        c1_3d = (-c_bar*c_bar)*ddensity    +  (ONE)*dpressure
        c2_3d = (density_bar*c_bar)*dvel1
        c3_3d = (density_bar*c_bar)*dvel2
        c4_3d = (density_bar*c_bar)*dvel3  +  (ONE)*dpressure

        ! Handle m=0 perturbation
        c1_1d = density_m
        c2_1d = density_m
        c3_1d = density_m
        c4_1d = density_m
        c5_1d = density_m
        c1_1d = ZERO
        c2_1d = ZERO
        c3_1d = ZERO
        c4_1d = ZERO
        c5_1d = ZERO

        ! Compute 1-4 characteristics from extrapolation
        ddensity  = density_bar  - density_avg 
        dvel1     = vel1_bar     - vel1_avg
        dvel2     = vel2_bar     - vel2_avg
        dvel3     = vel3_bar     - vel3_avg
        dpressure = pressure_bar - pressure_avg
        do igq = 1,size(ddensity)
            c1_1d(igq) = -c_avg*c_avg*ddensity(igq)    +  dpressure(igq)
            c2_1d(igq) = density_avg*c_avg*dvel1(igq)
            c3_1d(igq) = density_avg*c_avg*dvel2(igq)
            c4_1d(igq) = density_avg*c_avg*dvel3(igq)  +  dpressure(igq)
        end do

        ! Compute characteristic 5 to achieve average pressure
        c5_1d          = -TWO*(pressure_avg - p_user(1))
        ddensity_mean  =  c5_1d(1)/(TWO*c_avg*c_avg)
        dvel3_mean     = -c5_1d(1)/(TWO*density_avg*c_avg)
        dpressure_mean = HALF*c5_1d(1)

        dvel1_mean = dvel3_mean
        dvel2_mean = dvel3_mean
        dvel1_mean = ZERO
        dvel2_mean = ZERO


        ! Contribute average part to boundary state
        density_bc  = density_m
        vel1_bc     = density_m
        vel2_bc     = density_m
        vel3_bc     = density_m
        pressure_bc = density_m


        ! Compose boundary state beginning with average
        density_bc  = density_avg
        vel1_bc     = vel1_avg
        vel2_bc     = vel2_avg
        vel3_bc     = vel3_avg
        pressure_bc = pressure_avg


        ! Contribute perturbation from boundary-global 1D characteristic update
        do igq = 1,size(c1_1d)
            density_bc(igq)  = density_bc(igq)   +  (-ONE/(c_avg*c_avg))*c1_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c4_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c5_1d(igq)
            vel1_bc(igq)     = vel1_bc(igq)      +  (ONE/(density_avg*c_avg))*c2_1d(igq)
            vel2_bc(igq)     = vel2_bc(igq)      +  (ONE/(density_avg*c_avg))*c3_1d(igq)
            vel3_bc(igq)     = vel3_bc(igq)      +  (ONE/(TWO*density_avg*c_avg))*c4_1d(igq)  -  (ONE/(TWO*density_avg*c_avg))*c5_1d(igq)
            pressure_bc(igq) = pressure_bc(igq)  +  HALF*c4_1d(igq)  +  HALF*c5_1d(igq)
        end do


        ! Contribute perturbation from perturbation about radius-local mean from 
        ! quasi-3d nrbc.
        density_bc  = density_bc   +  (-ONE/(c_bar*c_bar))*c1_3d  +  (ONE/(TWO*c_bar*c_bar))*c4_3d  +  (ONE/(TWO*c_bar*c_bar))*c5_3d
        vel1_bc     = vel1_bc      +  (ONE/(density_bar*c_bar))*c2_3d
        vel2_bc     = vel2_bc      +  (ONE/(density_bar*c_bar))*c3_3d
        vel3_bc     = vel3_bc      +  (ONE/(TWO*density_bar*c_bar))*c4_3d  -  (ONE/(TWO*density_bar*c_bar))*c5_3d
        pressure_bc = pressure_bc  +  HALF*c4_3d  +  HALF*c5_3d


        ! Form conserved variables
        density_bc = density_bc
        mom1_bc    = density_bc*vel1_bc
        mom2_bc    = density_bc*vel2_bc
        mom3_bc    = density_bc*vel3_bc
        energy_bc  = pressure_bc/(gam - ONE)  + (density_bc*HALF)*(vel1_bc*vel1_bc + vel2_bc*vel2_bc + vel3_bc*vel3_bc)


        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        ! Store boundary condition state. Q_bc
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


    end subroutine compute_bc_state
    !*********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_fourier_decomposition(self,worker,bc_comm,                   &
                                             density_hat_real,  density_hat_imag,   &
                                             vel1_hat_real,     vel1_hat_imag,      &
                                             vel2_hat_real,     vel2_hat_imag,      &
                                             vel3_hat_real,     vel3_hat_imag,      &
                                             pressure_hat_real, pressure_hat_imag,  &
                                             c1_hat_real,       c1_hat_imag,        &
                                             c2_hat_real,       c2_hat_imag,        &
                                             c3_hat_real,       c3_hat_imag,        &
                                             c4_hat_real,       c4_hat_imag,        &
                                             c5_hat_real,       c5_hat_imag)
        class(outlet_giles_quasi3d_steady_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_imag(:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp,   &
            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
            c1,         c2,         c3,         c4,         c5,                                 &
            ddensity,   dvel1,      dvel2,      dvel3,      dpressure

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik)             :: nmodes, nradius, ntheta, iradius, itheta, ncoeff, imode, ierr
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: pitch(:)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff

        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])

        ! Allocate storage in result
        allocate(density_hat_real( ncoeff,nradius), density_hat_imag( ncoeff,nradius),  &
                 vel1_hat_real(    ncoeff,nradius), vel1_hat_imag(    ncoeff,nradius),  &
                 vel2_hat_real(    ncoeff,nradius), vel2_hat_imag(    ncoeff,nradius),  &
                 vel3_hat_real(    ncoeff,nradius), vel3_hat_imag(    ncoeff,nradius),  &
                 pressure_hat_real(ncoeff,nradius), pressure_hat_imag(ncoeff,nradius),  &
                 c1_hat_real(      ncoeff,nradius), c1_hat_imag(      ncoeff,nradius),  &
                 c2_hat_real(      ncoeff,nradius), c2_hat_imag(      ncoeff,nradius),  &
                 c3_hat_real(      ncoeff,nradius), c3_hat_imag(      ncoeff,nradius),  &
                 c4_hat_real(      ncoeff,nradius), c4_hat_imag(      ncoeff,nradius),  &
                 c5_hat_real(      ncoeff,nradius), c5_hat_imag(      ncoeff,nradius), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius

            ! Interpolate solution to physical_nodes at current radial station
            density = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
            mom1    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
            mom2    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
            mom3    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
            energy  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))

            if (worker%coordinate_system() == 'Cylindrical') then
                mom2 = mom2/self%r(iradius)  ! convert to tangential momentum
            end if

            ! Compute velocities and pressure
            vel1 = mom1/density
            vel2 = mom2/density
            vel3 = mom3/density
            pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )

            ! Compute Fourier transform
            call dft(density,  ZERO*density,  density_real_tmp,  density_imag_tmp )
            call dft(vel1,     ZERO*vel1,     vel1_real_tmp,     vel1_imag_tmp    )
            call dft(vel2,     ZERO*vel2,     vel2_real_tmp,     vel2_imag_tmp    )
            call dft(vel3,     ZERO*vel3,     vel3_real_tmp,     vel3_imag_tmp    )
            call dft(pressure, ZERO*pressure, pressure_real_tmp, pressure_imag_tmp)

            density_hat_real( :,iradius) = density_real_tmp
            density_hat_imag( :,iradius) = density_imag_tmp
            vel1_hat_real(    :,iradius) = vel1_real_tmp
            vel1_hat_imag(    :,iradius) = vel1_imag_tmp
            vel2_hat_real(    :,iradius) = vel2_real_tmp
            vel2_hat_imag(    :,iradius) = vel2_imag_tmp
            vel3_hat_real(    :,iradius) = vel3_real_tmp
            vel3_hat_imag(    :,iradius) = vel3_imag_tmp
            pressure_hat_real(:,iradius) = pressure_real_tmp
            pressure_hat_imag(:,iradius) = pressure_imag_tmp

            ! Get average term
            density_bar  = density_hat_real( 1,iradius)
            vel1_bar     = vel1_hat_real(    1,iradius)
            vel2_bar     = vel2_hat_real(    1,iradius)
            vel3_bar     = vel3_hat_real(    1,iradius)
            pressure_bar = pressure_hat_real(1,iradius)
            c_bar = sqrt(gam*pressure_bar/density_bar)

            ! Compute perturbation
            ddensity  = density  - density_bar
            dvel1     = vel1     - vel1_bar
            dvel2     = vel2     - vel2_bar
            dvel3     = vel3     - vel3_bar
            dpressure = pressure - pressure_bar

            ! Convert perturbation to 1D characteristics
            c1 = ddensity
            c2 = ddensity
            c3 = ddensity
            c4 = ddensity
            c5 = ddensity
            do itheta = 1,size(ddensity)
                c1(itheta) = -(c_bar*c_bar)*ddensity(itheta)    +  (ONE)*dpressure(itheta)
                c2(itheta) = (density_bar*c_bar)*dvel1(itheta)
                c3(itheta) = (density_bar*c_bar)*dvel2(itheta)
                c4(itheta) = (density_bar*c_bar)*dvel3(itheta)  +  (ONE)*dpressure(itheta)
                c5(itheta) = -(density_bar*c_bar)*dvel3(itheta) +  (ONE)*dpressure(itheta)
            end do

            ! Compute Fourier transform of characteristic variables
            call dft(c1, ZERO*c1, c1_real_tmp, c1_imag_tmp)
            call dft(c2, ZERO*c2, c2_real_tmp, c2_imag_tmp)
            call dft(c3, ZERO*c3, c3_real_tmp, c3_imag_tmp)
            call dft(c4, ZERO*c4, c4_real_tmp, c4_imag_tmp)
            call dft(c5, ZERO*c5, c5_real_tmp, c5_imag_tmp)


            ! Adjust Fourier coefficients so their phase is relative to self%theta_ref
            ! instead of the minimum theta of the transform.
            !
            !       q(relative to theta_ref) = q(relative to theta_min) * e^(j 2pi imode delta_theta/pitch)
            !
            ! NOTE: self%theta(:,1) are defined to be the DFT-theta_min at each radius
            !
            do imode = 1,size(c1_hat_real,1)
                shift_r = realpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/pitch(1)))
                shift_i = imagpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/pitch(1)))

                c1_hat_real(imode,iradius) = c1_real_tmp(imode)*shift_r - c1_imag_tmp(imode)*shift_i
                c2_hat_real(imode,iradius) = c2_real_tmp(imode)*shift_r - c2_imag_tmp(imode)*shift_i
                c3_hat_real(imode,iradius) = c3_real_tmp(imode)*shift_r - c3_imag_tmp(imode)*shift_i
                c4_hat_real(imode,iradius) = c4_real_tmp(imode)*shift_r - c4_imag_tmp(imode)*shift_i
                c5_hat_real(imode,iradius) = c5_real_tmp(imode)*shift_r - c5_imag_tmp(imode)*shift_i

                c1_hat_imag(imode,iradius) = c1_imag_tmp(imode)*shift_r + c1_real_tmp(imode)*shift_i
                c2_hat_imag(imode,iradius) = c2_imag_tmp(imode)*shift_r + c2_real_tmp(imode)*shift_i
                c3_hat_imag(imode,iradius) = c3_imag_tmp(imode)*shift_r + c3_real_tmp(imode)*shift_i
                c4_hat_imag(imode,iradius) = c4_imag_tmp(imode)*shift_r + c4_real_tmp(imode)*shift_i
                c5_hat_imag(imode,iradius) = c5_imag_tmp(imode)*shift_r + c5_real_tmp(imode)*shift_i
            end do !imode

        end do !iradius

    end subroutine compute_fourier_decomposition
    !*********************************************************************************





    !>  Determine rmin, rmax, and average theta.
    !!
    !!  Initialize radial stations and reference theta for Fourier transform.
    !!  
    !!  Defines:
    !!      self%r(:)
    !!      self%theta_ref
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine analyze_bc_geometry(self,mesh,group_ID,bc_comm)
        class(outlet_giles_quasi3d_steady_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik) :: bc_idomain_l, bc_ielement_l, patch_ID, face_ID, iface, ierr, inode
        real(rk)    :: face_rmin, face_rmax, local_rmin, local_rmax, global_rmin, global_rmax,  &
                       face_thetamin, face_thetamax, local_thetamin, local_thetamax, global_thetamin, global_thetamax
        real(rk)    :: ref_nodes(8,3), physical_nodes(8,3)

        ! Search for min/max radius on local processor
        local_rmin     =  HUGE(1._rk) ! any radius will be smaller than this, so it is guarunteed to be reset.
        local_rmax     = -HUGE(1._rk) ! any radius will be larger than this, so it is guarunteed to be reset.
        local_thetamin =  HUGE(1._rk) ! any theta will be smaller than this, so it is guarunteed to be reset.
        local_thetamax = -HUGE(1._rk) ! any theta will be larger than this, so it is guarunteed to be reset.
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                iface = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                ! Pick points to evaluate coordinates
                if (iface == XI_MIN) then
                    ref_nodes(1,:) = [-ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [-ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [-ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [-ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [-ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [-ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [-ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [-ONE, ZERO, -ONE]
                else if (iface == XI_MAX) then
                    ref_nodes(1,:) = [ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [ONE, ZERO, -ONE]

                else if (iface == ETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, -ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, -ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, -ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, ONE, -ONE]

                else if (iface == ZETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, -ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, -ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, -ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, -ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, -ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ZETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, -ONE, ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, ONE]

                else
                    call chidg_signal(FATAL,"outlet_giles_quasi3d_steady: analyze_bc_geometry, invalid face indec.")
                end if

                ! Evaluate physical coordinates on face edges
                bc_idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                bc_ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                do inode = 1,size(ref_nodes,1)
                    physical_nodes(inode,:) = mesh%domain(bc_idomain_l)%elems(bc_ielement_l)%physical_point(ref_nodes(inode,:),'Deformed')
                end do !inode

                ! Get face min/max radius
                face_rmin     = minval(physical_nodes(:,1))
                face_rmax     = maxval(physical_nodes(:,1))
                face_thetamin = minval(physical_nodes(:,2))
                face_thetamax = maxval(physical_nodes(:,2))

                ! Update processor-local value if new min/max values were found on the face
                if (face_rmin < local_rmin)         local_rmin     = face_rmin
                if (face_rmax > local_rmax)         local_rmax     = face_rmax
                if (face_thetamin < local_thetamin) local_thetamin = face_thetamin
                if (face_thetamax > local_thetamax) local_thetamax = face_thetamax

            end do !face_ID
        end do !patch_ID

        ! Reduce processor local values to determine boundary-global min/max values
        call MPI_AllReduce(local_rmin,    global_rmin,    1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_rmax,    global_rmax,    1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        call MPI_AllReduce(local_thetamin,global_thetamin,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_thetamax,global_thetamax,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        
        ! Create radial stations
        self%r = linspace(global_rmin,global_rmax,self%nr)

        ! Compute theta_ref
        self%theta_ref = (global_thetamin + global_thetamax)/TWO

    end subroutine analyze_bc_geometry
    !********************************************************************************




    !>  For each radial station, a theta-discretization exists upon which the 
    !!  discrete Fourier transform operation is computed. Since the theta-discretization
    !!  spans the entire boundary, a general interpolation procedure is used to fill
    !!  solution values that could come from many different elements on the boundary.
    !!
    !!  This procedure finds donor elements for each node in the theta-discretization
    !!  at each radial station. Additionally, the location of the physical coordinate
    !!  within the donor element's reference coordinate system is also found.
    !!
    !!  This procedure
    !!
    !!  Initializes:
    !!  ------------------------------
    !!  self%theta(:)
    !!  self%donor(:,:)
    !!  self%donor_node(:,:,:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine initialize_fourier_discretization(self,mesh,group_ID,bc_comm)
        class(outlet_giles_quasi3d_steady_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik)                 :: nmodes, ncoeff, nradius, ntheta, idomain_l, ielement_l, iface, &
                                       iradius, itheta, ierr, noverset
        real(rk)                    :: dtheta, dtheta_n, midpoint(3), try_offset(3), node(3), z
        real(rk),       allocatable :: pitch(:)
        character(:),   allocatable :: user_msg
        logical                     :: donor_found


        ! Determine z-location of some face on the boundary and assume 
        ! entire boundary is constant-z
        idomain_l  = mesh%bc_patch_group(group_ID)%patch(1)%idomain_l()
        ielement_l = mesh%bc_patch_group(group_ID)%patch(1)%ielement_l(1)
        iface      = mesh%bc_patch_group(group_ID)%patch(1)%iface(1)
        if (iface == XI_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([-ONE,ZERO,ZERO],'Deformed')
        else if (iface == XI_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ ONE,ZERO,ZERO],'Deformed')
        else if (iface == ETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,-ONE,ZERO],'Deformed')
        else if (iface == ETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO, ONE,ZERO],'Deformed')
        else if (iface == ZETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO,-ONE],'Deformed')
        else if (iface == ZETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO, ONE],'Deformed')
        end if
        z = midpoint(3)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff
        
        ! Initialize theta discretization parameters
        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta

        ! Construct theta discretization at each radius
        allocate(self%theta(size(self%r),ntheta), stat=ierr)
        if (ierr /= 0) call AllocationError
        do itheta = 1,ntheta
            self%theta(:,itheta) = self%theta_ref + (itheta-1)*dtheta_n
        end do

        ! Donor search offset, if needed
        try_offset = [ZERO, -pitch(1), ZERO]

        ! For each radial station, initialized donor for each node in theta grid
        allocate(self%donor(nradius,ntheta), self%donor_node(nradius,ntheta,3), stat=ierr)
        if (ierr /= 0) call AllocationError
        do iradius = 1,size(self%r)
            noverset = 0
            do itheta = 1,ntheta

                node = [self%r(iradius), self%theta(iradius,itheta), z]

                ! Try processor-LOCAL elements
                call find_gq_donor(mesh,                                    &
                                   node,                                    &
                                   [ZERO,ZERO,ZERO],                        &
                                   face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                   self%donor(iradius,itheta),              &
                                   self%donor_node(iradius,itheta,1:3),    &
                                   donor_found)

                ! Try LOCAL elements with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                    &
                                       node,                                    &
                                       try_offset,                              &
                                       face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                       self%donor(iradius,itheta),              &
                                       self%donor_node(iradius,itheta,1:3),    &
                                       donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if

                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                [ZERO,ZERO,ZERO],                       &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),   &
                                                donor_found)
                end if

                
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                try_offset,                             &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),   &
                                                donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_state_outlet_giles_quasi3d_steady%initialize_fourier_discretization: &
                            no donor element found for Fourier discretization node."
                if (.not. donor_found) call chidg_signal(FATAL,user_msg)

            end do !itheta

            ! Shift arrays so that we start with the theta_min point
            self%donor(iradius,:)         = cshift(self%donor(iradius,:), -noverset, dim=1)
            self%donor_node(iradius,:,:) = cshift(self%donor_node(iradius,:,:), -noverset, dim=1)
            self%theta(iradius,:)         = cshift(self%theta(iradius,:), -noverset, dim=1)

        end do !iradius

    end subroutine initialize_fourier_discretization
    !**************************************************************************************




end module bc_state_outlet_giles_quasi3d_steady
