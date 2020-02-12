module bc_state_inlet_characteristic
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ONE, TWO, ZERO, HALF
    use mod_fluid,                  only: Rgas, cp, gam
    use mod_dft,                    only: dft
    use type_bc_state,              only: bc_state_t
    use bc_state_fluid_averaging,   only: bc_fluid_averaging_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use type_point,                 only: point_t
    use mpi_f08,                    only: mpi_comm
    use ieee_arithmetic
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_fluid_averaging_t) :: inlet_characteristic_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state
        procedure   :: compute_temporal_dft

    end type inlet_characteristic_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(inlet_characteristic_t),   intent(inout) :: self
        
        ! Set name, family
        call self%set_name("Inlet - Characteristic")
        call self%set_family("Inlet")

        ! Add functions
        call self%bcproperties%add('Total Pressure',    'Required')
        call self%bcproperties%add('Total Temperature', 'Required')
        call self%bcproperties%add('Normal-1',          'Required')
        call self%bcproperties%add('Normal-2',          'Required')
        call self%bcproperties%add('Normal-3',          'Required')

        ! Set default values
        call self%set_fcn_option('Total Pressure',    'val', 110000._rk)
        call self%set_fcn_option('Total Temperature', 'val', 300._rk)
        call self%set_fcn_option('Normal-1',          'val', 1._rk)
        call self%set_fcn_option('Normal-2',          'val', 0._rk)
        call self%set_fcn_option('Normal-3',          'val', 0._rk)

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
        class(inlet_characteristic_t),   intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,  p_m,                      &
            density_mean, vel1_mean, vel2_mean, vel3_mean, p_mean, T_mean, mom1_mean, mom2_mean, mom3_mean, energy_mean, &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,                     &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_m,    v_m,    w_m,                                                        &
            vel1_bc,   vel2_bc,   vel3_bc,                                                       &
            T_bc,   vmag2_m, vmag, f, df, dT, T, vel, veln, rminus, asp_ext, asp_int, M, ddensity, dvel1, dvel2, dvel3, dp, denergy, dmom1, dmom2, dmom3, &
            density_real, vel1_real, vel2_real, vel3_real, pressure_real,   &
            density_imag, vel1_imag, vel2_imag, vel3_imag, pressure_imag, c1, c2, c3, c4, c5, &
            PT, TT, n1, n2, n3, nmag, r, unorm_1, unorm_2, unorm_3


        type(AD_D)  :: vel1_avg, vel2_avg, vel3_avg, density_avg, pressure_avg, c_avg

        real(rk)    :: K0, u_axial, amp, omega
            
        integer(ik) :: ierr, igq

        logical :: converged


        call self%compute_temporal_dft(worker,bc_comm,density_real,  density_imag,    &
                                                      vel1_real,     vel1_imag,       &
                                                      vel2_real,     vel2_imag,       &
                                                      vel3_real,     vel3_imag,       &
                                                      pressure_real, pressure_imag)

        ! Get spatio-temporal averages
        density_avg = density_real(1)
        vel1_avg    = vel1_real(1)
        vel2_avg    = vel2_real(1)
        vel3_avg    = vel3_real(1)
        pressure_avg = pressure_real(1)
        !call self%compute_averages(worker,bc_COMM,vel1_avg, vel2_avg, vel3_avg, density_avg, pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)


        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        PT   = self%bcproperties%compute('Total Pressure',   worker%time(),worker%coords())
        TT   = self%bcproperties%compute('Total Temperature',worker%time(),worker%coords())

        ! Get user-input normal vector and normalize
        n1 = self%bcproperties%compute('Normal-1', worker%time(), worker%coords())
        n2 = self%bcproperties%compute('Normal-2', worker%time(), worker%coords())
        n3 = self%bcproperties%compute('Normal-3', worker%time(), worker%coords())

        !   Explicit allocation to handle GCC bug:
        !       GCC/GFortran Bugzilla Bug 52162 
        allocate(nmag(size(n1)), stat=ierr)
        if (ierr /= 0) call AllocationError

        nmag = sqrt(n1*n1 + n2*n2 + n3*n3)
        n1 = n1/nmag
        n2 = n2/nmag
        n3 = n3/nmag


        ! Interpolate interior solution to quadrature nodes
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')
        p_m       = worker%get_field('Pressure'  , 'value', 'face interior')



        grad1_density_m = worker%get_field('Density'   ,'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density'   ,'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density'   ,'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field('Momentum-1','grad1', 'face interior')
        grad2_mom1_m    = worker%get_field('Momentum-1','grad2', 'face interior')
        grad3_mom1_m    = worker%get_field('Momentum-1','grad3', 'face interior')

        grad1_mom2_m    = worker%get_field('Momentum-2','grad1', 'face interior')
        grad2_mom2_m    = worker%get_field('Momentum-2','grad2', 'face interior')
        grad3_mom2_m    = worker%get_field('Momentum-2','grad3', 'face interior')

        grad1_mom3_m    = worker%get_field('Momentum-3','grad1', 'face interior')
        grad2_mom3_m    = worker%get_field('Momentum-3','grad2', 'face interior')
        grad3_mom3_m    = worker%get_field('Momentum-3','grad3', 'face interior')
        
        grad1_energy_m  = worker%get_field('Energy'    ,'grad1', 'face interior')
        grad2_energy_m  = worker%get_field('Energy'    ,'grad2', 'face interior')
        grad3_energy_m  = worker%get_field('Energy'    ,'grad3', 'face interior')


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if

        ! Compute velocity components
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m

        ! Compute velocity magnitude squared from interior state
        vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
        vmag = sqrt(vmag2_m)

        ! Compute boundary condition velocity components from imposed direction
        vel1_mean = vmag*n1
        vel2_mean = vmag*n2
        vel3_mean = vmag*n3


        ! Allocate perturbation storage
        ddensity = ZERO*density_m
        dvel1    = ZERO*density_m
        dvel2    = ZERO*density_m
        dvel3    = ZERO*density_m
        dp       = ZERO*density_m


        ! Compute perturbation from time-mean
        do igq = 1,size(ddensity)
            ddensity(igq) = density_m(igq) - density_avg
            dvel1(igq)    = u_m(igq)       - vel1_avg
            dvel2(igq)    = v_m(igq)       - vel2_avg
            dvel3(igq)    = w_m(igq)       - vel3_avg
            dp(igq)       = p_m(igq)       - pressure_avg
        end do

        ! Compute 1D characteristics
        c1 = ZERO*density_m
        c2 = ZERO*density_m
        c3 = ZERO*density_m
        c4 = ZERO*density_m
        c5 = ZERO*density_m
        do igq = 1,size(ddensity)
            c1(igq) = -(   c_avg*c_avg   )*ddensity(igq) + dp(igq)
            c2(igq) =  (density_avg*c_avg)*dvel1(igq)
            c3(igq) =  (density_avg*c_avg)*dvel2(igq)
            c4(igq) =  (density_avg*c_avg)*dvel3(igq)    + dp(igq)
            c5(igq) = -(density_avg*c_avg)*dvel3(igq)    + dp(igq)
        end do

        ! From the characteristics computed, 1-4 should be 0 unless we want to set an incoming 
        ! perturbation. c5 is kept from the interior.
        c1 = ZERO
        c2 = ZERO
        c3 = ZERO
        c4 = ZERO

        ! Incoming pressure wave
        amp = 10._rk
        omega = worker%time_manager%freqs(1)
        do igq = 1,size(ddensity)
            c4(igq) = amp*cos(omega*worker%time())
        end do


        ! Reconstruct perturbation field
        do igq = 1,size(ddensity)
            ddensity(igq) = -(ONE/(c_avg*c_avg))*c1(igq) + (ONE/(TWO*c_avg*c_avg))*c4(igq) + (ONE/(TWO*c_avg*c_avg))*c5(igq)
            dvel1(igq)    = (ONE/(density_avg*c_avg))*c2(igq)
            dvel2(igq)    = (ONE/(density_avg*c_avg))*c3(igq)
            dvel3(igq)    = (ONE/(TWO*density_avg*c_avg))*c4(igq) - (ONE/(TWO*density_avg*c_avg))*c5(igq)
            dp(igq)       = HALF*c4(igq) + HALF*c5(igq)
        end do


        p_mean = ZERO*density_m
        p_mean(:) = pressure_avg
        T_mean = TT*(p_mean/PT)**((gam-ONE)/gam)
        density_mean = p_mean/(T_mean*Rgas)

        vmag = sqrt(TWO*cp*(TT-T_mean))
        vel1_mean = n1*vmag
        vel2_mean = n2*vmag
        vel3_mean = n3*vmag

        
        ! Assemble mean + perturbation
        density_bc = density_mean + ddensity
        vel1_bc    = vel1_mean    + dvel1
        vel2_bc    = vel2_mean    + dvel2
        vel3_bc    = vel3_mean    + dvel3
        p_bc       = p_mean       + dp

        ! Form conserved variables
        density_bc = density_bc
        mom1_bc    = density_bc*vel1_bc
        mom2_bc    = density_bc*vel2_bc
        mom3_bc    = density_bc*vel3_bc
        energy_bc  = p_bc/(gam - ONE) + (density_bc*HALF)*(vel1_bc*vel1_bc + &
                                                           vel2_bc*vel2_bc + &
                                                           vel3_bc*vel3_bc)


        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        ! Store computed boundary state
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


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



    end subroutine compute_bc_state
    !******************************************************************************************



    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag)
        class(inlet_characteristic_t),      intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        type(mpi_comm),                     intent(in)      :: bc_comm
        type(AD_D),     allocatable,        intent(inout)   :: density_Ft_real(:)
        type(AD_D),     allocatable,        intent(inout)   :: density_Ft_imag(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel1_Ft_real(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel1_Ft_imag(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel2_Ft_real(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel2_Ft_imag(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel3_Ft_real(:)
        type(AD_D),     allocatable,        intent(inout)   :: vel3_Ft_imag(:)
        type(AD_D),     allocatable,        intent(inout)   :: pressure_Ft_real(:)
        type(AD_D),     allocatable,        intent(inout)   :: pressure_Ft_imag(:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp


        type(AD_D), allocatable, dimension(:)     ::  &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr
        real(rk)    :: physical_nodes(1,3)

        physical_nodes(1,:) = [0.5_rk, 0.5_rk, ZERO]

        !ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime
        ntime   = worker%time_manager%ntime

        ! Allocate storage for discrete time instances
        allocate(density(ntime),  &
                 mom1(ntime),  &
                 mom2(ntime),  &
                 mom3(ntime),  &
                 vel1(ntime),  &
                 vel2(ntime),  &
                 vel3(ntime),  &
                 energy(ntime),  &
                 pressure(ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        
        ! Allocate storage for temporal dft
        allocate(density_Ft_real( ntime), density_Ft_imag( ntime),  &
                 vel1_Ft_real(    ntime), vel1_Ft_imag(    ntime),  &
                 vel2_Ft_real(    ntime), vel2_Ft_imag(    ntime),  &
                 vel3_Ft_real(    ntime), vel3_Ft_imag(    ntime),  &
                 pressure_Ft_real(ntime), pressure_Ft_imag(ntime), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Construct theta discretization
        do itime = 1,ntime
            ! Interpolate solution to physical_nodes at current radial station: [ntheta]
            density(itime:itime) = worker%interpolate_field_general('Density',    physical_nodes=physical_nodes, itime=itime)
            mom1(itime:itime)    = worker%interpolate_field_general('Momentum-1', physical_nodes=physical_nodes, itime=itime)
            mom2(itime:itime)    = worker%interpolate_field_general('Momentum-2', physical_nodes=physical_nodes, itime=itime)
            mom3(itime:itime)    = worker%interpolate_field_general('Momentum-3', physical_nodes=physical_nodes, itime=itime)
            energy(itime:itime)  = worker%interpolate_field_general('Energy',     physical_nodes=physical_nodes, itime=itime)

            !if (worker%coordinate_system() == 'Cylindrical') then
            !    mom2(:,itime) = mom2(:,itime)/self%r(iradius)  ! convert to tangential momentum
            !end if
        end do

        ! Compute velocities and pressure at each time
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density
        pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )


        call dft(density(:),  ZERO*density(:),  density_real_tmp,  density_imag_tmp )
        call dft(vel1(:),     ZERO*vel1(:),     vel1_real_tmp,     vel1_imag_tmp    )
        call dft(vel2(:),     ZERO*vel2(:),     vel2_real_tmp,     vel2_imag_tmp    )
        call dft(vel3(:),     ZERO*vel3(:),     vel3_real_tmp,     vel3_imag_tmp    )
        call dft(pressure(:), ZERO*pressure(:), pressure_real_tmp, pressure_imag_tmp)

        density_Ft_real(:)  = density_real_tmp
        density_Ft_imag(:)  = density_imag_tmp
        vel1_Ft_real(:)     = vel1_real_tmp
        vel1_Ft_imag(:)     = vel1_imag_tmp
        vel2_Ft_real(:)     = vel2_real_tmp
        vel2_Ft_imag(:)     = vel2_imag_tmp
        vel3_Ft_real(:)     = vel3_real_tmp
        vel3_Ft_imag(:)     = vel3_imag_tmp
        pressure_Ft_real(:) = pressure_real_tmp
        pressure_Ft_imag(:) = pressure_imag_tmp


    end subroutine compute_temporal_dft






end module bc_state_inlet_characteristic
