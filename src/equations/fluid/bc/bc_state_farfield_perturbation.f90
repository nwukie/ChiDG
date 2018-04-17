module bc_state_farfield_perturbation
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ONE, RKTOL, FOUR
    use mod_fluid,              only: Rgas, gam
    use mod_dft,                only: dft

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
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: farfield_perturbation_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state
        procedure   :: compute_temporal_dft

    end type farfield_perturbation_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(farfield_perturbation_t),   intent(inout) :: self
        
        ! Set operator name, family
        call self%set_name("Farfield Perturbation")
        call self%set_family("Farfield")

        ! Add boundary condition parameters
        call self%bcproperties%add('Density',   'Required')
        call self%bcproperties%add('Pressure',  'Required')
        call self%bcproperties%add('Velocity-1','Required')
        call self%bcproperties%add('Velocity-2','Required')
        call self%bcproperties%add('Velocity-3','Required')

        ! Set default parameter values
        call self%set_fcn_option('Density',    'val', 1.2_rk)
        call self%set_fcn_option('Pressure',   'val', 100000._rk)
        call self%set_fcn_option('Velocity-1', 'val', 0._rk)
        call self%set_fcn_option('Velocity-2', 'val', 0._rk)
        call self%set_fcn_option('Velocity-3', 'val', 0._rk)

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
        class(farfield_perturbation_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,     &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            normal_momentum, normal_velocity, R_inf, R_extrapolated,    &
            u_bc_norm, v_bc_norm, w_bc_norm, u_bc_tang, v_bc_tang,      &
            w_bc_tang, entropy_bc, c_bc, c_m, p_m, T_m, u_m, v_m, w_m,  &
            density_real, vel1_real, vel2_real, vel3_real, pressure_real,   &
            density_imag, vel1_imag, vel2_imag, vel3_imag, pressure_imag,   &
            ddensity, dvel1, dvel2, dvel3, dpressure, dmom1, dmom2, dmom3, denergy, c1, c2, c3, c4, c5, vel1_bc, vel2_bc, vel3_bc, &
            density_avg, vel1_avg, vel2_avg, vel3_avg, pressure_avg, c_avg

        real(rk), allocatable, dimension(:) ::  &
            unorm_1, unorm_2, unorm_3, r,       &
            rho_input, p_input, u_input,        &
            v_input, w_input, T_input, c_input, &
            u_grid, v_grid, w_grid

        real(rk)    :: amp, omega
        integer(ik) :: igq

        logical, allocatable, dimension(:)  ::  &
            inflow, outflow


        !
        ! Get boundary condition input parameters
        !
        rho_input = self%bcproperties%compute('Density',    worker%time(), worker%coords())
        p_input   = self%bcproperties%compute('Pressure',   worker%time(), worker%coords())
        u_input   = self%bcproperties%compute('Velocity-1', worker%time(), worker%coords())
        v_input   = self%bcproperties%compute('Velocity-2', worker%time(), worker%coords())
        w_input   = self%bcproperties%compute('Velocity-3', worker%time(), worker%coords())
        T_input   = p_input/(rho_input*Rgas)
        c_input   = T_input ! allocate to silence error on DEBUG build
        c_input   = sqrt(gam*Rgas*T_input)



        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if



        !
        ! Get Pressure, Temperature from interior
        !
        p_m = worker%get_field('Pressure',    'value', 'face interior')
        T_m = worker%get_field('Temperature', 'value', 'face interior')
        c_m = sqrt(gam*Rgas*T_m)


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



        ! Initialize arrays
        density_bc = density_m
        mom1_bc    = mom1_m
        mom2_bc    = mom2_m
        mom3_bc    = mom3_m
        energy_bc  = energy_m
        R_inf      = density_m
        R_extrapolated = density_m
        u_bc_norm = density_m
        v_bc_norm = density_m
        w_bc_norm = density_m
        u_bc_tang = density_m
        v_bc_tang = density_m
        w_bc_tang = density_m
        entropy_bc = density_m


        !
        ! Get unit normal vector
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)




        !
        ! Dot momentum with normal vector
        !
        normal_momentum = mom1_m*unorm_1 + mom2_m*unorm_2 + mom3_m*unorm_3


        !
        ! Determine which nodes are inflow/outflow
        !
        inflow  = ( normal_momentum <= RKTOL )
        outflow = ( normal_momentum >  RKTOL )


        !
        ! Compute internal velocities
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m



        !
        ! Compute Riemann invariants
        !
        R_inf          = (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3) - TWO*c_input/(gam - ONE)
        R_extrapolated = (u_m*unorm_1     + v_m*unorm_2     + w_m*unorm_3    ) + TWO*c_m/(gam - ONE)


        !
        ! Compute boundary velocities
        !
        c_bc = ((gam - ONE)/FOUR)*(R_extrapolated - R_inf)

        u_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_1
        v_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_2
        w_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_3



        !
        ! Compute tangential velocities
        !
        where (inflow)

            u_bc_tang = u_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_1
            v_bc_tang = v_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_2
            w_bc_tang = w_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_3

            entropy_bc = p_input/(rho_input**gam)

        elsewhere !outflow

            u_bc_tang = u_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_1
            v_bc_tang = v_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_2
            w_bc_tang = w_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_3

            entropy_bc = p_m/(density_m**gam)

        end where


        ! Compute boundary state
        density_bc  = (c_bc*c_bc/(entropy_bc*gam))**(ONE/(gam-ONE))
        mom1_bc     = (u_bc_norm + u_bc_tang)*density_bc
        mom2_bc     = (v_bc_norm + v_bc_tang)*density_bc
        mom3_bc     = (w_bc_norm + w_bc_tang)*density_bc

        p_bc   = (density_bc**gam)*entropy_bc
        energy_bc = (p_bc/(gam - ONE)) + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc



        ! Override with strict inputs
        density_bc = rho_input
        mom1_bc    = rho_input*u_input
        mom2_bc    = rho_input*v_input
        mom3_bc    = rho_input*w_input
        energy_bc  = (p_input/(gam - ONE)) + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc



!        amp = 10._rk
!        omega = worker%time_manager%freqs(1)
!
!        c1 = ZERO*density_m
!        c2 = ZERO*density_m
!        c3 = ZERO*density_m
!        c4 = ZERO*density_m
!        c5 = ZERO*density_m
!
!        c1 = -c_input*c_input*rho_input + p_input
!        c2 =  rho_input*c_input*u_input
!        c3 =  rho_input*c_input*v_input
!        c4 =  rho_input*c_input*w_input + p_input
!        c5 = -rho_input*c_input*w_input + p_input
!        !c4 = amp*cos(omega*worker%time())
!        c4 = ZERO
!
!        density_bc = (-ONE/(c_input*c_input))*c1  +  (ONE/(TWO*c_input*c_input))*c4  +  (ONE/(TWO*c_input*c_input))*c5
!        vel1_bc    = (ONE/(rho_input*c_input))*c2
!        vel2_bc    = (ONE/(rho_input*c_input))*c3
!        vel3_bc    = (ONE/(TWO*rho_input*c_input))*c4 - (ONE/(TWO*rho_input*c_input))*c5
!        mom1_bc = density_bc*vel1_bc
!        mom2_bc = density_bc*vel2_bc
!        mom3_bc = density_bc*vel3_bc
!        p_bc       = HALF*c4 + HALF*c5
!        energy_bc  = (p_bc/(gam - ONE)) + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc
!
!        print*, density_bc(:)%x_ad_
!        print*, vel1_bc(:)%x_ad_
!        print*, vel3_bc(:)%x_ad_
!        print*, p_bc(:)%x_ad_







!        call self%compute_temporal_dft(worker,bc_comm,density_real,  density_imag,    &
!                                                 vel1_real,     vel1_imag,       &
!                                                 vel2_real,     vel2_imag,       &
!                                                 vel3_real,     vel3_imag,       &
!                                                 pressure_real, pressure_imag)
!
!        density_avg = density_real(1)
!        vel1_avg    = vel1_real(1)
!        vel2_avg    = vel2_real(1)
!        vel3_avg    = vel3_real(1)
!        pressure_avg = pressure_real(1)
!        c_avg = sqrt(gam*pressure_avg/density_avg)

        density_avg = ZERO*density_m
        vel1_avg    = ZERO*density_m
        vel2_avg    = ZERO*density_m
        vel3_avg    = ZERO*density_m
        pressure_avg = ZERO*density_m
        density_avg = rho_input
        vel1_avg = u_input
        vel2_avg = v_input
        vel3_avg = w_input
        pressure_avg = p_input
        c_avg = sqrt(gam*pressure_avg/density_avg)

        ddensity  = ZERO*density_m
        dvel1     = ZERO*density_m
        dvel2     = ZERO*density_m
        dvel3     = ZERO*density_m
        dmom1     = ZERO*density_m
        dmom2     = ZERO*density_m
        dmom3     = ZERO*density_m
        denergy   = ZERO*density_m
        dpressure = ZERO*density_m

        ! Compute perturbation
        amp = 10._rk
        omega = worker%time_manager%freqs(1)
        !do igq = 1,size(ddensity)
        !    ddensity(igq)  = (ONE/(TWO*c_avg*c_avg))*amp*sin(omega*worker%time())
        !    dvel3(igq)     = (ONE/(TWO*density_avg*c_avg))*amp*sin(omega*worker%time())
        !    dpressure(igq) = HALF*amp*sin(omega*worker%time())
        !end do
        ddensity  = (ONE/(TWO*c_avg*c_avg))*amp*cos(omega*worker%time())
        dvel3     = (ONE/(TWO*density_avg*c_avg))*amp*cos(omega*worker%time())
        dpressure = HALF*amp*cos(omega*worker%time())

        dmom1 = ddensity*dvel1
        dmom2 = ddensity*dvel2
        dmom3 = ddensity*dvel3
        denergy = dpressure/(gam-ONE) + HALF*(dmom1*dmom1 + dmom2*dmom2 + dmom3*dmom3)/ddensity

        print*, density_bc(:)%x_ad_
        print*, mom1_bc(:)%x_ad_
        print*, mom3_bc(:)%x_ad_
        print*, p_bc(:)%x_ad_

        density_bc = density_bc + ddensity
        mom1_bc    = mom1_bc    + dmom1
        mom2_bc    = mom2_bc    + dmom2
        mom3_bc    = mom3_bc    + dmom3
        energy_bc  = energy_bc  + denergy


!        c5 = ZERO*density_m
!        do igq = 1,size(ddensity)
!            dvel3(igq) = mom3_m(igq)/density_m(igq) - vel3_avg
!            dpressure(igq) = p_m(igq) - pressure_avg
!            c5(igq) = -density_avg*c_avg*dvel3(igq) + dpressure(igq)
!            ddensity(igq)  = (ONE/(TWO*c_avg*c_avg))*c5(igq)
!            dvel1(igq)     = ZERO
!            dvel2(igq)     = ZERO
!            dvel3(igq)     = -(ONE/(TWO*density_avg*c_avg))*c5(igq)
!            dpressure(igq) = HALF*c5(igq)
!        end do
!        dmom1 = ddensity*dvel1
!        dmom2 = ddensity*dvel2
!        dmom3 = ddensity*dvel3
!        denergy = dpressure/(gam-ONE) + HALF*(dmom1*dmom1 + dmom2*dmom2 + dmom3*dmom3)/ddensity
!
!
!        density_bc = density_bc + ddensity
!        mom1_bc    = mom1_bc    + dmom1
!        mom2_bc    = mom2_bc    + dmom2
!        mom3_bc    = mom3_bc    + dmom3
!        energy_bc  = energy_bc  + denergy


        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
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
    !*****************************************************************************************



    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag)
        class(farfield_perturbation_t),     intent(inout)   :: self
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

        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime

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





end module bc_state_farfield_perturbation
