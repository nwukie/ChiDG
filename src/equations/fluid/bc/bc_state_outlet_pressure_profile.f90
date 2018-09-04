module bc_state_outlet_pressure_profile
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO, FOUR
    use mod_fluid,              only: gam, Rgas
    use mod_interpolation,      only: interpolate

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use ieee_arithmetic
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_pressure_profile_t

        real(rk), allocatable   :: r(:)
        real(rk), allocatable   :: p(:)
        logical :: profile_initialized

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: read_profile
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type outlet_pressure_profile_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_pressure_profile_t),   intent(inout) :: self
        
        call self%set_name("Outlet - Pressure Profile")
        call self%set_family("Outlet")
        self%profile_initialized = .false.

    end subroutine init
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/4/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine read_profile(self,worker,prop,bc_comm)
        class(outlet_pressure_profile_t),   intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM

        integer(ik)                 :: nr
        real(rk), dimension(1000)   :: r, p
        integer                     :: unit, msg
        logical                     :: file_exists

        namelist /profile/ nr, r, p


        !
        ! Check if input from 'outlet_profile.nml' is available.
        !   1: if available: read 
        !   2: if not available: error
        !
        inquire(file='outlet_profile.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='outlet_profile.nml')
            read(unit,nml=profile,iostat=msg)
            if (msg /= 0) call chidg_signal_one(FATAL,'outlet_profile%init: error reading outlet_profile.nml',msg)
            close(unit)
        else
            call chidg_signal(FATAL,'outlet_profile%init: could not find outlet_profile.nml')
        end if

        self%r = r(1:nr)        
        self%p = p(1:nr)        
        self%profile_initialized = .true.

    end subroutine read_profile
    !********************************************************************************





    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_pressure_profile_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_bc,   v_bc,    w_bc,  T_m, T_bc, p_bc, mom_norm
            

        real(rk),   allocatable, dimension(:) :: r, unorm1, unorm2, unorm3
        real(rk),   allocatable, dimension(:) :: p_input

        real(rk)    :: rho0, K0
        integer :: igq

        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        p_input = r
        do igq = 1,size(r)
            p_input(igq) = interpolate('linear',self%r,self%p,r(igq))
        end do


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density_m = worker%get_field('Density'    , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1' , 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2' , 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3' , 'value', 'face interior')
        energy_m  = worker%get_field('Energy'     , 'value', 'face interior')
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
        end if



!        K0   = 100._rk
!        rho0 = 1.2_rk
!        p_input = p_input(1)  +  0.5_rk*rho0*K0*K0*(ONE/(FOUR*FOUR)  -  ONE/(r*r))




        ! Extrapolate temperature and velocity
        T_bc = T_m
        u_bc = mom1_m/density_m
        v_bc = mom2_m/density_m
        w_bc = mom3_m/density_m


        ! Set user-specified constant pressure
        p_bc = density_m
        p_bc = p_input



        ! Compute density, momentum, energy
        density_bc = p_bc/(Rgas*T_bc)
        mom1_bc    = u_bc*density_bc
        mom2_bc    = v_bc*density_bc
        mom3_bc    = w_bc*density_bc
        energy_bc  = p_bc/(gam - ONE) + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
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
    !**************************************************************************************






end module bc_state_outlet_pressure_profile
