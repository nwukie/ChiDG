module bc_state_inlet_total_profile
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, ZERO, HALF
    use mod_fluid,              only: Rgas, cp, gam
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use mod_interpolation,      only: interpolate
    use ieee_arithmetic
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   9/2/2018
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: inlet_total_profile_t

        real(rk), allocatable   :: r(:)
        real(rk), allocatable   :: p0(:)
        real(rk), allocatable   :: T0(:)
        real(rk), allocatable   :: n1(:)
        real(rk), allocatable   :: n2(:)
        real(rk), allocatable   :: n3(:)
        logical :: profile_initialized

    contains

        procedure   :: init
        procedure   :: read_profile
        procedure   :: compute_bc_state

    end type inlet_total_profile_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/2/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(inlet_total_profile_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name("Inlet - Total Profile")
        call self%set_family("Inlet")

        self%profile_initialized = .false.

    end subroutine init
    !********************************************************************************


    
    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/3/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine read_profile(self,worker,prop,bc_comm)
        class(inlet_total_profile_t),   intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        integer(ik)                 :: nr
        real(rk), dimension(1000)   :: r, p0, T0, n1, n2, n3
        integer                     :: unit, msg
        logical                     :: file_exists

        namelist /profile/ nr, r, p0, T0, n1, n2, n3


        !
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%mu
        !   2: if not available, do nothing and mu retains default value
        !
        inquire(file='inlet_total_profile.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='inlet_total_profile.nml')
            read(unit,nml=profile,iostat=msg)
            if (msg /= 0) call chidg_signal_one(FATAL,'inlet_total_profile%init: error reading inlet_total_profile.nml',msg)
            close(unit)
        else
            call chidg_signal(FATAL,'inlet_total_profile%init: could not find inlet_total_profile.nml')
        end if

        self%r  = r(1:nr)        
        self%p0 = p0(1:nr)        
        self%T0 = T0(1:nr)        
        self%n1 = n1(1:nr)        
        self%n2 = n2(1:nr)        
        self%n3 = n3(1:nr)        

        self%profile_initialized = .true.

    end subroutine read_profile
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/2/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(inlet_total_profile_t),   intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,  p_m,                      &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,                     &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_m,    v_m,    w_m,                                                        &
            u_bc,   v_bc,   w_bc,                                                       &
            T_bc,   vmag2_m, vmag, f, df, dT, T, vel, veln, rminus, asp_ext, asp_int, M

        real(rk),       allocatable, dimension(:)   ::  &
            PT, TT, DRHO, n1, n2, n3, nmag, alpha, r, unorm_1, unorm_2, unorm_3


        real(rk)    :: K0, u_axial
            
        integer(ik) :: ierr, igq, inewton, nmax

        logical :: converged

        if (.not. self%profile_initialized) call self%read_profile(worker,prop,bc_COMM)

        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        PT = r
        TT = r
        n1 = r
        n2 = r
        n3 = r
        do igq = 1,size(r)
            PT(igq) = interpolate('linear',self%r,self%p0,r(igq))
            TT(igq) = interpolate('linear',self%r,self%T0,r(igq))
            n1(igq) = interpolate('linear',self%r,self%n1,r(igq))
            n2(igq) = interpolate('linear',self%r,self%n2,r(igq))
            n3(igq) = interpolate('linear',self%r,self%n3,r(igq))
        end do


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
        u_bc = vmag*n1
        v_bc = vmag*n2
        w_bc = vmag*n3


        ! Compute boundary condition temperature and pressure
        T_bc = TT - (vmag2_m)/(TWO*cp)
        p_bc = PT*((T_bc/TT)**(gam/(gam-ONE)))


        ! Compute boundary condition density from ideal gas law
        density_bc = p_bc/(T_bc*Rgas)


        ! Compute bc momentum
        mom1_bc = density_bc * u_bc
        mom2_bc = density_bc * v_bc
        mom3_bc = density_bc * w_bc


        ! Compute bc energy
        energy_bc = p_bc/(gam - ONE) + HALF*((mom1_bc*mom1_bc) + (mom2_bc*mom2_bc) + (mom3_bc*mom3_bc))/density_bc


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









end module bc_state_inlet_total_profile
