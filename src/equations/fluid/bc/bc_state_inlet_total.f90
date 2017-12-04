module bc_state_inlet_total
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, ZERO, HALF
    use mod_fluid,              only: Rgas, cp, gam
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
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: inlet_total_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type inlet_total_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(inlet_total_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("Inlet - Total")
        call self%set_family("Inlet")



        !
        ! Add functions
        !
        call self%bcproperties%add('Total Pressure',       'Required')
        call self%bcproperties%add('Total Temperature',    'Required')
        !call self%bcproperties%add('Density Perturbation', 'Required')

        call self%bcproperties%add('Normal-1',         'Required')
        call self%bcproperties%add('Normal-2',         'Required')
        call self%bcproperties%add('Normal-3',         'Required')


        !
        ! Set default values
        !
        call self%set_fcn_option('Total Pressure',    'val', 110000._rk)
        call self%set_fcn_option('Total Temperature', 'val', 300._rk)
        call self%set_fcn_option('Normal-1', 'val', 1._rk)
        call self%set_fcn_option('Normal-2', 'val', 0._rk)
        call self%set_fcn_option('Normal-3', 'val', 0._rk)

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
        class(inlet_total_t),   intent(inout)   :: self
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

!<<<<<<< HEAD
!            TT, n1, n2, n3, nmag, alpha, r, PT
!
!        real(rk)    :: K0, u_x
!=======
!
!        real(rk)    :: ntol, rafac
!>>>>>>> sulu/dev
            
        integer(ik) :: ierr, igq, inewton, nmax

        logical :: converged


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        PT   = self%bcproperties%compute('Total Pressure',        worker%time(), worker%coords())
        TT   = self%bcproperties%compute('Total Temperature',     worker%time(), worker%coords())
        !DRHO = self%bcproperties%compute('Density Perturbation',  worker%time(), worker%coords())


        !
        ! Get user-input normal vector and normalize
        !
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




        !
        ! Interpolate interior solution to quadrature nodes
        !
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



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if



!        !
!        ! Compute normal vector
!        !
!        K0 = 10._rk
!        u_x = 40._rk
!        alpha = atan2(K0/r,u_x)
!        
!        n1 = ZERO
!        n2 = sin(alpha)
!        n3 = cos(alpha)





        !
        ! Compute velocity components
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m



        !
        ! Compute velocity magnitude squared from interior state
        !
        vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
        vmag = sqrt(vmag2_m)


        !
        ! Compute boundary condition velocity components from imposed direction
        !
        u_bc = vmag*n1
        v_bc = vmag*n2
        w_bc = vmag*n3




        !
        ! Compute boundary condition temperature and pressure
        !
        T_bc = TT - (vmag2_m)/(TWO*cp)
        p_bc = PT*((T_bc/TT)**(gam/(gam-ONE)))

!        print*, 'initial: ', T_bc(1)%x_ad_
!
!        ! Reversed to be inward facing for formulation
!        unorm_1 = -worker%unit_normal(1)
!        unorm_2 = -worker%unit_normal(2)
!        unorm_3 = -worker%unit_normal(3)
!
!
!        veln = u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3
!        T = worker%get_field('Temperature', 'value', 'face interior')
!        asp_int = sqrt(gam*Rgas*T)
!        rminus = veln - TWO*asp_int/(gam-ONE)
!        rafac = ONE
!
!        converged = .false.
!        nmax = 20
!        ntol = 1.e-6_rk
!        inewton = 1
!        do while((inewton < nmax) .and. (.not. converged))
!            asp_ext = sqrt(gam*Rgas*T_bc)
!            vel = rafac * (rminus + TWO*asp_ext/(gam-ONE))
!
!            f  = TT - T_bc - HALF*(gam-ONE)*vel*vel/(gam*Rgas)
!            df = -(ONE + vel*rafac/asp_ext)
!            dT = -f/df
!            T_bc = T_bc + dT
!            converged = sqrt(sum(dT*dT)) < ntol
!            inewton = inewton + 1
!        end do
!
!        if (.not. converged) print*, 'Total Inlet boundary condition newton iteration did not converge.'
!
!        print*, 'final: ', T_bc(1)%x_ad_, inewton
!
!
!        !p_bc = PT*((T_bc/TT)**(gam/(gam-ONE)))
!
!        M = vmag / sqrt(gam*Rgas*T_bc)
!        p_bc = PT*(ONE + ((gam-ONE)/TWO)*M*M)**(-gam/(gam-ONE))


        !
        ! Compute boundary condition density from ideal gas law
        !
        density_bc = p_bc/(T_bc*Rgas)

        !
        ! Compute perturbation quantities
        !
        !density_bc = density_bc + drho



        !
        ! Compute bc momentum
        !
        mom1_bc = density_bc * u_bc
        mom2_bc = density_bc * v_bc
        mom3_bc = density_bc * w_bc



        !
        ! Compute bc energy
        !
        energy_bc = p_bc/(gam - ONE) + HALF*((mom1_bc*mom1_bc) + (mom2_bc*mom2_bc) + (mom3_bc*mom3_bc))/density_bc


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if




        !
        ! Store computed boundary state
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
    !******************************************************************************************









end module bc_state_inlet_total
