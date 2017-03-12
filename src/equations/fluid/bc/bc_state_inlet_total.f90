module bc_state_inlet_total
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
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
        call self%bcproperties%add('Total Pressure',   'Required')
        call self%bcproperties%add('Total Temperature','Required')

        call self%bcproperties%add('Normal-1',         'Required')
        call self%bcproperties%add('Normal-2',         'Required')
        call self%bcproperties%add('Normal-3',         'Required')


        !
        ! Set default angle
        !
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
    subroutine compute_bc_state(self,worker,prop)
        class(inlet_total_t),   intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,  p_m,      &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,     &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            u_m,    v_m,    w_m,                                        &
            u_bc,   v_bc,   w_bc,                                       &
            T_bc,   vmag2_m, vmag, H_bc


        integer(ik)                                 :: ierr
        real(rk)                                    :: gam_m, cp_m, M
        real(rk),       allocatable, dimension(:)   ::  &
            TT, n1, n2, n3, nmag, alpha, r, PT
            


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        PT = self%bcproperties%compute('Total Pressure',     worker%time(), worker%coords())
        TT = self%bcproperties%compute('Total Temperature',  worker%time(), worker%coords())


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
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')



        drho_dx_m  = worker%get_primary_field_face('Density'   ,'grad1', 'face interior')
        drho_dy_m  = worker%get_primary_field_face('Density'   ,'grad2', 'face interior')
        drho_dz_m  = worker%get_primary_field_face('Density'   ,'grad3', 'face interior')

        drhou_dx_m = worker%get_primary_field_face('Momentum-1','grad1', 'face interior')
        drhou_dy_m = worker%get_primary_field_face('Momentum-1','grad2', 'face interior')
        drhou_dz_m = worker%get_primary_field_face('Momentum-1','grad3', 'face interior')

        drhov_dx_m = worker%get_primary_field_face('Momentum-2','grad1', 'face interior')
        drhov_dy_m = worker%get_primary_field_face('Momentum-2','grad2', 'face interior')
        drhov_dz_m = worker%get_primary_field_face('Momentum-2','grad3', 'face interior')

        drhow_dx_m = worker%get_primary_field_face('Momentum-3','grad1', 'face interior')
        drhow_dy_m = worker%get_primary_field_face('Momentum-3','grad2', 'face interior')
        drhow_dz_m = worker%get_primary_field_face('Momentum-3','grad3', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face('Energy'    ,'grad1', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face('Energy'    ,'grad2', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face('Energy'    ,'grad3', 'face interior')



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
        end if



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
        !& HARDCODED GAMMA. HARDCODED CP
        gam_m = 1.4_rk
        cp_m  = 287.15_rk*(gam_m/(gam_m-ONE))


        T_bc = TT - (vmag2_m)/(TWO*cp_m)
        p_bc = PT*((T_bc/TT)**(gam_m/(gam_m-ONE)))


        !
        ! Compute boundary condition density from ideal gas law
        !
        density_bc = p_bc/(T_bc*287.15_rk)


        !
        ! Compute bc momentum
        !
        mom1_bc = density_bc * u_bc
        mom2_bc = density_bc * v_bc
        mom3_bc = density_bc * w_bc


        !
        ! Compute bc energy
        !
        energy_bc = p_bc/(gam_m - ONE) + (density_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * worker%coordinate('1','boundary')
        end if


        !
        ! Store computed boundary state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')




        call worker%store_bc_state('Density'   , drho_dx_m,  'grad1')
        call worker%store_bc_state('Density'   , drho_dy_m,  'grad2')
        call worker%store_bc_state('Density'   , drho_dz_m,  'grad3')
                                                
        call worker%store_bc_state('Momentum-1', drhou_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-1', drhou_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-1', drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-2', drhov_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-2', drhov_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-2', drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-3', drhow_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-3', drhow_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-3', drhow_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Energy'    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state('Energy'    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state('Energy'    , drhoE_dz_m, 'grad3')



    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_inlet_total
