module bc_state_totalinlet
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
    type, public, extends(bc_state_t) :: totalinlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type totalinlet_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(totalinlet_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("Total Inlet")
        call self%set_family("Inlet")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("Momentum-1")
!        call self%set_equation("Momentum-2")
!        call self%set_equation("Momentum-3")
!        call self%set_equation("Energy"    )


        !
        ! Add functions
        !
        call self%bcproperties%add('Total Pressure',   'Required')
        call self%bcproperties%add('Total Temperature','Required')

        call self%bcproperties%add('X-Normal',         'Required')
        call self%bcproperties%add('Y-Normal',         'Required')
        call self%bcproperties%add('Z-Normal',         'Required')
        !call self%bcproperties%add('nr',              'Required')
        !call self%bcproperties%add('nt',              'Required')
        !call self%bcproperties%add('normal_direction','Required')


        !
        ! Set default angle
        !
        call self%set_fcn_option('X-Normal', 'val', 1._rk)
        call self%set_fcn_option('Y-Normal', 'val', 0._rk)
        call self%set_fcn_option('Z-Normal', 'val', 0._rk)

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
        class(totalinlet_t),    intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop



        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  p_m,    &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, p_bc,   &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            flux_x, flux_y, flux_z, integrand,                  &
            u_m,    v_m,    w_m,                                &
            u_bc,   v_bc,   w_bc,                               &
            T_bc,   vmag2_m, vmag, H_bc


        real(rk)                                    :: gam_m, cp_m, M, time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            TT, PT, nx, ny, nz, nr, nt, theta,          &
            normal_direction, x, y, r, normx, normy, normz



        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        coords = worker%coords()
        time   = worker%time()
        PT = self%bcproperties%compute("Total Pressure",     time, coords)
        TT = self%bcproperties%compute("Total Temperature",  time, coords)


        !normal_direction = self%bcproperties%compute("normal_direction",    time, coords)
        nx = self%bcproperties%compute("X-Normal", time, coords)
        ny = self%bcproperties%compute("Y-Normal", time, coords)
        nz = self%bcproperties%compute("Z-Normal", time, coords)



        !
        ! Interpolate interior solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face("Density"   , 'value', 'face interior')
        rhou_m = worker%get_primary_field_face("Momentum-1", 'value', 'face interior')
        rhov_m = worker%get_primary_field_face("Momentum-2", 'value', 'face interior')
        rhow_m = worker%get_primary_field_face("Momentum-3", 'value', 'face interior')
        rhoE_m = worker%get_primary_field_face("Energy"    , 'value', 'face interior')




        drho_dx_m  = worker%get_primary_field_face("Density"   ,'grad1', 'face interior')
        drho_dy_m  = worker%get_primary_field_face("Density"   ,'grad2', 'face interior')
        drho_dz_m  = worker%get_primary_field_face("Density"   ,'grad3', 'face interior')

        drhou_dx_m = worker%get_primary_field_face("Momentum-1",'grad1', 'face interior')
        drhou_dy_m = worker%get_primary_field_face("Momentum-1",'grad2', 'face interior')
        drhou_dz_m = worker%get_primary_field_face("Momentum-1",'grad3', 'face interior')

        drhov_dx_m = worker%get_primary_field_face("Momentum-2",'grad1', 'face interior')
        drhov_dy_m = worker%get_primary_field_face("Momentum-2",'grad2', 'face interior')
        drhov_dz_m = worker%get_primary_field_face("Momentum-2",'grad3', 'face interior')

        drhow_dx_m = worker%get_primary_field_face("Momentum-3",'grad1', 'face interior')
        drhow_dy_m = worker%get_primary_field_face("Momentum-3",'grad2', 'face interior')
        drhow_dz_m = worker%get_primary_field_face("Momentum-3",'grad3', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face("Energy"    ,'grad1', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face("Energy"    ,'grad2', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face("Energy"    ,'grad3', 'face interior')




        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute velocity components
        !
        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m



        !
        ! Compute velocity magnitude squared from interior state
        !
        vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
        vmag = sqrt(vmag2_m)


        !
        ! Compute boundary condition velocity components from imposed direction
        !
        u_bc = vmag*nx
        v_bc = vmag*ny
        w_bc = vmag*nz




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
        rho_bc = p_bc/(T_bc*287.15_rk)


        !
        ! Compute bc momentum
        !
        rhou_bc = rho_bc * u_bc
        rhov_bc = rho_bc * v_bc
        rhow_bc = rho_bc * w_bc


        !
        ! Compute bc energy
        !
        rhoE_bc = p_bc/(gam_m - ONE) + (rho_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )



        !
        ! Store computed boundary state
        !
        call worker%store_bc_state("Density"   ,rho_bc, 'value')
        call worker%store_bc_state("Momentum-1",rhou_bc,'value')
        call worker%store_bc_state("Momentum-2",rhov_bc,'value')
        call worker%store_bc_state("Momentum-3",rhow_bc,'value')
        call worker%store_bc_state("Energy"    ,rhoE_bc,'value')




        call worker%store_bc_state("Density"   , drho_dx_m, 'grad1')
        call worker%store_bc_state("Density"   , drho_dy_m, 'grad2')
        call worker%store_bc_state("Density"   , drho_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-1", drhou_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-1", drhou_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-1", drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-2", drhov_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-2", drhov_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-2", drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-3", drhow_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-3", drhow_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-3", drhow_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Energy"    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state("Energy"    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state("Energy"    , drhoE_dz_m, 'grad3')



    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_totalinlet
