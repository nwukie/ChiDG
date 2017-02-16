module bc_state_fluid_extrapolate
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF

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
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: fluid_extrapolate_t

    contains

        procedure   :: init                 !< Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     !< boundary condition function implementation

    end type fluid_extrapolate_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_extrapolate_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Fluid Extrapolate")
        call self%set_family("Symmetry")


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


    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(fluid_extrapolate_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::      &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            flux_x, flux_y,  flux_z,  integrand,        &
            u_bc,   v_bc,    w_bc,                      &
            H_bc


        real(rk)                                    :: time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            p_bc, normx, normy, normz


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        rho_m  = worker%get_primary_field_face("Density"   ,'value', 'face interior')
        rhou_m = worker%get_primary_field_face("Momentum-1",'value', 'face interior')
        rhov_m = worker%get_primary_field_face("Momentum-2",'value', 'face interior')
        rhow_m = worker%get_primary_field_face("Momentum-3",'value', 'face interior')
        rhoE_m = worker%get_primary_field_face("Energy"    ,'value', 'face interior')



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




        !
        ! Store boundary condition state
        !
        call worker%store_bc_state("Density"   ,rho_m, 'value')
        call worker%store_bc_state("Momentum-1",rhou_m,'value')
        call worker%store_bc_state("Momentum-2",rhov_m,'value')
        call worker%store_bc_state("Momentum-3",rhow_m,'value')
        call worker%store_bc_state("Energy"    ,rhoE_m,'value')




        call worker%store_bc_state("Density"   , drho_dx_m,  'grad1')
        call worker%store_bc_state("Density"   , drho_dy_m,  'grad2')
        call worker%store_bc_state("Density"   , drho_dz_m,  'grad3')

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
    !**********************************************************************************************






end module bc_state_fluid_extrapolate
