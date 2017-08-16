module euler_ale_bc_operator
    use mod_constants,      only: HALF,ONE, TWO
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: euler_ale_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler ALE BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

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
    subroutine compute(self,worker,prop)
        class(euler_ale_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            density_bc,  mom1_bc, mom2_bc, mom3_bc, energy_bc,  &
            u_bc,    v_bc,    w_bc,                             &
            H_bc,    p_bc,                                      &
            flux_x,  flux_y,  flux_z,  integrand

        type(AD_D), allocatable, dimension(:,:) :: flux

        real(rk),   allocatable, dimension(:)   :: norm_1, norm_2, norm_3
            
        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        density_bc = worker%get_primary_field_value_ale_face("Density"   , 'boundary')
        mom1_bc    = worker%get_primary_field_value_ale_face("Momentum-1", 'boundary')
        mom2_bc    = worker%get_primary_field_value_ale_face("Momentum-2", 'boundary')
        mom3_bc    = worker%get_primary_field_value_ale_face("Momentum-3", 'boundary')
        energy_bc  = worker%get_primary_field_value_ale_face("Energy"    , 'boundary')


        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        !
        ! Compute gamma
        !
        p_bc = worker%get_model_field_face('Pressure','value','boundary')


        !
        ! Compute velocity components
        !
        u_bc = mom1_bc/density_bc
        v_bc = mom2_bc/density_bc
        w_bc = mom3_bc/density_bc



        !
        ! Compute boundary condition energy and enthalpy
        !
        H_bc = (energy_bc + p_bc)/density_bc




        !=================================================
        ! Mass flux
        !=================================================
        flux_x = (density_bc * u_bc)
        flux_y = (density_bc * v_bc)
        flux_z = (density_bc * w_bc)
        flux = worker%post_process_boundary_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=density_bc, interp_source='face interior')

        integrand = flux(:,1)*norm_1 + flux(:,2)*norm_2 + flux(:,3)*norm_3

        call worker%integrate_boundary('Density',integrand)

        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = (density_bc * u_bc * u_bc) + p_bc
        flux_y = (density_bc * u_bc * v_bc) 
        flux_z = (density_bc * u_bc * w_bc) 
        flux = worker%post_process_boundary_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=mom1_bc, interp_source='face interior')

        integrand = flux(:,1)*norm_1 + flux(:,2)*norm_2 + flux(:,3)*norm_3

        call worker%integrate_boundary('Momentum-1',integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (density_bc * v_bc * u_bc) 
        flux_y = (density_bc * v_bc * v_bc) + p_bc
        flux_z = (density_bc * v_bc * w_bc) 
        flux = worker%post_process_boundary_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=mom2_bc, interp_source='face interior')

        integrand = flux(:,1)*norm_1 + flux(:,2)*norm_2 + flux(:,3)*norm_3

        call worker%integrate_boundary('Momentum-2',integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (density_bc * w_bc * u_bc) 
        flux_y = (density_bc * w_bc * v_bc) 
        flux_z = (density_bc * w_bc * w_bc) + p_bc
        flux = worker%post_process_boundary_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=mom3_bc, interp_source='face interior')

        integrand = flux(:,1)*norm_1 + flux(:,2)*norm_2 + flux(:,3)*norm_3

        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (density_bc * u_bc * H_bc) 
        flux_y = (density_bc * v_bc * H_bc) 
        flux_z = (density_bc * w_bc * H_bc) 
        flux = worker%post_process_boundary_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=energy_bc, interp_source='face interior')

        integrand = flux(:,1)*norm_1 + flux(:,2)*norm_2 + flux(:,3)*norm_3

        call worker%integrate_boundary('Energy',integrand)

    end subroutine compute
    !**********************************************************************************************























end module euler_ale_bc_operator
