module euler_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_fluid,          only: omega
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: euler_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler BC Flux")

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
        class(euler_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop

        ! data at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,   &
            u_bc, v_bc, w_bc, H_bc, p_bc, u_t, v_t, w_t,        &
            flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3, r
            
        real(rk) :: gam_bc



        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        density_bc = worker%get_primary_field_face('Density'   , 'value', 'boundary')
        mom1_bc    = worker%get_primary_field_face('Momentum-1', 'value', 'boundary')
        mom2_bc    = worker%get_primary_field_face('Momentum-2', 'value', 'boundary')
        mom3_bc    = worker%get_primary_field_face('Momentum-3', 'value', 'boundary')
        energy_bc  = worker%get_primary_field_face('Energy'    , 'value', 'boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc / worker%coordinate('1','boundary')
        end if





        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)



        !
        ! Compute gamma
        !
        p_bc   = worker%get_model_field_face('Pressure','value','boundary')
        gam_bc = 1.4_rk



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



        !
        ! Compute transport velocity
        !
        r = worker%coordinate('1','boundary') 
        u_t = u_bc
        v_t = v_bc - omega*r
        w_t = w_bc


        !=================================================
        ! mass flux
        !=================================================
        flux_1 = (density_bc * u_t )
        flux_2 = (density_bc * v_t )
        flux_3 = (density_bc * w_t )

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3


        call worker%integrate_boundary('Density',integrand)

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (density_bc * u_bc * u_t) + p_bc
        flux_2 = (density_bc * u_bc * v_t)
        flux_3 = (density_bc * u_bc * w_t)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Momentum-1',integrand)

        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = (density_bc * v_bc * u_t)
        flux_2 = (density_bc * v_bc * v_t) + p_bc
        flux_3 = (density_bc * v_bc * w_t)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * worker%coordinate('1','boundary')
        end if

        call worker%integrate_boundary('Momentum-2',integrand)

        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (density_bc * w_bc * u_t)
        flux_2 = (density_bc * w_bc * v_t)
        flux_3 = (density_bc * w_bc * w_t) + p_bc

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3


        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! energy flux
        !=================================================
        flux_1 = (density_bc * H_bc * u_t)
        flux_2 = (density_bc * H_bc * v_t)  +  r*omega*p_bc
        flux_3 = (density_bc * H_bc * w_t)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Energy',integrand)

    end subroutine compute
    !**********************************************************************************************























end module euler_bc_operator
