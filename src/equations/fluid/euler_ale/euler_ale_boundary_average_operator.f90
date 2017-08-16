module euler_ale_boundary_average_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_ale_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler ALE Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

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



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_ale_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            density_m,  density_p,                  &
            mom1_m,     mom1_p,                     &
            mom2_m,     mom2_p,                     &
            mom3_m,     mom3_p,                     &
            energy_m,   energy_p,                   &
            p_m,        p_p,                        &
            enthalpy_m, enthalpy_p,                 &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p,       &
            flux_1,     flux_2,     flux_3,         &
            invdensity_m,   invdensity_p,           &
            integrand

        type(AD_D), allocatable,    dimension(:,:)  :: flux_ref_m, flux_ref_p


        !
        ! Interpolate solution to quadrature nodes
        !
        density_m = worker%get_primary_field_value_ale_face('Density'   , 'face interior')
        density_p = worker%get_primary_field_value_ale_face('Density'   , 'face exterior')

        mom1_m    = worker%get_primary_field_value_ale_face('Momentum-1', 'face interior')
        mom1_p    = worker%get_primary_field_value_ale_face('Momentum-1', 'face exterior')

        mom2_m    = worker%get_primary_field_value_ale_face('Momentum-2', 'face interior')
        mom2_p    = worker%get_primary_field_value_ale_face('Momentum-2', 'face exterior')

        mom3_m    = worker%get_primary_field_value_ale_face('Momentum-3', 'face interior')
        mom3_p    = worker%get_primary_field_value_ale_face('Momentum-3', 'face exterior')

        energy_m  = worker%get_primary_field_value_ale_face('Energy'    , 'face interior')
        energy_p  = worker%get_primary_field_value_ale_face('Energy'    , 'face exterior')

        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p




        !
        ! Compute pressure and total enthalpy
        !
        p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')


        enthalpy_m = (energy_m + p_m)*invdensity_m
        enthalpy_p = (energy_p + p_p)*invdensity_p



        !================================
        !       MASS FLUX
        !================================
!        flux_1_m = mom1_m
!        flux_2_m = mom2_m
!        flux_3_m = mom3_m
!        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_1_m,flux_2_m,flux_3_m, advected_quantity=density_m, interp_source='face interior')
!
!        flux_1_p = mom1_p
!        flux_2_p = mom2_p
!        flux_3_p = mom3_p
!        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_1_p,flux_2_p,flux_3_p, advected_quantity=density_p, interp_source='face exterior')
!
!        flux_1 = (flux_ref_m(:,1) + flux_ref_p(:,1))
!        flux_2 = (flux_ref_m(:,2) + flux_ref_p(:,2))
!        flux_3 = (flux_ref_m(:,3) + flux_ref_p(:,3))
!
!
!
!        ! dot with normal vector
!        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)
!
!        call worker%integrate_boundary('Density',integrand)


        flux_1_m = mom1_m
        flux_2_m = mom2_m
        flux_3_m = mom3_m

        flux_1_p = mom1_p
        flux_2_p = mom2_p
        flux_3_p = mom3_p

        call worker%integrate_boundary_average('Density','Advective',       &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !       X-MOMENTUM FLUX
        !================================
!        flux_1_m = (mom1_m*mom1_m)*invdensity_m + p_m
!        flux_2_m = (mom1_m*mom2_m)*invdensity_m
!        flux_3_m = (mom1_m*mom3_m)*invdensity_m
!        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_1_m,flux_2_m,flux_3_m, advected_quantity=mom1_m, interp_source='face interior')
!
!        flux_1_p = (mom1_p*mom1_p)*invdensity_p + p_p
!        flux_2_p = (mom1_p*mom2_p)*invdensity_p
!        flux_3_p = (mom1_p*mom3_p)*invdensity_p
!        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_1_p,flux_2_p,flux_3_p, advected_quantity=mom1_p, interp_source='face exterior')
!
!        flux_1 = (flux_ref_m(:,1) + flux_ref_p(:,1))
!        flux_2 = (flux_ref_m(:,2) + flux_ref_p(:,2))
!        flux_3 = (flux_ref_m(:,3) + flux_ref_p(:,3))
!
!
!        ! dot with normal vector
!        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)
!
!        call worker%integrate_boundary('Momentum-1',integrand)

        flux_1_m = (mom1_m*mom1_m)*invdensity_m + p_m
        flux_2_m = (mom1_m*mom2_m)*invdensity_m
        flux_3_m = (mom1_m*mom3_m)*invdensity_m

        flux_1_p = (mom1_p*mom1_p)*invdensity_p + p_p
        flux_2_p = (mom1_p*mom2_p)*invdensity_p
        flux_3_p = (mom1_p*mom3_p)*invdensity_p

        call worker%integrate_boundary_average('Momentum-1','Advective',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !================================
        !       Y-MOMENTUM FLUX
        !================================
!        flux_1_m = (mom2_m*mom1_m)*invdensity_m
!        flux_2_m = (mom2_m*mom2_m)*invdensity_m + p_m
!        flux_3_m = (mom2_m*mom3_m)*invdensity_m
!        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_1_m,flux_2_m,flux_3_m, advected_quantity=mom2_m, interp_source='face interior')
!
!        flux_1_p = (mom2_p*mom1_p)*invdensity_p
!        flux_2_p = (mom2_p*mom2_p)*invdensity_p + p_p
!        flux_3_p = (mom2_p*mom3_p)*invdensity_p
!        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_1_p,flux_2_p,flux_3_p, advected_quantity=mom2_p, interp_source='face exterior')
!
!        flux_1 = (flux_ref_m(:,1) + flux_ref_p(:,1))
!        flux_2 = (flux_ref_m(:,2) + flux_ref_p(:,2))
!        flux_3 = (flux_ref_m(:,3) + flux_ref_p(:,3))
!
!        ! dot with normal vector
!        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)
!
!        call worker%integrate_boundary('Momentum-2',integrand)

        flux_1_m = (mom2_m*mom1_m)*invdensity_m
        flux_2_m = (mom2_m*mom2_m)*invdensity_m + p_m
        flux_3_m = (mom2_m*mom3_m)*invdensity_m

        flux_1_p = (mom2_p*mom1_p)*invdensity_p
        flux_2_p = (mom2_p*mom2_p)*invdensity_p + p_p
        flux_3_p = (mom2_p*mom3_p)*invdensity_p

        call worker%integrate_boundary_average('Momentum-2','Advective',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !       Z-MOMENTUM FLUX
        !================================
!        flux_1_m = (mom3_m*mom1_m)*invdensity_m
!        flux_2_m = (mom3_m*mom2_m)*invdensity_m
!        flux_3_m = (mom3_m*mom3_m)*invdensity_m + p_m
!        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_1_m,flux_2_m,flux_3_m, advected_quantity=mom3_m, interp_source='face interior')
!
!        flux_1_p = (mom3_p*mom1_p)*invdensity_p
!        flux_2_p = (mom3_p*mom2_p)*invdensity_p
!        flux_3_p = (mom3_p*mom3_p)*invdensity_p + p_p
!        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_1_p,flux_2_p,flux_3_p, advected_quantity=mom3_p, interp_source='face exterior')
!
!        flux_1 = (flux_ref_m(:,1) + flux_ref_p(:,1))
!        flux_2 = (flux_ref_m(:,2) + flux_ref_p(:,2))
!        flux_3 = (flux_ref_m(:,3) + flux_ref_p(:,3))
!
!        ! dot with normal vector
!        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)
!
!        call worker%integrate_boundary('Momentum-3',integrand)

        flux_1_m = (mom3_m*mom1_m)*invdensity_m
        flux_2_m = (mom3_m*mom2_m)*invdensity_m
        flux_3_m = (mom3_m*mom3_m)*invdensity_m + p_m

        flux_1_p = (mom3_p*mom1_p)*invdensity_p
        flux_2_p = (mom3_p*mom2_p)*invdensity_p
        flux_3_p = (mom3_p*mom3_p)*invdensity_p + p_p

        call worker%integrate_boundary_average('Momentum-3','Advective',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !          ENERGY FLUX
        !================================
!        flux_1_m = mom1_m*enthalpy_m
!        flux_2_m = mom2_m*enthalpy_m
!        flux_3_m = mom3_m*enthalpy_m
!        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_1_m,flux_2_m,flux_3_m, advected_quantity=energy_m, interp_source='face interior')
!
!        flux_1_p = mom1_p*enthalpy_p
!        flux_2_p = mom2_p*enthalpy_p
!        flux_3_p = mom3_p*enthalpy_p
!        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_1_p,flux_2_p,flux_3_p, advected_quantity=energy_p, interp_source='face exterior')
!
!        flux_1 = (flux_ref_m(:,1) + flux_ref_p(:,1))
!        flux_2 = (flux_ref_m(:,2) + flux_ref_p(:,2))
!        flux_3 = (flux_ref_m(:,3) + flux_ref_p(:,3))
!
!        ! dot with normal vector
!        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)
!
!        call worker%integrate_boundary('Energy',integrand)

        flux_1_m = mom1_m*enthalpy_m
        flux_2_m = mom2_m*enthalpy_m
        flux_3_m = mom3_m*enthalpy_m

        flux_1_p = mom1_p*enthalpy_p
        flux_2_p = mom2_p*enthalpy_p
        flux_3_p = mom3_p*enthalpy_p

        call worker%integrate_boundary_average('Energy','Advective',        &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

    end subroutine compute
    !*********************************************************************************************************












end module euler_ale_boundary_average_operator
