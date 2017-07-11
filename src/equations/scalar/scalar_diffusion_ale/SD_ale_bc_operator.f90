module SD_ale_bc_operator
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: SD_ale_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_ale_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_ale_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Scalar Diffusion ALE BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('u')

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
        class(SD_ale_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u, grad2_u, grad3_u,              &
            flux_1,  flux_2,  flux_3,               &
            integrand, mu

        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3

        type(AD_D), allocatable, dimension(:,:)   ::  &
            gradu, &
            flux_ref



       

        gradu = worker%get_primary_field_grad_ale_face('u', 'boundary')
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        !
        ! Compute scalar coefficient
        !
        mu = worker%get_model_field_face('Scalar Diffusion Coefficient', 'value', 'boundary')



        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = -mu*gradu(:,1)
        flux_2 = -mu*gradu(:,2)
        flux_3 = -mu*gradu(:,3)

        flux_ref = worker%post_process_boundary_diffusive_flux_ale(flux_1, flux_2, flux_3)
        !
        ! Compute boundary average flux
        !
        integrand = flux_ref(:,1)*norm_1 + flux_ref(:,2)*norm_2 + flux_ref(:,3)*norm_3




        if (any(ieee_is_nan(integrand(:)%x_ad_))) then
            print*, 'BC OP: ', integrand(:)%x_ad_
        end if


        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !**********************************************************************************************










end module SD_ale_bc_operator
