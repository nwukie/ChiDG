module SA_ale_bc_operator
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: SA_ale_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_ale_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_ale_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Scalar Advection ALE BC Operator')

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
        class(SA_ale_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            u, c1, c2, c3,                          &
            flux_1, flux_2, flux_3, integrand, flux_1_ref, flux_2_ref, flux_3_ref

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid, testx


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid


        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3

 
        u_grid = worker%get_grid_velocity_face("u_grid")
        v_grid = worker%get_grid_velocity_face("v_grid")
        w_grid = worker%get_grid_velocity_face("w_grid")

!        print *, 'u_grid'
!        print *, u_grid
!        print *, 'v_grid'
!        print *, v_grid
!        print *, 'w_grid'
!        print *, w_grid
        jacobian_grid = worker%get_inv_jacobian_grid_face()
        det_jacobian_grid = worker%get_det_jacobian_grid_face()


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        u  = worker%get_primary_field_face('u', 'value', 'boundary')


        !
        ! Get model coefficients
        !
        c1 = worker%get_model_field_face('Scalar Advection Velocity-1', 'value', 'boundary')
        c2 = worker%get_model_field_face('Scalar Advection Velocity-2', 'value', 'boundary')
        c3 = worker%get_model_field_face('Scalar Advection Velocity-3', 'value', 'boundary')

        
        !
        ! Get normal vectors
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = (c1-u_grid)*u
        flux_2 = (c2-v_grid)*u
        flux_3 = (c3-w_grid)*u
        flux_1_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1 + jacobian_grid(:,1,2)*flux_2 + jacobian_grid(:,1,3)*flux_3)
        flux_2_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1 + jacobian_grid(:,2,2)*flux_2 + jacobian_grid(:,2,3)*flux_3)
        flux_3_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1 + jacobian_grid(:,3,2)*flux_2 + jacobian_grid(:,3,3)*flux_3)

        integrand = flux_1_ref*norm_1 + flux_2_ref*norm_2 + flux_3_ref*norm_3


        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !**********************************************************************************************










end module SA_ale_bc_operator
