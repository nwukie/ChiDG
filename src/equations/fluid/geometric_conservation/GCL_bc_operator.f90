module GCL_bc_operator
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
    type, public, extends(operator_t) :: GCL_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type GCL_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/10/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(GCL_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Geometric Conservation BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('g_bar')

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/10/2017 
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(GCL_bc_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            g_bar, flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            det_jacobian_grid,   &
            norm_1, norm_2, norm_3

        real(rk),   allocatable, dimension(:,:) :: grid_velocity

        real(rk),   allocatable, dimension(:,:,:) :: &
            inv_jacobian_grid



        ! Just to set up AD variables
        g_bar = worker%get_primary_field_face('g_bar', 'value', 'boundary')


        !
        ! Get model coefficients
        !
!        grid_velocity_1   = worker%get_grid_velocity_face('u_grid','face interior')
!        grid_velocity_2   = worker%get_grid_velocity_face('v_grid','face interior')
!        grid_velocity_3   = worker%get_grid_velocity_face('w_grid','face interior')
        grid_velocity     = worker%get_grid_velocity_face('face interior')
        det_jacobian_grid = worker%get_det_jacobian_grid_face('value','face interior')
        inv_jacobian_grid = worker%get_inv_jacobian_grid_face('face interior')


        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)




        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = g_bar  !just to initialize AD allocation
        flux_2 = g_bar  !just to initialize AD allocation
        flux_3 = g_bar  !just to initialize AD allocation

        flux_1 = (inv_jacobian_grid(:,1,1)*grid_velocity(:,1) + inv_jacobian_grid(:,1,2)*grid_velocity(:,2) + inv_jacobian_grid(:,1,3)*grid_velocity(:,3))*det_jacobian_grid
        flux_2 = (inv_jacobian_grid(:,2,1)*grid_velocity(:,1) + inv_jacobian_grid(:,2,2)*grid_velocity(:,2) + inv_jacobian_grid(:,2,3)*grid_velocity(:,3))*det_jacobian_grid
        flux_3 = (inv_jacobian_grid(:,3,1)*grid_velocity(:,1) + inv_jacobian_grid(:,3,2)*grid_velocity(:,2) + inv_jacobian_grid(:,3,3)*grid_velocity(:,3))*det_jacobian_grid

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        !
        ! Integrate volume flux
        !
        integrand = -integrand
        call worker%integrate_boundary('g_bar',integrand)


    end subroutine compute
    !**********************************************************************************************










end module GCL_bc_operator
