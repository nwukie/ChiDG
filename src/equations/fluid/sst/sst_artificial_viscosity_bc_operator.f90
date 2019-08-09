module sst_artificial_viscosity_bc_operator
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use mod_sst, only: sst_avc
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: sst_artificial_viscosity_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_artificial_viscosity_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_artificial_viscosity_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('SST Artificial Viscosity BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density * k'   )
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
        class(sst_artificial_viscosity_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            sigma_k, mu_neg, grad1, grad2, grad3, pgrad1, pgrad2, pgrad3, av, &
            flux_1, flux_2, flux_3, integrand


        integer(ik) :: p
        real(rk) :: avc, h(3)
        real(rk), allocatable :: h_s(:,:)

        mu_neg = worker%get_field('SST Artificial Viscosity k-neg'   , 'value', 'boundary')
        sigma_k = worker%get_field('SST sigma_k'   , 'value', 'boundary')
        
        h_s = worker%h_smooth()
        
        
        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        h = h/real(p+1,rk)

        
        h_s = h_s/real(p+1,rk)

        grad1 = worker%get_field('Density * k'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * k'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * k'   , 'grad3', 'boundary')
        flux_1 = -(h_s(:,1)**TWO*sst_avc+mu_neg*sigma_k)*grad1
        flux_2 = -(h_s(:,2)**TWO*sst_avc+mu_neg*sigma_k)*grad2
        flux_3 = -(h_s(:,3)**TWO*sst_avc+mu_neg*sigma_k)*grad3
        
        


        call worker%integrate_boundary_condition('Density * k','Diffusion',flux_1,flux_2,flux_3)
    end subroutine compute
    !**********************************************************************************************









end module sst_artificial_viscosity_bc_operator
