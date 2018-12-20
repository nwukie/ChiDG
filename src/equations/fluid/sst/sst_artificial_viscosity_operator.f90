module sst_artificial_viscosity_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use mod_sst, only: sst_avc
    use ieee_arithmetic
    implicit none

    private

    
    !> Volume flux for Fluid laplacian_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_artificial_viscosity_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_artificial_viscosity_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_artificial_viscosity_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('SST Artificial Viscosity Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density * k'   )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid laplacian_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_artificial_viscosity_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            sigma_k, mu_neg, grad1, grad2, grad3, &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r

        integer(ik) :: p
        real(rk)    :: h(3), avc


       
        !
        ! Get Model fields:
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !


        
        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        h = h/real(p+1,rk)

        mu_neg = worker%get_field('SST Artificial Viscosity k-neg'   , 'value', 'element')
        sigma_k = worker%get_field('SST sigma_k'   , 'value', 'element')
        grad1 = worker%get_field('Density * k'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * k'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * k'   , 'grad3', 'element')
        flux_1 = -(h(1)**TWO*sst_avc-mu_neg*sigma_k)*grad1
        flux_2 = -(h(2)**TWO*sst_avc-mu_neg*sigma_k)*grad2
        flux_3 = -(h(3)**TWO*sst_avc-mu_neg*sigma_k)*grad3
        
        
        call worker%integrate_volume_flux('Density * k','Diffusion',flux_1,flux_2,flux_3)


    end subroutine compute
    !*********************************************************************************************************






end module sst_artificial_viscosity_operator
