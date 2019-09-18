module fluid_laplacian_anisotropic_av_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private

    
    !> Volume flux for Fluid laplacian_anisotropic_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: fluid_laplacian_anisotropic_av_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_laplacian_anisotropic_av_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_laplacian_anisotropic_av_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Fluid Laplacian Anisotropic AV Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid laplacian_anisotropic_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_laplacian_anisotropic_av_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            grad1, grad2, grad3, pgrad1, pgrad2, pgrad3, av1, av2, av3, &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r


       
        !
        ! Get Model fields:
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !
        av1 = worker%get_field('Artificial Viscosity - 1', 'value', 'element')
        av2 = worker%get_field('Artificial Viscosity - 2', 'value', 'element')
        av3 = worker%get_field('Artificial Viscosity - 3', 'value', 'element')


        




        !----------------------------------
        !            mass flux
        !----------------------------------
        grad1 = worker%get_field('Density'   , 'grad1', 'element')
        grad2 = worker%get_field('Density'   , 'grad2', 'element')
        grad3 = worker%get_field('Density'   , 'grad3', 'element')
        flux_1 = -av1*grad1
        flux_2 = -av2*grad2
        flux_3 = -av3*grad3
        
        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'density flux 1 is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'density flux 2 is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'density flux 3 is nan'
        
        call worker%integrate_volume_flux('Density','Diffusion',flux_1,flux_2,flux_3)



        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        grad1 = worker%get_field('Momentum-1'   , 'grad1', 'element')
        grad2 = worker%get_field('Momentum-1'   , 'grad2', 'element')
        grad3 = worker%get_field('Momentum-1'   , 'grad3', 'element')
        flux_1 = -av1*grad1
        flux_2 = -av2*grad2
        flux_3 = -av3*grad3
        
        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'mom1 flux 1 is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'mom1 flux 2 is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'mom1 flux 3 is nan'
        call worker%integrate_volume_flux('Momentum-1','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        grad1 = worker%get_field('Momentum-2'   , 'grad1', 'element')
        grad2 = worker%get_field('Momentum-2'   , 'grad2', 'element')
        grad3 = worker%get_field('Momentum-2'   , 'grad3', 'element')
        flux_1 = -av1*grad1
        flux_2 = -av2*grad2
        flux_3 = -av3*grad3
        

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if

        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'mom2 flux 1 is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'mom2 flux 2 is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'mom2 flux 3 is nan'
        call worker%integrate_volume_flux('Momentum-2','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        grad1 = worker%get_field('Momentum-3'   , 'grad1', 'element')
        grad2 = worker%get_field('Momentum-3'   , 'grad2', 'element')
        grad3 = worker%get_field('Momentum-3'   , 'grad3', 'element')

        if (any(ieee_is_nan(grad1(:)%x_ad_))) print *, 'mom3 grad 1 is nan'
        if (any(ieee_is_nan(grad2(:)%x_ad_))) print *, 'mom3 grad 2 is nan'
        if (any(ieee_is_nan(grad3(:)%x_ad_))) print *, 'mom3 grad 3 is nan'
        flux_1 = -av1*grad1
        flux_2 = -av2*grad2
        flux_3 = -av3*grad3
 
        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'mom3 flux 1 is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'mom3 flux 2 is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'mom3 flux 3 is nan'
        call worker%integrate_volume_flux('Momentum-3','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !           energy flux
        !----------------------------------
        grad1 = worker%get_field('Energy'   , 'grad1', 'element')
        grad2 = worker%get_field('Energy'   , 'grad2', 'element')
        grad3 = worker%get_field('Energy'   , 'grad3', 'element')

        if (any(ieee_is_nan(grad1(:)%x_ad_))) print *, 'energy grad 1 is nan'
        if (any(ieee_is_nan(grad2(:)%x_ad_))) print *, 'energy grad 2 is nan'
        if (any(ieee_is_nan(grad3(:)%x_ad_))) print *, 'energy grad 3 is nan'
        pgrad1 = worker%get_field('Pressure Gradient - 1'   , 'value', 'element')
        pgrad2 = worker%get_field('Pressure Gradient - 2'   , 'value', 'element')
        pgrad3 = worker%get_field('Pressure Gradient - 3'   , 'value', 'element')
        if (any(ieee_is_nan(pgrad1(:)%x_ad_))) print *, 'pressure grad 1 is nan'
        if (any(ieee_is_nan(pgrad2(:)%x_ad_))) print *, 'pressure grad 2 is nan'
        if (any(ieee_is_nan(pgrad3(:)%x_ad_))) print *, 'pressure grad 3 is nan'
        flux_1 = -av1*(grad1+pgrad1)
        flux_2 = -av2*(grad2+pgrad2)
        flux_3 = -av3*(grad3+pgrad3)
 
        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'energy flux 1 is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'energy flux 2 is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'energy flux 3 is nan'
        call worker%integrate_volume_flux('Energy','Diffusion',flux_1,flux_2,flux_3)

    end subroutine compute
    !*********************************************************************************************************






end module fluid_laplacian_anisotropic_av_volume_operator
