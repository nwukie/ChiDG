module fluid_laplacian_anisotropic_av_boundary_average_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use ieee_arithmetic
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
    type, extends(operator_t), public :: fluid_laplacian_anisotropic_av_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_laplacian_anisotropic_av_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_laplacian_anisotropic_av_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Fluid Laplacian Anisotropic AV Boundary Average Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_laplacian_anisotropic_av_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                        &
            grad1_m, grad2_m, grad3_m, pgrad1_m, pgrad2_m, pgrad3_m,    &
            grad1_p, grad2_p, grad3_p, pgrad1_p, pgrad2_p, pgrad3_p,    &
            av1_m, av1_p,                                                 &
            av2_m, av2_p,                                                 &
            av3_m, av3_p,                                                 &
            flux_1_m, flux_2_m, flux_3_m,                               &
            flux_1_p, flux_2_p, flux_3_p


        real(rk), allocatable, dimension(:) :: r
!
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
        end if



        av1_m  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 1'    , 'value', 'face interior')
        av1_p  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 1'    , 'value', 'face exterior')

        av2_m  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 2'    , 'value', 'face interior')
        av2_p  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 2'    , 'value', 'face exterior')

        av3_m  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 3'    , 'value', 'face interior')
        av3_p  = worker%get_field('Smoothed Anisotropic Artificial Viscosity - 3'    , 'value', 'face exterior')

        !----------------------------------
        !            mass flux
        !----------------------------------
        grad1_m = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density'   , 'grad3', 'face interior')
        flux_1_m = -av1_m*grad1_m
        flux_2_m = -av2_m*grad2_m
        flux_3_m = -av3_m*grad3_m

        grad1_p = worker%get_field('Density'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density'   , 'grad3', 'face exterior')
        flux_1_p = -av1_p*grad1_p
        flux_2_p = -av2_p*grad2_p
        flux_3_p = -av3_p*grad3_p


        call worker%integrate_boundary_average('Density','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        grad1_m = worker%get_field('Momentum-1'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Momentum-1'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Momentum-1'   , 'grad3', 'face interior')
        flux_1_m = -av1_m*grad1_m
        flux_2_m = -av2_m*grad2_m
        flux_3_m = -av3_m*grad3_m

        grad1_p = worker%get_field('Momentum-1'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Momentum-1'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Momentum-1'   , 'grad3', 'face exterior')
        flux_1_p = -av1_p*grad1_p
        flux_2_p = -av2_p*grad2_p
        flux_3_p = -av3_p*grad3_p


        call worker%integrate_boundary_average('Momentum-1','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        grad1_m = worker%get_field('Momentum-2'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Momentum-2'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Momentum-2'   , 'grad3', 'face interior')
        flux_1_m = -av1_m*grad1_m
        flux_2_m = -av2_m*grad2_m
        flux_3_m = -av3_m*grad3_m

        grad1_p = worker%get_field('Momentum-2'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Momentum-2'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Momentum-2'   , 'grad3', 'face exterior')
        flux_1_p = -av1_p*grad1_p
        flux_2_p = -av2_p*grad2_p
        flux_3_p = -av3_p*grad3_p
        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1_m = flux_1_m * r
            flux_2_m = flux_2_m * r
            flux_3_m = flux_3_m * r

            flux_1_p = flux_1_p * r
            flux_2_p = flux_2_p * r
            flux_3_p = flux_3_p * r
        end if

        call worker%integrate_boundary_average('Momentum-2','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        grad1_m = worker%get_field('Momentum-3'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Momentum-3'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Momentum-3'   , 'grad3', 'face interior')
        flux_1_m = -av1_m*grad1_m
        flux_2_m = -av2_m*grad2_m
        flux_3_m = -av3_m*grad3_m

        grad1_p = worker%get_field('Momentum-3'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Momentum-3'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Momentum-3'   , 'grad3', 'face exterior')
        flux_1_p = -av1_p*grad1_p
        flux_2_p = -av2_p*grad2_p
        flux_3_p = -av3_p*grad3_p

        
        call worker%integrate_boundary_average('Momentum-3','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           energy flux
        !----------------------------------
        grad1_m = worker%get_field('Energy'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Energy'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Energy'   , 'grad3', 'face interior')
        pgrad1_m = worker%get_field('Pressure Gradient - 1'   , 'value', 'face interior')
        pgrad2_m = worker%get_field('Pressure Gradient - 2'   , 'value', 'face interior')
        pgrad3_m = worker%get_field('Pressure Gradient - 3'   , 'value', 'face interior')

        if (any(ieee_is_nan(pgrad1_m(:)%x_ad_))) print *, 'pressure grad 1_m is nan'
        if (any(ieee_is_nan(pgrad2_m(:)%x_ad_))) print *, 'pressure grad 2_m is nan'
        if (any(ieee_is_nan(pgrad3_m(:)%x_ad_))) print *, 'pressure grad 3_m is nan'
        flux_1_m = -av1_m*(grad1_m+pgrad1_m)
        flux_2_m = -av2_m*(grad2_m+pgrad2_m)
        flux_3_m = -av3_m*(grad3_m+pgrad3_m)

        grad1_p = worker%get_field('Energy'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Energy'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Energy'   , 'grad3', 'face exterior')

        pgrad1_p = worker%get_field('Pressure Gradient - 1'   , 'value', 'face exterior')
        pgrad2_p = worker%get_field('Pressure Gradient - 2'   , 'value', 'face exterior')
        pgrad3_p = worker%get_field('Pressure Gradient - 3'   , 'value', 'face exterior')
        if (any(ieee_is_nan(pgrad1_p(:)%x_ad_))) print *, 'pressure grad 1_p is nan'
        if (any(ieee_is_nan(pgrad2_p(:)%x_ad_))) print *, 'pressure grad 2_p is nan'
        if (any(ieee_is_nan(pgrad3_p(:)%x_ad_))) print *, 'pressure grad 3_p is nan'
        flux_1_p = -av1_p*(grad1_p+pgrad1_p)
        flux_2_p = -av2_p*(grad2_p+pgrad2_p)
        flux_3_p = -av3_p*(grad3_p+pgrad3_p)
        if (any(ieee_is_nan(flux_1_p(:)%x_ad_))) print *, 'energy flux 1_p is nan'
        if (any(ieee_is_nan(flux_2_p(:)%x_ad_))) print *, 'energy flux 2_p is nan'
        if (any(ieee_is_nan(flux_3_p(:)%x_ad_))) print *, 'energy flux 3_p is nan'

        if (any(ieee_is_nan(flux_1_m(:)%x_ad_))) print *, 'energy flux 1_m is nan'
        if (any(ieee_is_nan(flux_2_m(:)%x_ad_))) print *, 'energy flux 2_m is nan'
        if (any(ieee_is_nan(flux_3_m(:)%x_ad_))) print *, 'energy flux 3_m is nan'
       
        call worker%integrate_boundary_average('Energy','Diffusion',        &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

    end subroutine compute
    !*********************************************************************************************************












end module fluid_laplacian_anisotropic_av_boundary_average_operator
