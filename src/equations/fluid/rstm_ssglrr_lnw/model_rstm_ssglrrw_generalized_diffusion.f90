module model_rstm_ssglrrw_generalized_diffusion
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
    use mod_rstm_ssglrrw
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_generalized_diffusion_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_generalized_diffusion_t
    !***************************************************************************************





contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_generalized_diffusion_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Generalized Diffusion')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Diffusion-111')
        call self%add_model_field('Diffusion-221')
        call self%add_model_field('Diffusion-331')
        call self%add_model_field('Diffusion-121')
        call self%add_model_field('Diffusion-131')
        call self%add_model_field('Diffusion-231')

        call self%add_model_field('Diffusion-112')
        call self%add_model_field('Diffusion-222')
        call self%add_model_field('Diffusion-332')
        call self%add_model_field('Diffusion-122')
        call self%add_model_field('Diffusion-132')
        call self%add_model_field('Diffusion-232')

        call self%add_model_field('Diffusion-113')
        call self%add_model_field('Diffusion-223')
        call self%add_model_field('Diffusion-333')
        call self%add_model_field('Diffusion-123')
        call self%add_model_field('Diffusion-133')
        call self%add_model_field('Diffusion-233')


    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_generalized_diffusion_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            omega, density, mu_t, mu, D, dissipation_coefficient,              &
            grad1_reynolds_11, grad2_reynolds_11, grad3_reynolds_11,    &
            grad1_reynolds_22, grad2_reynolds_22, grad3_reynolds_22,    &
            grad1_reynolds_33, grad2_reynolds_33, grad3_reynolds_33,    &
            grad1_reynolds_12, grad2_reynolds_12, grad3_reynolds_12,    &
            grad1_reynolds_13, grad2_reynolds_13, grad3_reynolds_13,    &
            grad1_reynolds_23, grad2_reynolds_23, grad3_reynolds_23,    &
            density_reynolds_11, density_reynolds_22, density_reynolds_33, &
            density_reynolds_12, density_reynolds_13, density_reynolds_23, &
            diffusion_111, diffusion_112, diffusion_113,                &
            diffusion_221, diffusion_222, diffusion_223,                &
            diffusion_331, diffusion_332, diffusion_333,                &
            diffusion_121, diffusion_122, diffusion_123,                &
            diffusion_131, diffusion_132, diffusion_133,                &
            diffusion_231, diffusion_232, diffusion_233                

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')

        omega       = worker%get_field('Omega', 'value')
        mu          = worker%get_field('Laminar Viscosity', 'value')
        D           = worker%get_field('RSTMSSGLRRW D_GD', 'value')

        dissipation_coefficient =  mu

        
        grad1_reynolds_11 = worker%get_field('Reynolds-11 - Gradient 1',    'value')
        grad2_reynolds_11 = worker%get_field('Reynolds-11 - Gradient 2',    'value')
        grad3_reynolds_11 = worker%get_field('Reynolds-11 - Gradient 3',    'value')
        grad1_reynolds_22 = worker%get_field('Reynolds-22 - Gradient 1',    'value')
        grad2_reynolds_22 = worker%get_field('Reynolds-22 - Gradient 2',    'value')
        grad3_reynolds_22 = worker%get_field('Reynolds-22 - Gradient 3',    'value')
        grad1_reynolds_33 = worker%get_field('Reynolds-33 - Gradient 1',    'value')
        grad2_reynolds_33 = worker%get_field('Reynolds-33 - Gradient 2',    'value')
        grad3_reynolds_33 = worker%get_field('Reynolds-33 - Gradient 3',    'value')
        grad1_reynolds_12 = worker%get_field('Reynolds-12 - Gradient 1',    'value')
        grad2_reynolds_12 = worker%get_field('Reynolds-12 - Gradient 2',    'value')
        grad3_reynolds_12 = worker%get_field('Reynolds-12 - Gradient 3',    'value')
        grad1_reynolds_13 = worker%get_field('Reynolds-13 - Gradient 1',    'value')
        grad2_reynolds_13 = worker%get_field('Reynolds-13 - Gradient 2',    'value')
        grad3_reynolds_13 = worker%get_field('Reynolds-13 - Gradient 3',    'value')
        grad1_reynolds_23 = worker%get_field('Reynolds-23 - Gradient 1',    'value')
        grad2_reynolds_23 = worker%get_field('Reynolds-23 - Gradient 2',    'value')
        grad3_reynolds_23 = worker%get_field('Reynolds-23 - Gradient 3',    'value')

        density_reynolds_11 = worker%get_field('Density * Reynolds-11', 'value')
        density_reynolds_22 = worker%get_field('Density * Reynolds-22', 'value')
        density_reynolds_33 = worker%get_field('Density * Reynolds-33', 'value')
        density_reynolds_12 = worker%get_field('Density * Reynolds-12', 'value')
        density_reynolds_13 = worker%get_field('Density * Reynolds-13', 'value')
        density_reynolds_23 = worker%get_field('Density * Reynolds-23', 'value')


        diffusion_111   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_11 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_11 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_11 

        diffusion_112   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_11 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_11 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_11 

        diffusion_113   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_11 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_11 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_11 

        diffusion_221   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_22 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_22 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_22 

        diffusion_222   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_22 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_22 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_22 

        diffusion_223   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_22 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_22 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_22 

        diffusion_331   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_33 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_33 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_33 

        diffusion_332   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_33 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_33 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_33 

        diffusion_333   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_33 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_33 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_33 

        diffusion_121   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_12 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_12 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_12 

        diffusion_122   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_12 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_12 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_12 

        diffusion_123   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_12 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_12 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_12 

        diffusion_131   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_13 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_13 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_13 

        diffusion_132   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_13 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_13 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_13 

        diffusion_133   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_13 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_13 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_13 

        diffusion_231   =   (dissipation_coefficient + D*density_reynolds_11/(SSG_LRRW_cmu*omega) )*grad1_reynolds_23 + & 
                            (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad2_reynolds_23 + & 
                            (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad3_reynolds_23 

        diffusion_232   =   (D*density_reynolds_12/(SSG_LRRW_cmu*omega) )*grad1_reynolds_23 + & 
                            (dissipation_coefficient + D*density_reynolds_22/(SSG_LRRW_cmu*omega) )*grad2_reynolds_23 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad3_reynolds_23 

        diffusion_233   =   (D*density_reynolds_13/(SSG_LRRW_cmu*omega) )*grad1_reynolds_23 + & 
                            (D*density_reynolds_23/(SSG_LRRW_cmu*omega) )*grad2_reynolds_23 + & 
                            (dissipation_coefficient + D*density_reynolds_33/(SSG_LRRW_cmu*omega) )*grad3_reynolds_23 



        call worker%store_model_field('Diffusion-111', 'value', diffusion_111)
        call worker%store_model_field('Diffusion-112', 'value', diffusion_112)
        call worker%store_model_field('Diffusion-113', 'value', diffusion_113)

        call worker%store_model_field('Diffusion-221', 'value', diffusion_221)
        call worker%store_model_field('Diffusion-222', 'value', diffusion_222)
        call worker%store_model_field('Diffusion-223', 'value', diffusion_223)

        call worker%store_model_field('Diffusion-331', 'value', diffusion_331)
        call worker%store_model_field('Diffusion-332', 'value', diffusion_332)
        call worker%store_model_field('Diffusion-333', 'value', diffusion_333)

        call worker%store_model_field('Diffusion-121', 'value', diffusion_121)
        call worker%store_model_field('Diffusion-122', 'value', diffusion_122)
        call worker%store_model_field('Diffusion-123', 'value', diffusion_123)

        call worker%store_model_field('Diffusion-131', 'value', diffusion_131)
        call worker%store_model_field('Diffusion-132', 'value', diffusion_132)
        call worker%store_model_field('Diffusion-133', 'value', diffusion_133)

        call worker%store_model_field('Diffusion-231', 'value', diffusion_231)
        call worker%store_model_field('Diffusion-232', 'value', diffusion_232)
        call worker%store_model_field('Diffusion-233', 'value', diffusion_233)



    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_generalized_diffusion
