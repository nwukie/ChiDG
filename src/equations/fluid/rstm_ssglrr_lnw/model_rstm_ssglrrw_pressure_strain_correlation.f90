!>
!! Pressure-Strain Correlation Model
!!              
!!
!! @author Eric M. Wolf
!! @date   01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_rstm_ssglrrw_pressure_strain_correlation
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
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
    type, extends(model_t)  :: rstm_ssglrrw_pressure_strain_correlation_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_pressure_strain_correlation_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric Wolf 
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_pressure_strain_correlation_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Pressure-Strain Correlation')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Pressure-Strain-11')
        call self%add_model_field('Pressure-Strain-22')
        call self%add_model_field('Pressure-Strain-33')
        call self%add_model_field('Pressure-Strain-12')
        call self%add_model_field('Pressure-Strain-13')
        call self%add_model_field('Pressure-Strain-23')



    end subroutine init
    !***************************************************************************************




    !>
    !! Computes the Pressure-Strain Correlation tensor
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_pressure_strain_correlation_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density, mu, nu, density_nutilde, mu_t, nutilde, lamda_t,   &
            c1, c1_star, c2, c3, c3_star, c4, c5,                       &
            epsilon_t, k_t,                                             &
            pressure_strain_11, pressure_strain_22, pressure_strain_33, &
            pressure_strain_12, pressure_strain_13, pressure_strain_23, &
            strain_rate_11, strain_rate_22, strain_rate_33,             &
            strain_rate_12, strain_rate_13, strain_rate_23,             &
            anisotropy_11, anisotropy_22, anisotropy_33,             &
            anisotropy_12, anisotropy_13, anisotropy_23,             &
            production_11, production_22, production_33,             &
            production_12, production_13, production_23,             &
            production_trace, strain_rate_trace,                        &
            anisotropy_product_sum, anisotropy_strain_rate_product_sum,   &
            rotation_rate_12, rotation_rate_13, rotation_rate_23

        !
        ! Interpolate solution to quadrature nodes
        !
        
        density = worker%get_field('Density',    'value')

        ! Get blended coefficients
        c1      = worker%get_field('RSTMSSGLRRW C1',        'value')
        c1_star = worker%get_field('RSTMSSGLRRW C1_star',   'value')
        c2      = worker%get_field('RSTMSSGLRRW C2',        'value')
        c3      = worker%get_field('RSTMSSGLRRW C3',        'value')
        c3_star = worker%get_field('RSTMSSGLRRW C3_star',   'value')
        c4      = worker%get_field('RSTMSSGLRRW C4',        'value')
        c5      = worker%get_field('RSTMSSGLRRW C5',        'value')


        ! Get turbulence quantities
        k_t         = worker%get_field('Turbulence Kinetic Energy',     'value')
        epsilon_t   = worker%get_field('Turbulence Isotropic Dissipation Rate',   'value')


        
        strain_rate_11 = worker%get_field('Strain Rate-11',    'value')
        strain_rate_22 = worker%get_field('Strain Rate-22',    'value')
        strain_rate_33 = worker%get_field('Strain Rate-33',    'value')
        strain_rate_12 = worker%get_field('Strain Rate-12',    'value')
        strain_rate_13 = worker%get_field('Strain Rate-13',    'value')
        strain_rate_23 = worker%get_field('Strain Rate-23',    'value')
        strain_rate_trace = worker%get_field('Strain Rate-Trace',    'value')

        anisotropy_11 = worker%get_field('Anisotropy-11',    'value')
        anisotropy_22 = worker%get_field('Anisotropy-22',    'value')
        anisotropy_33 = worker%get_field('Anisotropy-33',    'value')
        anisotropy_12 = worker%get_field('Anisotropy-12',    'value')
        anisotropy_13 = worker%get_field('Anisotropy-13',    'value')
        anisotropy_23 = worker%get_field('Anisotropy-23',    'value')

        rotation_rate_12 = worker%get_field('Rotation Rate-12',    'value')
        rotation_rate_13 = worker%get_field('Rotation Rate-13',    'value')
        rotation_rate_23 = worker%get_field('Rotation Rate-23',    'value')
        
        production_11 = worker%get_field('Production-11',    'value')
        production_22 = worker%get_field('Production-22',    'value')
        production_33 = worker%get_field('Production-33',    'value')
        production_12 = worker%get_field('Production-12',    'value')
        production_13 = worker%get_field('Production-13',    'value')
        production_23 = worker%get_field('Production-23',    'value')
        production_trace = worker%get_field('Production-Trace',    'value')



        anisotropy_product_sum = &
                                    anisotropy_11*anisotropy_11 +anisotropy_12*anisotropy_12 +anisotropy_13*anisotropy_13 + &
                                    anisotropy_12*anisotropy_12 +anisotropy_22*anisotropy_22 +anisotropy_23*anisotropy_23 + &
                                    anisotropy_13*anisotropy_13 +anisotropy_23*anisotropy_23 +anisotropy_33*anisotropy_33

        anisotropy_strain_rate_product_sum = &
                                    anisotropy_11*strain_rate_11 +anisotropy_12*strain_rate_12 +anisotropy_13*strain_rate_13 + &
                                    anisotropy_12*strain_rate_12 +anisotropy_22*strain_rate_22 +anisotropy_23*strain_rate_23 + &
                                    anisotropy_13*strain_rate_13 +anisotropy_23*strain_rate_23 +anisotropy_33*strain_rate_33



        !
        ! NOTE: The anisotropy and strain rate tensors are symmetric, while the rotation rate tensor is antisymmetric.
        !       The antisymmetry of the rotation rate tensor is reflected in the negative signs appearing in the terms
        !       multiplied by c5.
        !

        ! Diagonal terms
        pressure_strain_11 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_11 + &
            c2*density*epsilon_t*(anisotropy_11*anisotropy_11 + anisotropy_12*anisotropy_12 + anisotropy_13*anisotropy_13 - & 
            (ONE/THREE)*anisotropy_product_sum) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*(strain_rate_11 -(ONE/THREE)*strain_rate_trace)+ &
            c4*density*k_t*(anisotropy_11*strain_rate_11 + anisotropy_12*strain_rate_12 + anisotropy_13*strain_rate_13 + &
                            anisotropy_11*strain_rate_11 + anisotropy_12*strain_rate_12 + anisotropy_13*strain_rate_13 - &
                            (TWO/THREE)*anisotropy_strain_rate_product_sum) +&
            c5*density*k_t*( anisotropy_12*rotation_rate_12 + anisotropy_13*rotation_rate_13 + &
                                 anisotropy_12*rotation_rate_12 + anisotropy_13*rotation_rate_13)

        pressure_strain_22 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_22 + &
            c2*density*epsilon_t*(anisotropy_12*anisotropy_12 + anisotropy_22*anisotropy_22 + anisotropy_23*anisotropy_23 - & 
            (ONE/THREE)*anisotropy_product_sum) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*(strain_rate_22 -(ONE/THREE)*strain_rate_trace) + &
            c4*density*k_t*(anisotropy_12*strain_rate_12 + anisotropy_22*strain_rate_22 + anisotropy_23*strain_rate_23 + &
                            anisotropy_12*strain_rate_12 + anisotropy_22*strain_rate_22 + anisotropy_23*strain_rate_23 - &
                            (TWO/THREE)*anisotropy_strain_rate_product_sum) +&
            c5*density*k_t*(-anisotropy_12*rotation_rate_12  + anisotropy_23*rotation_rate_23  &
                                -anisotropy_12*rotation_rate_12 +  anisotropy_23*rotation_rate_23)

        pressure_strain_33 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_33 + &
            c2*density*epsilon_t*(anisotropy_13*anisotropy_13 + anisotropy_23*anisotropy_23 + anisotropy_33*anisotropy_33 - & 
            (ONE/THREE)*anisotropy_product_sum) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*(strain_rate_33 -(ONE/THREE)*strain_rate_trace) + &
            c4*density*k_t*(anisotropy_13*strain_rate_13 + anisotropy_23*strain_rate_23 + anisotropy_33*strain_rate_33 + &
                            anisotropy_13*strain_rate_13 + anisotropy_23*strain_rate_23 + anisotropy_33*strain_rate_33 - &
                            (TWO/THREE)*anisotropy_strain_rate_product_sum) +&
            c5*density*k_t*(-anisotropy_13*rotation_rate_13  - anisotropy_23*rotation_rate_23   &
                                -anisotropy_13*rotation_rate_13 - anisotropy_23*rotation_rate_23 )
        ! Off-diagonal terms 
        pressure_strain_12 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_12 + &
            c2*density*epsilon_t*(anisotropy_11*anisotropy_12 + anisotropy_12*anisotropy_22 + anisotropy_13*anisotropy_23  & 
            ) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*strain_rate_12 + &
            c4*density*k_t*(anisotropy_11*strain_rate_12 + anisotropy_12*strain_rate_22 + anisotropy_13*strain_rate_23 + &
                            anisotropy_12*strain_rate_11 + anisotropy_22*strain_rate_12 + anisotropy_23*strain_rate_13  &
                            ) +&
            c5*density*k_t*(-anisotropy_11*rotation_rate_12  + anisotropy_13*rotation_rate_23 + &
                                 anisotropy_22*rotation_rate_12 + anisotropy_23*rotation_rate_13)

        pressure_strain_13 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_13 + &
            c2*density*epsilon_t*(anisotropy_11*anisotropy_13 + anisotropy_12*anisotropy_23 + anisotropy_13*anisotropy_33  & 
            ) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*strain_rate_13 + &
            c4*density*k_t*(anisotropy_11*strain_rate_13 + anisotropy_12*strain_rate_23 + anisotropy_13*strain_rate_33 + &
                            anisotropy_13*strain_rate_11 + anisotropy_23*strain_rate_12 + anisotropy_33*strain_rate_13  &
                            ) +&
            c5*density*k_t*(-anisotropy_11*rotation_rate_13 - anisotropy_12*rotation_rate_23  + &
                               anisotropy_23*rotation_rate_12 + anisotropy_33*rotation_rate_13)

        pressure_strain_23 = -(c1*density*epsilon_t + HALF*c1_star*production_trace)*anisotropy_23 + &
            c2*density*epsilon_t*(anisotropy_12*anisotropy_13 + anisotropy_22*anisotropy_23 + anisotropy_23*anisotropy_33  & 
            ) + &
            (c3-c3_star*sqrt(anisotropy_product_sum))*density*k_t*strain_rate_23 + &
            c4*density*k_t*(anisotropy_12*strain_rate_13 + anisotropy_22*strain_rate_23 + anisotropy_23*strain_rate_33 + &
                            anisotropy_13*strain_rate_12 + anisotropy_23*strain_rate_22 + anisotropy_33*strain_rate_23  &
                            ) +&
            c5*density*k_t*(-anisotropy_12*rotation_rate_13 - anisotropy_22*rotation_rate_23   &
                                -anisotropy_13*rotation_rate_12 + anisotropy_33*rotation_rate_23)


        call worker%store_model_field('Pressure-Strain-11', 'value', pressure_strain_11)
        call worker%store_model_field('Pressure-Strain-22', 'value', pressure_strain_22)
        call worker%store_model_field('Pressure-Strain-33', 'value', pressure_strain_33)
        call worker%store_model_field('Pressure-Strain-12', 'value', pressure_strain_12)
        call worker%store_model_field('Pressure-Strain-13', 'value', pressure_strain_13)
        call worker%store_model_field('Pressure-Strain-23', 'value', pressure_strain_23)




    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_pressure_strain_correlation
