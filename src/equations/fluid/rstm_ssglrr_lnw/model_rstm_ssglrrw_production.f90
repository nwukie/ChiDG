module model_rstm_ssglrrw_production
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
    type, extends(model_t)  :: rstm_ssglrrw_production_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_production_t
    !***************************************************************************************





contains




    !>  Model to compute the production term in the Reynolds Stress transport equation.
    !!  Note that this term is computed exactly from mean flow quantities with no modelling.
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_production_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Production')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Production-11')
        call self%add_model_field('Production-22')
        call self%add_model_field('Production-33')
        call self%add_model_field('Production-12')
        call self%add_model_field('Production-13')
        call self%add_model_field('Production-23')
        call self%add_model_field('Production-Trace')



    end subroutine init
    !***************************************************************************************






    !>  
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_production_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density_reynolds_11, density_reynolds_22, density_reynolds_33, &
            density_reynolds_12, density_reynolds_13, density_reynolds_23, &
            grad1_vel1, grad2_vel1, grad3_vel1,                             &
            grad1_vel2, grad2_vel2, grad3_vel2,                             &
            grad1_vel3, grad2_vel3, grad3_vel3,                             &
            production_11, production_22, production_33,                    &
            production_12, production_13, production_23,                    &
            production_trace,                                               &
            density, mu, nu, density_nutilde, mu_t, nutilde, lamda_t,   &
            chi, f_v1, k_t

        !
        ! Interpolate solution to quadrature nodes
        !
        
        density = worker%get_field('Density',    'value')
        density_reynolds_11 = worker%get_field('Reynolds-11',    'value')*density
        density_reynolds_22 = worker%get_field('Reynolds-22',    'value')*density
        density_reynolds_33 = worker%get_field('Reynolds-33',    'value')*density
        density_reynolds_12 = worker%get_field('Reynolds-12',    'value')*density
        density_reynolds_13 = worker%get_field('Reynolds-13',    'value')*density
        density_reynolds_23 = worker%get_field('Reynolds-23',    'value')*density
        !density_reynolds_11 = worker%get_field('Density * Reynolds-11',    'value')
        !density_reynolds_22 = worker%get_field('Density * Reynolds-22',    'value')
        !density_reynolds_33 = worker%get_field('Density * Reynolds-33',    'value')
        !density_reynolds_12 = worker%get_field('Density * Reynolds-12',    'value')
        !density_reynolds_13 = worker%get_field('Density * Reynolds-13',    'value')
        !density_reynolds_23 = worker%get_field('Density * Reynolds-23',    'value')


        grad1_vel1  = worker%get_field('Velocity 1 - Gradient 1',   'value')
        grad2_vel1  = worker%get_field('Velocity 1 - Gradient 2',   'value')
        grad3_vel1  = worker%get_field('Velocity 1 - Gradient 3',   'value')

        grad1_vel2  = worker%get_field('Velocity 2 - Gradient 1',   'value')
        grad2_vel2  = worker%get_field('Velocity 2 - Gradient 2',   'value')
        grad3_vel2  = worker%get_field('Velocity 2 - Gradient 3',   'value')

        grad1_vel3  = worker%get_field('Velocity 3 - Gradient 1',   'value')
        grad2_vel3  = worker%get_field('Velocity 3 - Gradient 2',   'value')
        grad3_vel3  = worker%get_field('Velocity 3 - Gradient 3',   'value')


        production_11 = -TWO*(density_reynolds_11*grad1_vel1 + density_reynolds_12*grad2_vel1 + density_reynolds_13*grad3_vel1)
        production_22 = -TWO*(density_reynolds_12*grad1_vel2 + density_reynolds_22*grad2_vel2 + density_reynolds_23*grad3_vel2)
        production_33 = -TWO*(density_reynolds_13*grad1_vel3 + density_reynolds_23*grad2_vel3 + density_reynolds_33*grad3_vel3)

        production_12 = -(density_reynolds_11*grad1_vel2 + density_reynolds_12*grad2_vel2 + density_reynolds_13*grad3_vel2 + &
                            density_reynolds_12*grad1_vel1 + density_reynolds_22*grad2_vel1 + density_reynolds_23*grad3_vel1)

        production_13 = -(density_reynolds_11*grad1_vel3 + density_reynolds_12*grad2_vel3 + density_reynolds_13*grad3_vel3 + &
                            density_reynolds_13*grad1_vel1 + density_reynolds_23*grad2_vel1 + density_reynolds_33*grad3_vel1)

        production_23 = -(density_reynolds_12*grad1_vel3 + density_reynolds_22*grad2_vel3 + density_reynolds_23*grad3_vel3 + &
                            density_reynolds_13*grad1_vel2 + density_reynolds_23*grad2_vel2 + density_reynolds_33*grad3_vel2)

        production_trace = production_11 + production_22 + production_33
        call worker%store_model_field('Production-11', 'value', production_11)
        call worker%store_model_field('Production-22', 'value', production_22)
        call worker%store_model_field('Production-33', 'value', production_33)
        call worker%store_model_field('Production-12', 'value', production_12)
        call worker%store_model_field('Production-13', 'value', production_13)
        call worker%store_model_field('Production-23', 'value', production_23)
        call worker%store_model_field('Production-Trace', 'value', production_trace)




    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_production
