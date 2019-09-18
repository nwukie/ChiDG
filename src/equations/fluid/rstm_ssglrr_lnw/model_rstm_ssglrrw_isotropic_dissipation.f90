module model_rstm_ssglrrw_isotropic_dissipation
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
    type, extends(model_t)  :: rstm_ssglrrw_isotropic_dissipation_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_isotropic_dissipation_t
    !***************************************************************************************





contains




    !>  Model to compute the isotropic_dissipation term in the Reynolds Stress transport equation.
    !!  Note that this term is computed exactly from mean flow quantities with no modelling.
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_isotropic_dissipation_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Isotropic Dissipation')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Turbulence Dissipation-11')
        call self%add_model_field('Turbulence Dissipation-22')
        call self%add_model_field('Turbulence Dissipation-33')
        call self%add_model_field('Turbulence Dissipation-12')
        call self%add_model_field('Turbulence Dissipation-13')
        call self%add_model_field('Turbulence Dissipation-23')



    end subroutine init
    !***************************************************************************************






    !>  
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_isotropic_dissipation_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
                density, epsilon_t,                             &
                dissipation_rate_11, dissipation_rate_22, dissipation_rate_33, &
                dissipation_rate_12, dissipation_rate_13, dissipation_rate_23
                

        density     = worker%get_field('Density',    'value')
        epsilon_t   = worker%get_field('Turbulence Isotropic Dissipation Rate', 'value')
        dissipation_rate_11 = (TWO/THREE)*density*epsilon_t
        dissipation_rate_22 = (TWO/THREE)*density*epsilon_t
        dissipation_rate_33 = (TWO/THREE)*density*epsilon_t
        dissipation_rate_12 = ZERO*density
        dissipation_rate_13 = ZERO*density
        dissipation_rate_23 = ZERO*density

        call worker%store_model_field('Turbulence Dissipation-11', 'value', dissipation_rate_11)
        call worker%store_model_field('Turbulence Dissipation-22', 'value', dissipation_rate_22)
        call worker%store_model_field('Turbulence Dissipation-33', 'value', dissipation_rate_33)
        call worker%store_model_field('Turbulence Dissipation-12', 'value', dissipation_rate_12)
        call worker%store_model_field('Turbulence Dissipation-13', 'value', dissipation_rate_13)
        call worker%store_model_field('Turbulence Dissipation-23', 'value', dissipation_rate_23)




    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_isotropic_dissipation
