!>
!! Description: This model computes the Reynolds stress tensor for the mean flow equations
!!              from the computed values of density*Rij.
!!
!! @author Eric M. Wolf
!! @date   01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_rstm_ssglrrw_reynolds_stress
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
    type, extends(model_t)  :: rstm_ssglrrw_reynolds_stress_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_reynolds_stress_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric Wolf 
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_reynolds_stress_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Reynolds Stress')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Reynolds-11')
        call self%add_model_field('Reynolds-22')
        call self%add_model_field('Reynolds-33')
        call self%add_model_field('Reynolds-12')
        call self%add_model_field('Reynolds-13')
        call self%add_model_field('Reynolds-23')



    end subroutine init
    !***************************************************************************************




    !>
    !! Description: Computes the Reynolds stress tensor
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_reynolds_stress_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density,   &
            density_reynolds_11, density_reynolds_22, density_reynolds_33, &
            density_reynolds_12, density_reynolds_13, density_reynolds_23

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        
        density_reynolds_11 = worker%get_field('Density * Reynolds-11',    'value')
        density_reynolds_22 = worker%get_field('Density * Reynolds-22',    'value')
        density_reynolds_33 = worker%get_field('Density * Reynolds-33',    'value')
        density_reynolds_12 = worker%get_field('Density * Reynolds-12',    'value')
        density_reynolds_13 = worker%get_field('Density * Reynolds-13',    'value')
        density_reynolds_23 = worker%get_field('Density * Reynolds-23',    'value')

        call worker%store_model_field('Reynolds-11', 'value', density_reynolds_11/density)
        call worker%store_model_field('Reynolds-22', 'value', density_reynolds_22/density)
        call worker%store_model_field('Reynolds-33', 'value', density_reynolds_33/density)
        call worker%store_model_field('Reynolds-12', 'value', density_reynolds_12/density)
        call worker%store_model_field('Reynolds-13', 'value', density_reynolds_13/density)
        call worker%store_model_field('Reynolds-23', 'value', density_reynolds_23/density)




    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_reynolds_stress
