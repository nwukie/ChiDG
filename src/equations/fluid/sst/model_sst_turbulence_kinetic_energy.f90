!>
!! TKE k computed on its own because it is needed for the SST ideal gas model.
!!
!! @author  Eric M. Wolf
!! @date    10/23/2018 
!!
!--------------------------------------------------------------------------------
module model_sst_turbulence_kinetic_energy
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, FOUR, PI
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_sst
    use mod_fluid

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: sst_turbulence_kinetic_energy_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_turbulence_kinetic_energy_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric M. Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(sst_turbulence_kinetic_energy_t), intent(inout)   :: self

        call self%set_name('SST Turbulence Kinetic Energy')
        call self%set_dependency('f(Q-)')


        ! k
 
        call self%add_model_field('k')

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
        class(sst_turbulence_kinetic_energy_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density, density_k, k, k_smooth

        real(rk) :: b, c, k_infty
        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
               
        density_k = worker%get_field('Density * k',    'value')
        k = density_k/density

        !b = 100.0_rk
        !c = HALF - atan(b)/PI
        !k_smooth = k*(atan(b*k)/PI+HALF)+c
        !print *, 'k_smooth - v1'
        !print *, k_smooth(:)%x_ad_
        !print *, 'k max arg'
        !print *, k(:)%x_ad_
        !k_smooth = smax(k, ZERO*k)
        !k_infty = (sst_tu_infty*69.0_rk)**TWO
        k_smooth = k
        k_smooth = k*sin_ramp(k, ZERO, 10.0_rk*sst_k_infty)
        !;print *, 'k_smooth - v2'
        !;print *, k_smooth(:)%x_ad_

        call worker%store_model_field('k', 'value', k_smooth)



    end subroutine compute
    !***************************************************************************************




end module model_sst_turbulence_kinetic_energy
