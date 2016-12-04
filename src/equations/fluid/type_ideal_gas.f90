module type_ideal_gas
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Pressure
    !!      - Temperature
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: ideal_gas_t

        real(rk)    :: gam = 1.4_rk     ! ratio of specific heats
        real(rk)    :: R   = 287.15_rk  ! ideal gas constant [J/(kg*K)]

    contains

        procedure   :: init
        procedure   :: compute

    end type ideal_gas_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(ideal_gas_t), intent(inout)   :: self

        call self%set_name('Ideal Gas')

        call self%add_model_field('Pressure')
        call self%add_model_field('Temperature')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(ideal_gas_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            rho, rhou, rhov, rhow, rhoE, P, T


        !
        ! Interpolate solution to quadrature nodes
        !
        rho  = worker%get_primary_field_general('Density',    'value')
        rhou = worker%get_primary_field_general('X-Momentum', 'value')
        rhov = worker%get_primary_field_general('Y-Momentum', 'value')
        rhow = worker%get_primary_field_general('Z-Momentum', 'value')
        rhoE = worker%get_primary_field_general('Energy',     'value')


        P = (self%gam-ONE)*(rhoE - HALF*( (rhou*rhou) + (rhov*rhov) + (rhow*rhow) )/rho )
        T = P/(rho*self%R)



        call worker%store_model_field('Pressure',    'value', P)
        call worker%store_model_field('Temperature', 'value', T)


    end subroutine compute
    !***************************************************************************************




end module type_ideal_gas
