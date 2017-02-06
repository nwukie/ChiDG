module type_artificial_viscosity_resolution_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO, THREE, NFACES, PI
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing a sensor quantity for the artificial viscosity 
    !!  equation. Sensor quantity is computed based on the ratio of energy contained
    !!  in higher-order modes to the lower-order modes.
    !!
    !!  Model Fields:
    !!      - Artificial Viscosity Sensor
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/02/2017
    !!
    !------------------------------------------------------------------------------
    type, extends(model_t)  :: artificial_viscosity_resolution_sensor_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_resolution_sensor_t
    !******************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !------------------------------------------------------------------------------
    subroutine init(self)   
        class(artificial_viscosity_resolution_sensor_t), intent(inout)   :: self

        call self%set_name('Artificial Viscosity Resolution Sensor')

        call self%add_model_field('Artificial Viscosity Sensor')

    end subroutine init
    !******************************************************************************






    !>  Routine for computing the 'Artificial Viscosity Sensor'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(artificial_viscosity_resolution_sensor_t),    intent(in)      :: self
        type(chidg_worker_t),                               intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::  &
            turb_total,     turb_low,   turb_gen,   &
            diff_squared,   total_squared,          &
            sensor

        type(AD_D)  :: diff_integral, total_integral, Fk, tmp

        integer(ik) :: order
        real(rk)    :: theta, psi, dpsi

        real(rk),   allocatable,    dimension(:)    ::  &
            weights, jinv


        !
        ! Get polynomial order
        !
        order = worker%solution_order('interior')
        

        !
        ! Interpolate solution to quadrature nodes
        !
        !turb_total = worker%get_primary_field_element('Density * NuTilde', 'value')
        turb_total = worker%get_primary_field_element('Density', 'value')




        if (order > 0) then

            !turb_low   = worker%get_primary_field_element('Density * NuTilde', 'value', Pmin=0,Pmax=(order-1))
            turb_low   = worker%get_primary_field_element('Density', 'value', Pmin=0,Pmax=(order-1))



            !
            ! Compute integrands
            !
            diff_squared = (turb_total - turb_low)**TWO
            total_squared = (turb_total)**TWO


            !
            ! Get quadrature weights
            !
            weights = worker%quadrature_weights('element')
            jinv    = worker%inverse_jacobian('element')

            !

            !
            ! Integrate
            !
            diff_integral  = sum(diff_squared  * weights * jinv)
            total_integral = sum(total_squared * weights * jinv)



            !
            ! Sensor
            !
            Fk = log10(diff_integral/total_integral)



            !
            ! Scale the sensor to be [0,maxval]
            !
            theta = 100000.0_rk
            psi  = -(3.0_rk + 4.25_rk*log10(real(order+1,rk)))
            dpsi = 0.5_rk



            turb_gen = worker%get_primary_field_general('Density * NuTilde', 'value')
            sensor = turb_gen
            sensor = ZERO
            if ( Fk <= (psi - dpsi) ) then
                sensor = ZERO
            else if ( Fk >= (psi + dpsi) ) then
                sensor = theta
            else if ( abs(Fk - psi) < dpsi ) then
                sensor = (HALF*theta)*(ONE + sin(PI*(Fk-psi)/(TWO*dpsi)))
            end if



        else

            !turb_gen = worker%get_primary_field_general('Density * NuTilde', 'value')
            turb_gen = worker%get_primary_field_general('Density', 'value')
            sensor = turb_gen
            sensor = ZERO

        end if



        call worker%store_model_field('Artificial Viscosity Sensor', 'value', sensor)


    end subroutine compute
    !******************************************************************************




end module type_artificial_viscosity_resolution_sensor
