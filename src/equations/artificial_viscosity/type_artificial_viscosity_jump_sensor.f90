module type_artificial_viscosity_jump_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO, THREE, NFACES, PI
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing a sensor quantity for the artificial viscosity 
    !!  equation. Sensor quantity is computed based on the solution jump between 
    !!  elements.
    !!
    !!  Model Fields:
    !!      - Artificial Viscosity Sensor
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !------------------------------------------------------------------------------
    type, extends(model_t)  :: artificial_viscosity_jump_sensor_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_jump_sensor_t
    !******************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !------------------------------------------------------------------------------
    subroutine init(self)   
        class(artificial_viscosity_jump_sensor_t), intent(inout)   :: self

        call self%set_name('Artificial Viscosity Jump Sensor')

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
        class(artificial_viscosity_jump_sensor_t),  intent(in)      :: self
        type(chidg_worker_t),                       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::          &
            turb_m,         turb_p,                         &
            jump_x,         jump_y,         jump_z,         &
            integrand_x,    integrand_y,    integrand_z,    &
            integrand,      average,        sensor, turb_element

        type(AD_D)  :: integral

        integer(ik) :: original_face, iface, order
        real(rk)    :: total_area, face_area, theta, psi, dpsi

        real(rk),   allocatable,    dimension(:)    ::  &
            normx_m,  normy_m,  normz_m,                &
            normx_p,  normy_p,  normz_p,                &
            unormx_m, unormy_m, unormz_m,               &
            unormx_p, unormy_p, unormz_p,               &
            face_weights

        

        !
        ! Interpolate solution to quadrature nodes
        !
        !turb_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        turb_m = worker%get_field('Energy', 'value', 'face interior')


        !
        ! Get original face being worked on
        !
        original_face = worker%iface




        !
        ! Loop through faces, integrate jump in primary field
        !
        integral = turb_m(1)
        integral   = ZERO
        total_area = ZERO

        do iface = 1,NFACES

            !
            ! Update face for worker
            !
            call worker%set_face(iface)


            !
            ! Get normal vectors
            !
            normx_m = worker%normal(1)
            normy_m = worker%normal(2)
            normz_m = worker%normal(3)

            normx_p = -normx_m
            normy_p = -normy_m
            normz_p = -normz_m

            unormx_m = worker%unit_normal(1)
            unormy_m = worker%unit_normal(2)
            unormz_m = worker%unit_normal(3)

            unormx_p = -unormx_m
            unormy_p = -unormy_m
            unormz_p = -unormz_m


            !
            ! Get exterior field
            !
            !turb_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')
            turb_p = worker%get_field('Energy', 'value', 'face exterior')


            jump_x   = turb_p*unormx_p + turb_m*unormx_m
            jump_y   = turb_p*unormy_p + turb_m*unormy_m
            jump_z   = turb_p*unormz_p + turb_m*unormz_m
            average = HALF*(turb_p + turb_m)

            integrand_x = abs(jump_x/average)
            integrand_y = abs(jump_y/average)
            integrand_z = abs(jump_z/average)

            integrand = integrand_x*normx_m + integrand_y*normy_m + integrand_z*normz_m


            !
            ! Get quadrature weights
            !
            face_weights = worker%quadrature_weights('face')
            face_area    = worker%face_area()

            integral = integral + abs(sum(integrand * face_weights))

            total_area = total_area + face_area

        end do !iface


        !
        ! Normalize sensor by area
        !
        integral = integral/abs(total_area)


        !
        ! Reset face
        !
        call worker%set_face(original_face)


        !
        ! Get polynomial order
        !
        order = worker%solution_order('interior')


        !
        ! Scale the sensor to be [0,maxval]
        !
        theta = ONE
        !psi  = -(2.25_rk + THREE*log10(real(order+1,rk)))
        psi  = -(2.25_rk + THREE*log10(real(order+1,rk)))
        dpsi = HALF



        !turb_element = worker%get_primary_field_general('Density * NuTilde', 'value')
        turb_element = worker%get_field('Energy', 'value')
        sensor = turb_element
        sensor = ZERO

!        print*, 'integral: ', abs(integral%x_ad_)
!        print*, 'log10(integral): ', log10(abs(integral%x_ad_))

        print*, 'integral, p-dp, p+dp: ', integral%x_ad_, (psi-dpsi), (psi+dpsi)


        integral = log10(integral)

        if (abs(integral) < 1.e-12) then

            sensor = ZERO

        else

            !if ( log10(abs(integral)) <= (psi - dpsi) ) then
            !    sensor = ZERO
            !else if ( log10(abs(integral)) >= (psi + dpsi) ) then
            !    sensor = theta
            !else if ( abs(log10(abs(integral)) - psi) < dpsi ) then
            !    sensor = (HALF*theta)*(ONE + sin(PI*(log10(abs(integral))-psi)/(TWO*dpsi)))
            !end if
            if ( integral <= (psi - dpsi) ) then
                sensor = ZERO
            else if ( integral >= (psi + dpsi) ) then
                sensor = theta
            else if ( abs(integral - psi) < dpsi ) then
                sensor = (HALF*theta)*(ONE + sin(PI*(integral-psi)/(TWO*dpsi)))
            end if

        end if

        print*, 'sensor: ', sensor(1)%x_ad_

        call worker%store_model_field('Artificial Viscosity Sensor', 'value', sensor)


    end subroutine compute
    !******************************************************************************




end module type_artificial_viscosity_jump_sensor
