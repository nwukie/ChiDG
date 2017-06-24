module RANS_bc_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,THREE,HALF
    use mod_fluid,              only: omega

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: RANS_bc_advection_t

        real(rk)    :: gam = 1.4_rk
        real(rk)    :: R   = 287.15_rk
        real(rk)    :: Cp  = 1003.0_rk
        real(rk)    :: Pr  = 0.72_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type RANS_bc_advection_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(RANS_bc_advection_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("RANS BC Advection")

        ! Set operator type
        call self%set_operator_type("BC Advective Flux")

        ! Set operator equations
        call self%add_primary_field('Density'          )
        call self%add_primary_field('Momentum-1'       )
        call self%add_primary_field('Momentum-2'       )
        call self%add_primary_field('Momentum-3'       )
        call self%add_primary_field('Energy'           )
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(RANS_bc_advection_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u_a, v_a, w_a,                          &
            enthalpy, pressure,                     &
            density_nutilde,                        &
            flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   :: r, norm_1, norm_2, norm_3


        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_primary_field_face('Density'   ,        'value','boundary')
        mom1            = worker%get_primary_field_face('Momentum-1',        'value','boundary')
        mom2            = worker%get_primary_field_face('Momentum-2',        'value','boundary')
        mom3            = worker%get_primary_field_face('Momentum-3',        'value','boundary')
        energy          = worker%get_primary_field_face('Energy'    ,        'value','boundary')
        density_nutilde = worker%get_primary_field_face('Density * NuTilde', 'value','boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2       = mom2 / r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        !
        ! Get normal vector
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)




        !
        ! Get fluid advection velocity
        !
        u_a = worker%get_model_field_face('Advection Velocity-1', 'value','boundary')
        v_a = worker%get_model_field_face('Advection Velocity-2', 'value','boundary')
        w_a = worker%get_model_field_face('Advection Velocity-3', 'value','boundary')


        !
        ! Compute Pressure, Enthalpy
        !
        pressure = (self%gam-ONE)*(energy - HALF*( (mom1*mom1) + (mom2*mom2) + (mom3*mom3) )/density )
        enthalpy = (energy + pressure)/density



        !=================================================
        ! mass flux
        !=================================================
        flux_1 = (density * u_a)
        flux_2 = (density * v_a)
        flux_3 = (density * w_a)


        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Density',integrand)


        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (mom1 * u_a)  +  pressure
        flux_2 = (mom1 * v_a)
        flux_3 = (mom1 * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Momentum-1',integrand)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = (mom2 * u_a)
        flux_2 = (mom2 * v_a)  +  pressure
        flux_3 = (mom2 * w_a)

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Momentum-2',integrand)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (mom3 * u_a)
        flux_2 = (mom3 * v_a)
        flux_3 = (mom3 * w_a)  +  pressure

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Momentum-3',integrand)


        !=================================================
        ! energy flux
        !=================================================
        flux_1 = (density * enthalpy * u_a)
        flux_2 = (density * enthalpy * v_a)  +  omega*r*pressure
        flux_3 = (density * enthalpy * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Energy',integrand)


        !=================================================
        ! turbulence flux
        !=================================================
        flux_1 = (density_nutilde * u_a)
        flux_2 = (density_nutilde * v_a)
        flux_3 = (density_nutilde * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
        call worker%integrate_boundary('Density * NuTilde',integrand)

    

    end subroutine compute
    !*******************************************************************************




end module RANS_bc_advection
