module RANS_bc_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,THREE,HALF
    use mod_fluid,              only: omega, gam

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
        class(RANS_bc_advection_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w,                                &
            enthalpy, pressure,                     &
            density_nutilde,                        &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_field('Density'   ,        'value', 'boundary')
        mom1            = worker%get_field('Momentum-1',        'value', 'boundary')
        mom2            = worker%get_field('Momentum-2',        'value', 'boundary')
        mom3            = worker%get_field('Momentum-3',        'value', 'boundary')
        energy          = worker%get_field('Energy'    ,        'value', 'boundary')
        density_nutilde = worker%get_field('Density * NuTilde', 'value', 'boundary')


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
        ! Get fluid advection velocity
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity


        !
        ! Compute Pressure, Enthalpy
        !
        pressure = get_field('Pressure', 'value', 'boundary')
        enthalpy = (energy + pressure)/density



        !=================================================
        ! mass flux
        !=================================================
        flux_1 = (density * u)
        flux_2 = (density * v)
        flux_3 = (density * w)

        call worker%integrate_boundary_condition('Density','Advection', flux_1, flux_2, flux_3)


        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (mom1 * u)  +  pressure
        flux_2 = (mom1 * v)
        flux_3 = (mom1 * w)

        call worker%integrate_boundary_condition('Momentum-1','Advection', flux_1, flux_2, flux_3)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = (mom2 * u)
        flux_2 = (mom2 * v)  +  pressure
        flux_3 = (mom2 * w)

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

        call worker%integrate_boundary_condition('Momentum-2','Advection', flux_1, flux_2, flux_3)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (mom3 * u)
        flux_2 = (mom3 * v)
        flux_3 = (mom3 * w)  +  pressure

        call worker%integrate_boundary_condition('Momentum-3','Advection', flux_1, flux_2, flux_3)


        !=================================================
        ! energy flux
        !=================================================
        flux_1 = (density * enthalpy * u)
        flux_2 = (density * enthalpy * v)
        flux_3 = (density * enthalpy * w)
        !flux_2 = (density * enthalpy * v)  +  omega*r*pressure

        call worker%integrate_boundary_condition('Energy','Advection', flux_1, flux_2, flux_3)


        !=================================================
        ! turbulence flux
        !=================================================
        flux_1 = (density_nutilde * u)
        flux_2 = (density_nutilde * v)
        flux_3 = (density_nutilde * w)

        call worker%integrate_boundary_condition('Density * NuTilde','Advection', flux_1, flux_2, flux_3)

    

    end subroutine compute
    !*******************************************************************************




end module RANS_bc_advection
