module spalart_allmaras_volume_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_volume_advection_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_volume_advection_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_volume_advection_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("Spalart-Allmaras Volume Advection Operator")

        ! Set operator type
        call self%set_operator_type("Volume Advective Operator")

        ! Set operator equations
        call self%add_primary_field("Density * NuTilde")

    end subroutine init
    !***********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_volume_advection_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                                intent(inout)   :: worker
        class(properties_t),                                 intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            rho, rhou, rhov, rhow, rho_nutilde,     &
            invrho, u, v, w, flux_x, flux_y, flux_z

        real(rk),   dimension(:), allocatable   ::  &
            normx, normy, normz, unormx, unormy, unormz


        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_primary_field_element('Density',           'value')
        rhou        = worker%get_primary_field_element('Momentum-1',        'value')
        rhov        = worker%get_primary_field_element('Momentum-2',        'value')
        rhow        = worker%get_primary_field_element('Momentum-3',        'value')
        rho_nutilde = worker%get_primary_field_element('Density * NuTilde', 'value')



        !
        ! Compute velocities
        !
        invrho = ONE/rho
        u = rhou*invrho
        v = rhov*invrho
        w = rhow*invrho


        !
        ! Get normal vector
        !
        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)


        !
        ! Compute average flux and field difference.
        ! 
        flux_x = u*rho_nutilde
        flux_y = v*rho_nutilde
        flux_z = w*rho_nutilde


        !
        ! Integrate flux
        !
        call worker%integrate_volume('Density * NuTilde',flux_x,flux_y,flux_z)


    end subroutine compute
    !************************************************************************************************









end module spalart_allmaras_volume_advection
