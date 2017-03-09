module spalart_allmaras_bc_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: omega
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
    type, extends(operator_t), public :: spalart_allmaras_bc_advection_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_bc_advection_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_bc_advection_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Spalart-Allmaras BC Advection Operator')

        ! Set operator type
        call self%set_operator_type('BC Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !************************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_bc_advection_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                            intent(inout)   :: worker
        class(properties_t),                             intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            density, mom1, mom2, mom3, density_nutilde,     &
            invdensity, u, v, w, flux_1, flux_2, flux_3, integrand

        real(rk),   dimension(:), allocatable   ::  &
            norm_1, norm_2, norm_3, unorm_1, unorm_2, unorm_3, r



        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_primary_field_face('Density',           'value', 'boundary')
        mom1            = worker%get_primary_field_face('Momentum-1',        'value', 'boundary')
        mom2            = worker%get_primary_field_face('Momentum-2',        'value', 'boundary')
        mom3            = worker%get_primary_field_face('Momentum-3',        'value', 'boundary')
        density_nutilde = worker%get_primary_field_face('Density * NuTilde', 'value', 'boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2 = mom2 / worker%coordinate('1','boundary')
        end if


        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity


        !
        ! Compute transport velocity
        !
        r = worker%coordinate('1','boundary') 
        v = v - omega*r

        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)


        !
        ! Compute advection of spalart-allmaras working variable 
        ! 
        flux_1 = u*density_nutilde
        flux_2 = v*density_nutilde
        flux_3 = w*density_nutilde


        !
        ! Integrate flux
        !
        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !************************************************************************************************









end module spalart_allmaras_bc_advection
