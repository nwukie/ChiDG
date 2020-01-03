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


        type(AD_D), dimension(:), allocatable   ::          &
            density, mom1, mom2, mom3, density_nutilde,     &
            u, v, w, invdensity, flux_1, flux_2, flux_3, r


        ! Interpolate solution to quadrature nodes
        density         = worker%get_field('Density',           'value', 'boundary')
        mom1            = worker%get_field('Momentum-1',        'value', 'boundary')
        mom2            = worker%get_field('Momentum-2',        'value', 'boundary')
        mom3            = worker%get_field('Momentum-3',        'value', 'boundary')
        density_nutilde = worker%get_field('Density * NuTilde', 'value', 'boundary')


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary') 
            mom2 = mom2 / r
        end if


        ! Get fluid advection velocity
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity
    

        ! Compute advection of spalart-allmaras working variable 
        flux_1 = density_nutilde * u
        flux_2 = density_nutilde * v
        flux_3 = density_nutilde * w


        call worker%integrate_boundary_condition('Density * NuTilde','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !************************************************************************************************









end module spalart_allmaras_bc_advection
