module sst_bc_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,THREE,HALF
    use mod_fluid,              only: omega
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_bc_advection_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_bc_advection_operator_t
    !***********************************************************************************************

contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_bc_advection_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('SST BC Advection Operator')

        ! Set operator type
        call self%set_operator_type('BC Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy')
        

    end subroutine init
    !************************************************************************************************



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_bc_advection_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                            intent(inout)   :: worker
        class(properties_t),                             intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::          &
            density, mom1, mom2, mom3, density_omega, density_k,    &
            u, v, w, invdensity, flux_1, flux_2, flux_3

        real(rk),   allocatable,    dimension(:)    :: r

        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_field('Density',           'value', 'boundary')
        mom1            = worker%get_field('Momentum-1',        'value', 'boundary')
        mom2            = worker%get_field('Momentum-2',        'value', 'boundary')
        mom3            = worker%get_field('Momentum-3',        'value', 'boundary')
        density_omega   = worker%get_field('Density * Omega', 'value', 'boundary')
        density_k   = worker%get_field('Density * k', 'value', 'boundary')

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary') 
            mom2 = mom2 / r
        end if



        !
        ! Get fluid advection velocity
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity
    

        ! Density * Omega

        !
        ! Compute advection of spalart-allmaras working variable 
        ! 
        flux_1 = density_omega * u
        flux_2 = density_omega * v
        flux_3 = density_omega * w


        call worker%integrate_boundary_condition('Density * Omega','Advection',flux_1,flux_2,flux_3)

        ! Density * k

        !
        ! Compute advection of spalart-allmaras working variable 
        ! 
        flux_1 = density_k * u
        flux_2 = density_k * v
        flux_3 = density_k * w


        call worker%integrate_boundary_condition('Density * k','Advection',flux_1,flux_2,flux_3)



        !
        ! Momentum Turbulence Advection Flux
        !

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = (TWO/THREE)*density_k 
        flux_2 = ZERO*(TWO/THREE)*density_k 
        flux_3 = ZERO*(TWO/THREE)*density_k 


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_condition('Momentum-1','Advection',flux_1,flux_2,flux_3)

        flux_1 = ZERO*(TWO/THREE)*density_k 
        flux_2 = (TWO/THREE)*density_k 
        flux_3 = ZERO*(TWO/THREE)*density_k 


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_condition('Momentum-2','Advection',flux_1,flux_2,flux_3)

        flux_1 = ZERO*(TWO/THREE)*density_k 
        flux_2 = ZERO*(TWO/THREE)*density_k 
        flux_3 = (TWO/THREE)*density_k 


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_condition('Momentum-3','Advection',flux_1,flux_2,flux_3)


        flux_1 = (TWO/THREE)*density_k * u
        flux_2 = (TWO/THREE)*density_k * v
        flux_3 = (TWO/THREE)*density_k * w


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_condition('Energy','Advection',flux_1,flux_2,flux_3)





    end subroutine compute
    !************************************************************************************************









end module sst_bc_advection
