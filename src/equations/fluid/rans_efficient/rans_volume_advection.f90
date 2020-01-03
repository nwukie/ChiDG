module rans_volume_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
    use mod_spalart_allmaras,   only: SA_sigma, SA_c_n1
    use mod_rans_efficient

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: rans_volume_advection_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rans_volume_advection_t
    !******************************************************************************


contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_volume_advection_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('RANS Volume Advection')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rans_volume_advection_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w, invdensity, H, p, T,           &
            flux_1, flux_2, flux_3, r

        type(AD_D), allocatable, dimension(:)   ::  &
            density_nutilde

        ! Interpolate boundary condition state to face quadrature nodes
        density         = worker%get_field('Density',           'value', 'element')
        mom1            = worker%get_field('Momentum-1',        'value', 'element')
        mom2            = worker%get_field('Momentum-2',        'value', 'element')
        mom3            = worker%get_field('Momentum-3',        'value', 'element')
        energy          = worker%get_field('Energy',            'value', 'element')

        if (turbulence_model == 'Spalart Allmaras') then
            density_nutilde = worker%get_field('Density * NuTilde', 'value', 'element')
        else
            density_nutilde = ZERO*density
        end if


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','volume')
            mom2 = mom2 / r
        end if


        ! Compute velocities
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        ! Compute pressure
        call compute_pressure_temperature(density,mom1,mom2,mom3,energy,p,T)


        ! Compute boundary condition energy and enthalpy
        H = (energy + p)*invdensity



        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = (density*u)
        flux_2 = (density*v)
        flux_3 = (density*w)
        call worker%integrate_volume_flux('Density','Advection',flux_1,flux_2,flux_3)
        

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (mom1*u) + p
        flux_2 = (mom1*v) 
        flux_3 = (mom1*w) 
        call worker%integrate_volume_flux('Momentum-1','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = (mom2*u) 
        flux_2 = (mom2*v) + p
        flux_3 = (mom2*w) 
        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if
        call worker%integrate_volume_flux('Momentum-2','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (mom3*u) 
        flux_2 = (mom3*v) 
        flux_3 = (mom3*w) + p
        call worker%integrate_volume_flux('Momentum-3','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! Energy flux
        !=================================================
        flux_1 = (density*H*u) 
        flux_2 = (density*H*v) 
        flux_3 = (density*H*w) 
        call worker%integrate_volume_flux('Energy','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! Spalart-Allmaras: Density * NuTilde
        !=================================================
        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                flux_1 = (density_nutilde*u)
                flux_2 = (density_nutilde*v)
                flux_3 = (density_nutilde*w)
                call worker%integrate_volume_flux('Density * NuTilde','Advection',flux_1,flux_2,flux_3)
            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select


    end subroutine compute


end module rans_volume_advection
