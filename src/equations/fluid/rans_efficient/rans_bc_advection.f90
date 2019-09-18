module rans_bc_advection
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: HALF, ONE, TWO
    use mod_spalart_allmaras,   only: SA_sigma, SA_c_n1
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mod_rans_efficient
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: rans_bc_advection_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rans_bc_advection_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_bc_advection_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('RANS BC Advection')

        ! Set operator type
        call self%set_operator_type('BC Advective Flux')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )


    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rans_bc_advection_t),  intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w, invdensity, H, p, T,           &
            flux_1, flux_2, flux_3

        type(AD_D), allocatable, dimension(:)   ::  &
            density_nutilde, nutilde

        real(rk),   allocatable, dimension(:)   ::  r


        ! Interpolate boundary condition state to face quadrature nodes
        density         = worker%get_field('Density',           'value', 'boundary')
        mom1            = worker%get_field('Momentum-1',        'value', 'boundary')
        mom2            = worker%get_field('Momentum-2',        'value', 'boundary')
        mom3            = worker%get_field('Momentum-3',        'value', 'boundary')
        energy          = worker%get_field('Energy',            'value', 'boundary')

        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                density_nutilde = worker%get_field('Density * NuTilde', 'value', 'boundary')
            case('none')
                density_nutilde = ZERO*density
            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select




        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2 = mom2 / r
        end if


        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        !
        ! Compute pressure
        !
        call compute_pressure_temperature(density,mom1,mom2,mom3,energy,p,T)



        !
        ! Compute boundary condition energy and enthalpy
        !
        H = (energy + p)*invdensity




        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = (density*u)
        flux_2 = (density*v)
        flux_3 = (density*w)
        call worker%integrate_boundary_condition('Density','Advection',flux_1,flux_2,flux_3)
        

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (mom1*u) + p
        flux_2 = (mom1*v) 
        flux_3 = (mom1*w) 
        call worker%integrate_boundary_condition('Momentum-1','Advection',flux_1,flux_2,flux_3)


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
        call worker%integrate_boundary_condition('Momentum-2','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (mom3*u) 
        flux_2 = (mom3*v) 
        flux_3 = (mom3*w) + p
        call worker%integrate_boundary_condition('Momentum-3','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! Energy flux
        !=================================================
        flux_1 = (density*H*u) 
        flux_2 = (density*H*v) 
        flux_3 = (density*H*w) 
        call worker%integrate_boundary_condition('Energy','Advection',flux_1,flux_2,flux_3)


        !=================================================
        ! Spalart-Allmaras: Density * NuTilde
        !=================================================
        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                flux_1 = density_nutilde*u
                flux_2 = density_nutilde*v
                flux_3 = density_nutilde*w
                call worker%integrate_boundary_condition('Density * NuTilde','Advection',flux_1,flux_2,flux_3)
            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select





    end subroutine compute
    !**********************************************************************************************









end module rans_bc_advection
