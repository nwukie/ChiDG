module fluid_viscous_bc_operator
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: fluid_viscous_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Fluid Viscous BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Flux')

        !
        ! Set operator equations
        !
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
        class(fluid_viscous_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w, invdensity,                    &
            grad1_T, grad2_T, grad3_T,              &
            k,       k_l,     k_t,                  &
            tau_11, tau_22, tau_33,                 &
            tau_12, tau_13, tau_23,                 &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   ::  r

        integer(ik) :: idomain_l, ielement_l, iface, group_ID


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        density = worker%get_field('Density'   , 'value', 'boundary')
        mom1    = worker%get_field('Momentum-1', 'value', 'boundary')
        mom2    = worker%get_field('Momentum-2', 'value', 'boundary')
        mom3    = worker%get_field('Momentum-3', 'value', 'boundary')
        energy  = worker%get_field('Energy'    , 'value', 'boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2 = mom2 / r
        end if


        !
        ! Get Model fields:
        !   Pressure
        !   Thermal Conductivity
        !
        k_l = worker%get_field('Laminar Thermal Conductivity',   'value', 'boundary')
        k_t = worker%get_field('Turbulent Thermal Conductivity', 'value', 'boundary')



        !
        ! Compute effective conductivity. Laminar + Turbulent.
        !
        k = k_l + k_t



        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        !
        ! get temperature gradient
        !
        grad1_T = worker%get_field('Temperature Gradient - 1', 'value', 'boundary')
        grad2_T = worker%get_field('Temperature Gradient - 2', 'value', 'boundary')
        grad3_T = worker%get_field('Temperature Gradient - 3', 'value', 'boundary')


        !
        ! get shear stress components
        !
        tau_11 = worker%get_field('Shear-11', 'value', 'boundary')
        tau_22 = worker%get_field('Shear-22', 'value', 'boundary')
        tau_33 = worker%get_field('Shear-33', 'value', 'boundary')

        tau_12 = worker%get_field('Shear-12', 'value', 'boundary')
        tau_13 = worker%get_field('Shear-13', 'value', 'boundary')
        tau_23 = worker%get_field('Shear-23', 'value', 'boundary')

        !=================================================
        ! Mass flux
        !=================================================
        

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = -tau_11
        flux_2 = -tau_12
        flux_3 = -tau_13

        call worker%integrate_boundary_condition('Momentum-1','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = -tau_12
        flux_2 = -tau_22
        flux_3 = -tau_23

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if

        call worker%integrate_boundary_condition('Momentum-2','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = -tau_13
        flux_2 = -tau_23
        flux_3 = -tau_33

        call worker%integrate_boundary_condition('Momentum-3','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! Energy flux
        !=================================================
        flux_1 = -k*grad1_T  -  (u*tau_11 + v*tau_12 + w*tau_13)
        flux_2 = -k*grad2_T  -  (u*tau_12 + v*tau_22 + w*tau_23)
        flux_3 = -k*grad3_T  -  (u*tau_13 + v*tau_23 + w*tau_33)

        call worker%integrate_boundary_condition('Energy','Diffusion',flux_1,flux_2,flux_3)

    end subroutine compute
    !**********************************************************************************************









end module fluid_viscous_bc_operator
