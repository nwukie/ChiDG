module euler_bc_operator
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO
    use mod_fluid,          only: omega
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
    type, public, extends(operator_t)   :: euler_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

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
        class(euler_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop

        ! data at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,   &
            H_bc, p_bc, u_a, v_a, w_a,                          &
            flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3, r
            

        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        !density_bc = worker%get_primary_field_face('Density'   , 'value', 'boundary')
        !mom1_bc    = worker%get_primary_field_face('Momentum-1', 'value', 'boundary')
        !mom2_bc    = worker%get_primary_field_face('Momentum-2', 'value', 'boundary')
        !mom3_bc    = worker%get_primary_field_face('Momentum-3', 'value', 'boundary')
        !energy_bc  = worker%get_primary_field_face('Energy'    , 'value', 'boundary')

        density_bc = worker%get_field('Density'   , 'value', 'boundary')
        mom1_bc    = worker%get_field('Momentum-1', 'value', 'boundary')
        mom2_bc    = worker%get_field('Momentum-2', 'value', 'boundary')
        mom3_bc    = worker%get_field('Momentum-3', 'value', 'boundary')
        energy_bc  = worker%get_field('Energy'    , 'value', 'boundary')

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary') 
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc / r
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
        !u_a = worker%get_model_field_face('Advection Velocity-1', 'value', 'boundary')
        !v_a = worker%get_model_field_face('Advection Velocity-2', 'value', 'boundary')
        !w_a = worker%get_model_field_face('Advection Velocity-3', 'value', 'boundary')
        u_a = worker%get_field('Advection Velocity-1', 'value', 'boundary')
        v_a = worker%get_field('Advection Velocity-2', 'value', 'boundary')
        w_a = worker%get_field('Advection Velocity-3', 'value', 'boundary')

        !
        ! Get pressure
        !
        !p_bc = worker%get_model_field_face('Pressure','value','boundary')
        p_bc = worker%get_field('Pressure','value','boundary')


        !
        ! Compute boundary condition energy and enthalpy
        !
        H_bc = (energy_bc + p_bc)/density_bc


        !=================================================
        ! mass flux
        !=================================================
        flux_1 = (density_bc * u_a )
        flux_2 = (density_bc * v_a )
        flux_3 = (density_bc * w_a )

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Density',integrand)

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = (mom1_bc * u_a) + p_bc
        flux_2 = (mom1_bc * v_a)
        flux_3 = (mom1_bc * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Momentum-1',integrand)

        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = (mom2_bc * u_a)
        flux_2 = (mom2_bc * v_a) + p_bc
        flux_3 = (mom2_bc * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        call worker%integrate_boundary('Momentum-2',integrand)

        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = (mom3_bc * u_a)
        flux_2 = (mom3_bc * v_a)
        flux_3 = (mom3_bc * w_a) + p_bc

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! energy flux
        !=================================================
        flux_1 = (density_bc * H_bc * u_a)
        flux_2 = (density_bc * H_bc * v_a)  +  omega*r*p_bc
        flux_3 = (density_bc * H_bc * w_a)

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Energy',integrand)

    end subroutine compute
    !**********************************************************************************************























end module euler_bc_operator
