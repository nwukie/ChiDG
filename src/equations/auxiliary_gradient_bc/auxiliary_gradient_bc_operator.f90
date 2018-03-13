module auxiliary_gradient_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE, HALF
    use mod_fluid,          only: gam, Rgas
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: auxiliary_gradient_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type auxiliary_gradient_bc_operator_t
    !***************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine init(self)
        class(auxiliary_gradient_bc_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('Auxiliary Gradient BC Operator')

        ! Set operator type
        call self%set_operator_type('BC Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density_TEMP')
        call self%add_primary_field('Velocity-1_TEMP')
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !**************************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !---------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(auxiliary_gradient_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                  &
            grad1_p_bc,       grad2_p_bc,       grad3_p_bc,         &
            grad1_p,          grad2_p,          grad3_p,            &
            grad1_vel1_bc,    grad2_vel1_bc,    grad3_vel1_bc,      &
            grad1_vel1,       grad2_vel1,       grad3_vel1,         &
            grad1_density_bc, grad2_density_bc, grad3_density_bc,   &
            grad1_density,    grad2_density,    grad3_density,      &
            flux_1,    flux_2,    flux_3


        !-------------------------------------
        ! Auxiliary Density
        !-------------------------------------
        grad1_density_bc = worker%get_field('Density_TEMP', 'grad1', 'boundary')
        grad2_density_bc = worker%get_field('Density_TEMP', 'grad2', 'boundary')
        grad3_density_bc = worker%get_field('Density_TEMP', 'grad3', 'boundary')

        grad1_density    = worker%get_field('Density', 'grad1', 'boundary', override_lift=.true.)
        grad2_density    = worker%get_field('Density', 'grad2', 'boundary', override_lift=.true.)
        grad3_density    = worker%get_field('Density', 'grad3', 'boundary', override_lift=.true.)


        flux_1 = grad1_density_bc - grad1_density
        flux_2 = grad2_density_bc - grad2_density
        flux_3 = grad3_density_bc - grad3_density

        call worker%integrate_boundary_condition('Density_TEMP','Diffusion',flux_1,flux_2,flux_3)




        !-------------------------------------
        ! Auxiliary Velocity-1
        !-------------------------------------
        grad1_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad1', 'boundary')
        grad2_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad2', 'boundary')
        grad3_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad3', 'boundary')

        grad1_vel1    = worker%get_field('Velocity-1 Gradient-1', 'value', 'boundary')
        grad2_vel1    = worker%get_field('Velocity-1 Gradient-2', 'value', 'boundary')
        grad3_vel1    = worker%get_field('Velocity-1 Gradient-3', 'value', 'boundary')


        flux_1 = grad1_vel1_bc - grad1_vel1
        flux_2 = grad2_vel1_bc - grad2_vel1
        flux_3 = grad3_vel1_bc - grad3_vel1

        call worker%integrate_boundary_condition('Velocity-1_TEMP','Diffusion',flux_1,flux_2,flux_3)





        !-------------------------------------
        ! Auxiliary Pressure
        !-------------------------------------
        grad1_p_bc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary')
        grad2_p_bc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary')
        grad3_p_bc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary')

        grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
        grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
        grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')


        flux_1 = grad1_p_bc - grad1_p
        flux_2 = grad2_p_bc - grad2_p
        flux_3 = grad3_p_bc - grad3_p


        call worker%integrate_boundary_condition('Pressure_TEMP','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !***************************************************************************************










end module auxiliary_gradient_bc_operator
