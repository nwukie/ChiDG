module artificial_viscosity_bc_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: HALF, ONE, TWO, THREE, FIVE
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private



    !> Implementation of the artificial viscosity boundary averate flux.
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !-------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: artificial_viscosity_bc_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_bc_operator_t
    !********************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !--------------------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Artificial Viscosity BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Artificial Viscosity')


    end subroutine init
    !*******************************************************************************************



    !>  Boundary Flux routine for Artificial Viscosity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(artificial_viscosity_bc_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) ::         &
            lamda,  tau,                                    &
            diffusion_1,    diffusion_2,    diffusion_3,    &
            deps_dx,        deps_dy,        deps_dz,        &
            flux_1,         flux_2,         flux_3,         &
            integrand


        integer(ik) :: order
        real(rk)    :: C1, C2, h(3),  &
                       eta_1, eta_2, eta_3


        !
        ! Interpolate solution to quadrature nodes
        !
        deps_dx = worker%get_field('Artificial Viscosity', 'grad1', 'boundary')
        deps_dy = worker%get_field('Artificial Viscosity', 'grad2', 'boundary')
        deps_dz = worker%get_field('Artificial Viscosity', 'grad3', 'boundary')

        

        !
        ! Get model field for Maximum Wave Speed
        !
        lamda = worker%get_field('Maximum Wave Speed', 'value', 'boundary')


        !
        ! Get element resolution information
        !
        h     = worker%element_size('interior')
        order = worker%solution_order('interior')


        !
        ! Compute time constant and conductivity
        !
        C1 = THREE
        C2 = FIVE


        tau = minval(h)/(C1 * real(order+1,rk) * lamda)
        eta_1 = C2*(h(1)**TWO)
        eta_2 = C2*(h(2)**TWO)
        eta_3 = C2*(h(3)**TWO)



        !
        ! Compute diffusion coefficient tensor diagonal
        !
        diffusion_1 = eta_1/tau
        diffusion_2 = eta_2/tau
        diffusion_3 = eta_3/tau





        !================================
        !             FLUX
        !================================
        flux_1 = -(diffusion_1 * deps_dx)
        flux_2 = -(diffusion_2 * deps_dy)
        flux_3 = -(diffusion_3 * deps_dz)

        call worker%integrate_boundary_condition('Artificial Viscosity','Diffusion',flux_1,flux_2,flux_3)


    end subroutine compute
    !*********************************************************************************************












end module artificial_viscosity_bc_operator
