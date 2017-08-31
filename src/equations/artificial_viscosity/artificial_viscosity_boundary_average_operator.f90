module artificial_viscosity_boundary_average_operator
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
    !!  @date   1/29/2017
    !!
    !-------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: artificial_viscosity_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_boundary_average_operator_t
    !*******************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   01/31/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Artificial Viscosity Boundary Average Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Operator')

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
        class(artificial_viscosity_boundary_average_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                                       intent(inout)   :: worker
        class(properties_t),                                        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) ::         &
            lamda_m,        lamda_p,                        &
            tau_m,          tau_p,                          &
            diffusion_1_m,  diffusion_2_m,  diffusion_3_m,  &
            diffusion_1_p,  diffusion_2_p,  diffusion_3_p,  &
            grad1_eps_m,      grad2_eps_m,      grad3_eps_m,      &
            grad1_eps_p,      grad2_eps_p,      grad3_eps_p,      &
            flux_1_m,       flux_2_m,       flux_3_m,       &
            flux_1_p,       flux_2_p,       flux_3_p,       &
            flux_x,         flux_y,         flux_z,         &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz

        integer(ik) :: order_m, order_p
        real(rk)    :: C1, C2, h_m(3), h_p(3),      &
                       eta_1_m, eta_2_m, eta_3_m,   &
                       eta_1_p, eta_2_p, eta_3_p


        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_eps_m  = worker%get_field('Artificial Viscosity', 'grad1', 'face interior')
        grad1_eps_p  = worker%get_field('Artificial Viscosity', 'grad1', 'face exterior')

        grad2_eps_m  = worker%get_field('Artificial Viscosity', 'grad2', 'face interior')
        grad2_eps_p  = worker%get_field('Artificial Viscosity', 'grad2', 'face exterior')

        grad3_eps_m  = worker%get_field('Artificial Viscosity', 'grad3', 'face interior')
        grad3_eps_p  = worker%get_field('Artificial Viscosity', 'grad3', 'face exterior')



        !
        ! Get model field for Maximum Wave Speed
        !
        lamda_m = worker%get_field('Maximum Wave Speed', 'value', 'face interior')
        lamda_p = worker%get_field('Maximum Wave Speed', 'value', 'face exterior')




        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Get element resolution information
        !
        h_m     = worker%element_size(  'interior')
        order_m = worker%solution_order('interior')

        h_p     = worker%element_size(  'exterior')
        order_p = worker%solution_order('exterior')




        !
        ! Compute time constant and conductivity
        !
        C1 = THREE
        C2 = FIVE

        tau_m = minval(h_m)/(C1 * real(order_m+1,rk) * lamda_m)
        tau_p = minval(h_p)/(C1 * real(order_p+1,rk) * lamda_p)

        eta_1_m = C2*(h_m(1)**TWO)
        eta_2_m = C2*(h_m(2)**TWO)
        eta_3_m = C2*(h_m(3)**TWO)

        eta_1_p = C2*(h_p(1)**TWO)
        eta_2_p = C2*(h_p(2)**TWO)
        eta_3_p = C2*(h_p(3)**TWO)



        !
        ! Compute diffusion tensor diagonal
        !
        diffusion_1_m = eta_1_m/tau_m
        diffusion_2_m = eta_2_m/tau_m
        diffusion_3_m = eta_3_m/tau_m

        diffusion_1_p = eta_1_p/tau_p
        diffusion_2_p = eta_2_p/tau_p
        diffusion_3_p = eta_3_p/tau_p





        !================================
        !             FLUX
        !================================
        flux_1_m = -diffusion_1_m * grad1_eps_m
        flux_2_m = -diffusion_2_m * grad2_eps_m
        flux_3_m = -diffusion_3_m * grad3_eps_m

        flux_1_p = -diffusion_1_p * grad1_eps_p
        flux_2_p = -diffusion_2_p * grad2_eps_p
        flux_3_p = -diffusion_3_p * grad3_eps_p

        call worker%integrate_boundary_average('Artificial Viscosity','Diffusion',  &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)



    end subroutine compute
    !********************************************************************************************












end module artificial_viscosity_boundary_average_operator
