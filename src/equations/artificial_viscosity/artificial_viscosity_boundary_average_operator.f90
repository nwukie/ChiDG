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
            diffusion_x_m,  diffusion_y_m,  diffusion_z_m,  &
            diffusion_x_p,  diffusion_y_p,  diffusion_z_p,  &
            deps_dx_m,      deps_dy_m,      deps_dz_m,      &
            deps_dx_p,      deps_dy_p,      deps_dz_p,      &
            flux_x_m,       flux_y_m,       flux_z_m,       &
            flux_x_p,       flux_y_p,       flux_z_p,       &
            flux_x,         flux_y,         flux_z,         &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz

        integer(ik) :: order_m, order_p
        real(rk)    :: C1, C2, h_m(3), h_p(3),      &
                       eta_x_m, eta_y_m, eta_z_m,   &
                       eta_x_p, eta_y_p, eta_z_p


        !
        ! Interpolate solution to quadrature nodes
        !
        deps_dx_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddx+lift', 'face interior')
        deps_dx_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddx+lift', 'face exterior')

        deps_dy_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddy+lift', 'face interior')
        deps_dy_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddy+lift', 'face exterior')

        deps_dz_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddz+lift', 'face interior')
        deps_dz_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddz+lift', 'face exterior')



        !
        ! Get model field for Maximum Wave Speed
        !
        lamda_m = worker%get_model_field_face('Maximum Wave Speed', 'value', 'face interior')
        lamda_p = worker%get_model_field_face('Maximum Wave Speed', 'value', 'face exterior')




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

        eta_x_m = C2*(h_m(1)**TWO)
        eta_y_m = C2*(h_m(2)**TWO)
        eta_z_m = C2*(h_m(3)**TWO)

        eta_x_p = C2*(h_p(1)**TWO)
        eta_y_p = C2*(h_p(2)**TWO)
        eta_z_p = C2*(h_p(3)**TWO)



        !
        ! Compute diffusion tensor diagonal
        !
        diffusion_x_m = eta_x_m/tau_m
        diffusion_y_m = eta_y_m/tau_m
        diffusion_z_m = eta_z_m/tau_m

        diffusion_x_p = eta_x_p/tau_p
        diffusion_y_p = eta_y_p/tau_p
        diffusion_z_p = eta_z_p/tau_p





        !================================
        !             FLUX
        !================================
!        flux_x_m = -diffusion_x_m * deps_dx_m
!        flux_y_m = -diffusion_y_m * deps_dy_m
!        flux_z_m = -diffusion_z_m * deps_dz_m
!
!        flux_x_p = -diffusion_x_p * deps_dx_p
!        flux_y_p = -diffusion_y_p * deps_dy_p
!        flux_z_p = -diffusion_z_p * deps_dz_p
        flux_x_m = -0.20_rk*deps_dx_m
        flux_y_m = -0.20_rk*deps_dy_m
        flux_z_m = -0.20_rk*deps_dz_m

        flux_x_p = -0.20_rk*deps_dx_p
        flux_y_p = -0.20_rk*deps_dy_p
        flux_z_p = -0.20_rk*deps_dz_p


        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Artificial Viscosity',integrand)


    end subroutine compute
    !********************************************************************************************












end module artificial_viscosity_boundary_average_operator
