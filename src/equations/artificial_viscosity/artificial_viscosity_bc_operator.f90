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
            diffusion_x,    diffusion_y,    diffusion_z,    &
            deps_dx,        deps_dy,        deps_dz,        &
            flux_x,         flux_y,         flux_z,         &
            integrand


        integer(ik) :: order
        real(rk)    :: C1, C2, h(3),  &
                       eta_x, eta_y, eta_z

        real(rk),   allocatable, dimension(:)   :: &
            normx, normy, normz


        !
        ! Interpolate solution to quadrature nodes
        !
        deps_dx = worker%get_primary_field_face('Artificial Viscosity', 'grad1+lift', 'boundary')
        deps_dy = worker%get_primary_field_face('Artificial Viscosity', 'grad2+lift', 'boundary')
        deps_dz = worker%get_primary_field_face('Artificial Viscosity', 'grad3+lift', 'boundary')

        

        !
        ! Get model field for Maximum Wave Speed
        !
        lamda = worker%get_model_field_face('Maximum Wave Speed', 'value', 'boundary')



        !
        ! Get normal vector
        !
        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


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

        eta_x = C2*(h(1)**TWO)
        eta_y = C2*(h(2)**TWO)
        eta_z = C2*(h(3)**TWO)



        !
        ! Compute diffusion coefficient tensor diagonal
        !
        diffusion_x = eta_x/tau
        diffusion_y = eta_y/tau
        diffusion_z = eta_z/tau





        !================================
        !             FLUX
        !================================
        !flux_x = -(diffusion_x * deps_dx)
        !flux_y = -(diffusion_y * deps_dy)
        !flux_z = -(diffusion_z * deps_dz)
        flux_x = -0.20_rk*deps_dx
        flux_y = -0.20_rk*deps_dy
        flux_z = -0.20_rk*deps_dz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Artificial Viscosity',integrand)


    end subroutine compute
    !*********************************************************************************************












end module artificial_viscosity_bc_operator
