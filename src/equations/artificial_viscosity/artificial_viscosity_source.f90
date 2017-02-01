module artificial_viscosity_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: HALF, ONE, TWO, THREE, FIVE, THIRD
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
    !-----------------------------------------------------------------------------------------
    type, extends(operator_t), public :: artificial_viscosity_source_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_source_t
    !*****************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_source_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Artificial Viscosity Source')

        !
        ! Set operator type
        !
        call self%set_operator_type('Volume Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Artificial Viscosity')


    end subroutine init
    !*****************************************************************************************



    !>  Boundary Flux routine for Artificial Viscosity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(artificial_viscosity_source_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) ::         &
            eps,            lamda,          tau,            &
            diffusion_x,    diffusion_y,    diffusion_z,    &
            flux_x,         flux_y,         flux_z,         &
            integrand,      source,         sensor


        integer(ik) :: order
        real(rk)    :: C1, C2, h(3),  &
                       eta_x, eta_y, eta_z, hbar


        !
        ! Interpolate solution to quadrature nodes
        !
        eps = worker%get_primary_field_element('Artificial Viscosity', 'value')


        !
        ! Get model field for Maximum Wave Speed
        !
        lamda  = worker%get_model_field_element('Maximum Wave Speed', 'value')
        sensor = worker%get_model_field_element('Artificial Viscosity Sensor', 'value')



        !
        ! Get element resolution information
        !
        h     = worker%element_size(  'interior')
        order = worker%solution_order('interior')

        hbar = THIRD*(h(1) + h(2) + h(3))




        !
        ! Compute time constant and conductivity
        !
        C1 = THREE
        C2 = FIVE

        tau = minval(h)/(C1 * real(order+1,rk) * lamda)



        !================================
        !             FLUX
        !================================
        source = (ONE/tau)*( (hbar/real(order+1,rk))*lamda*sensor - eps )

        call worker%integrate_volume('Artificial Viscosity',source)


    end subroutine compute
    !******************************************************************************************












end module artificial_viscosity_source
