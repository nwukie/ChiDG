module artificial_viscosity_boundary_average_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF
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
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: artificial_viscosity_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Artificial Viscosity Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Artificial Viscosity")
        call self%add_primary_field("X-Momentum")
        call self%add_primary_field("Y-Momentum")
        call self%add_primary_field("Z-Momentum")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(artificial_viscosity_boundary_average_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                                       intent(inout)   :: worker
        class(properties_t),                                        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            eps_m,      eps_p,                      &
            deps_dx_m,  deps_dy_m,  deps_dz_m,      &
            deps_dx_p,  deps_dy_p,  deps_dz_p,      &
            flux_x_m,   flux_y_m,   flux_z_m,       &
            flux_x_p,   flux_y_p,   flux_z_p,       &
            flux_x,     flux_y,     flux_z,         &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz



        !
        ! Interpolate solution to quadrature nodes
        !
        eps_m      = worker%get_primary_field_face('Artificial Viscosity', 'value',    'face interior')
        eps_p      = worker%get_primary_field_face('Artificial Viscosity', 'value',    'face exterior')

        deps_dx_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddx+lift', 'face interior')
        deps_dx_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddx+lift', 'face exterior')

        deps_dy_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddy+lift', 'face interior')
        deps_dy_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddy+lift', 'face exterior')

        deps_dz_m  = worker%get_primary_field_face('Artificial Viscosity', 'ddz+lift', 'face interior')
        deps_dz_p  = worker%get_primary_field_face('Artificial Viscosity', 'ddz+lift', 'face exterior')


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)




        !
        ! Compute diffusion coefficient
        !
        eta = 
        tau = 



        !================================
        !             FLUX
        !================================
        flux_x_m = rhou_m*H_m
        flux_y_m = rhov_m*H_m
        flux_z_m = rhow_m*H_m

        flux_x_p = rhou_p*H_p
        flux_y_p = rhov_p*H_p
        flux_z_p = rhow_p*H_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Artificial Viscosity',integrand)


    end subroutine compute
    !*********************************************************************************************************












end module artificial_viscosity_boundary_average_operator
