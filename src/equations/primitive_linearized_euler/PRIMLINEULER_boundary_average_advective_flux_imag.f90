module PRIMLINEULER_boundary_average_advective_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use mod_primitive_linearized_euler
    implicit none
    private




    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: PRIMLINEULER_boundary_average_advective_flux_imag_t

    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_boundary_average_advective_flux_imag_t
    !*******************************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_boundary_average_advective_flux_imag_t), intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("PRIMLINEULER Boundary Average Imag")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density(imag)"   )
        call self%add_primary_field("Velocity-1(imag)")
        call self%add_primary_field("Velocity-2(imag)")
        call self%add_primary_field("Velocity-3(imag)")
        call self%add_primary_field("Pressure(imag)"  )

    end subroutine init
    !********************************************************************************











    !>   Boundary Flux routine for Euler
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_boundary_average_advective_flux_imag_t), intent(in)      :: self
        type(chidg_worker_t),                                       intent(inout)   :: worker
        class(properties_t),                                        intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,   rho_p,                         &
            u_m,     u_p,                           &
            v_m,     v_p,                           &
            w_m,     w_p,                           &
            p_m,     p_p,                           &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p


        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m = worker%interpolate('Density(imag)',    'value', 'face interior')
        rho_p = worker%interpolate('Density(imag)',    'value', 'face exterior')

        u_m   = worker%interpolate('Velocity-1(imag)', 'value', 'face interior')
        u_p   = worker%interpolate('Velocity-1(imag)', 'value', 'face exterior')

        v_m   = worker%interpolate('Velocity-2(imag)', 'value', 'face interior')
        v_p   = worker%interpolate('Velocity-2(imag)', 'value', 'face exterior')

        w_m   = worker%interpolate('Velocity-3(imag)', 'value', 'face interior')
        w_p   = worker%interpolate('Velocity-3(imag)', 'value', 'face exterior')

        p_m   = worker%interpolate('Pressure(imag)',   'value', 'face interior')
        p_p   = worker%interpolate('Pressure(imag)',   'value', 'face exterior')







        !================================
        !       MASS FLUX
        !================================
        flux_1_m = rho_1_rho * rho_m + &
                   rho_1_u   * u_m   + &
                   rho_1_v   * v_m   + &
                   rho_1_w   * w_m   + &
                   rho_1_p   * p_m
        flux_2_m = rho_2_rho * rho_m + &
                   rho_2_u   * u_m   + &
                   rho_2_v   * v_m   + &
                   rho_2_w   * w_m   + &
                   rho_2_p   * p_m 
        flux_3_m = rho_3_rho * rho_m + &
                   rho_3_u   * u_m   + &
                   rho_3_v   * v_m   + &
                   rho_3_w   * w_m   + &
                   rho_3_p   * p_m 


        flux_1_p = rho_1_rho * rho_p + &
                   rho_1_u   * u_p   + &
                   rho_1_v   * v_p   + &
                   rho_1_w   * w_p   + &
                   rho_1_p   * p_p
        flux_2_p = rho_2_rho * rho_p + &
                   rho_2_u   * u_p   + &
                   rho_2_v   * v_p   + &
                   rho_2_w   * w_p   + &
                   rho_2_p   * p_p 
        flux_3_p = rho_3_rho * rho_p + &
                   rho_3_u   * u_p   + &
                   rho_3_v   * v_p   + &
                   rho_3_w   * w_p   + &
                   rho_3_p   * p_p 


        call worker%integrate_boundary_average('Density(imag)','Advection', &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)
        


        !================================
        !       X-MOMENTUM FLUX
        !================================
        flux_1_m = u_1_rho * rho_m + &
                   u_1_u   * u_m   + &
                   u_1_v   * v_m   + &
                   u_1_w   * w_m   + &
                   u_1_p   * p_m
        flux_2_m = u_2_rho * rho_m + &
                   u_2_u   * u_m   + &
                   u_2_v   * v_m   + &
                   u_2_w   * w_m   + &
                   u_2_p   * p_m 
        flux_3_m = u_3_rho * rho_m + &
                   u_3_u   * u_m   + &
                   u_3_v   * v_m   + &
                   u_3_w   * w_m   + &
                   u_3_p   * p_m 


        flux_1_p = u_1_rho * rho_p + &
                   u_1_u   * u_p   + &
                   u_1_v   * v_p   + &
                   u_1_w   * w_p   + &
                   u_1_p   * p_p
        flux_2_p = u_2_rho * rho_p + &
                   u_2_u   * u_p   + &
                   u_2_v   * v_p   + &
                   u_2_w   * w_p   + &
                   u_2_p   * p_p 
        flux_3_p = u_3_rho * rho_p + &
                   u_3_u   * u_p   + &
                   u_3_v   * v_p   + &
                   u_3_w   * w_p   + &
                   u_3_p   * p_p 

        call worker%integrate_boundary_average('Velocity-1(imag)','Advection',  &
                                                flux_1_m,flux_2_m,flux_3_m,     &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        flux_1_m = v_1_rho * rho_m + &
                   v_1_u   * u_m   + &
                   v_1_v   * v_m   + &
                   v_1_w   * w_m   + &
                   v_1_p   * p_m
        flux_2_m = v_2_rho * rho_m + &
                   v_2_u   * u_m   + &
                   v_2_v   * v_m   + &
                   v_2_w   * w_m   + &
                   v_2_p   * p_m 
        flux_3_m = v_3_rho * rho_m + &
                   v_3_u   * u_m   + &
                   v_3_v   * v_m   + &
                   v_3_w   * w_m   + &
                   v_3_p   * p_m 


        flux_1_p = v_1_rho * rho_p + &
                   v_1_u   * u_p   + &
                   v_1_v   * v_p   + &
                   v_1_w   * w_p   + &
                   v_1_p   * p_p
        flux_2_p = v_2_rho * rho_p + &
                   v_2_u   * u_p   + &
                   v_2_v   * v_p   + &
                   v_2_w   * w_p   + &
                   v_2_p   * p_p 
        flux_3_p = v_3_rho * rho_p + &
                   v_3_u   * u_p   + &
                   v_3_v   * v_p   + &
                   v_3_w   * w_p   + &
                   v_3_p   * p_p 

        call worker%integrate_boundary_average('Velocity-2(imag)','Advection',  &
                                                flux_1_m,flux_2_m,flux_3_m,     &   
                                                flux_1_p,flux_2_p,flux_3_p)



        !================================
        !       Z-MOMENTUM FLUX
        !================================
        flux_1_m = w_1_rho * rho_m + &
                   w_1_u   * u_m   + &
                   w_1_v   * v_m   + &
                   w_1_w   * w_m   + &
                   w_1_p   * p_m
        flux_2_m = w_2_rho * rho_m + &
                   w_2_u   * u_m + &
                   w_2_v   * v_m + &
                   w_2_w   * w_m + &
                   w_2_p   * p_m 
        flux_3_m = w_3_rho * rho_m + &
                   w_3_u   * u_m + &
                   w_3_v   * v_m + &
                   w_3_w   * w_m + &
                   w_3_p   * p_m 


        flux_1_p = w_1_rho * rho_p + &
                   w_1_u   * u_p   + &
                   w_1_v   * v_p   + &
                   w_1_w   * w_p   + &
                   w_1_p   * p_p
        flux_2_p = w_2_rho * rho_p + &
                   w_2_u   * u_p   + &
                   w_2_v   * v_p   + &
                   w_2_w   * w_p   + &
                   w_2_p   * p_p 
        flux_3_p = w_3_rho * rho_p + &
                   w_3_u   * u_p   + &
                   w_3_v   * v_p   + &
                   w_3_w   * w_p   + &
                   w_3_p   * p_p 


        call worker%integrate_boundary_average('Velocity-3(imag)','Advection',  &
                                                flux_1_m,flux_2_m,flux_3_m,     &   
                                                flux_1_p,flux_2_p,flux_3_p)




        !================================
        !          ENERGY FLUX
        !================================
        flux_1_m = p_1_rho * rho_m + &
                   p_1_u   * u_m   + &
                   p_1_v   * v_m   + &
                   p_1_w   * w_m   + &
                   p_1_p   * p_m
        flux_2_m = p_2_rho * rho_m + &
                   p_2_u   * u_m   + &
                   p_2_v   * v_m   + &
                   p_2_w   * w_m   + &
                   p_2_p   * p_m 
        flux_3_m = p_3_rho * rho_m + &
                   p_3_u   * u_m   + &
                   p_3_v   * v_m   + &
                   p_3_w   * w_m   + &
                   p_3_p   * p_m 


        flux_1_p = p_1_rho * rho_p + &
                   p_1_u   * u_p   + &
                   p_1_v   * v_p   + &
                   p_1_w   * w_p   + &
                   p_1_p   * p_p
        flux_2_p = p_2_rho * rho_p + &
                   p_2_u   * u_p   + &
                   p_2_v   * v_p   + &
                   p_2_w   * w_p   + &
                   p_2_p   * p_p 
        flux_3_p = p_3_rho * rho_p + &
                   p_3_u   * u_p   + &
                   p_3_v   * v_p   + &
                   p_3_w   * w_p   + &
                   p_3_p   * p_p 


        call worker%integrate_boundary_average('Pressure(imag)','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m,     &   
                                                flux_1_p,flux_2_p,flux_3_p)

    end subroutine compute
    !************************************************************************************************************












end module PRIMLINEULER_boundary_average_advective_flux_imag
