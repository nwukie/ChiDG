module bc_primlineuler_extrapolate
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO, ME

    use type_bc,            only: bc_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    
    use PRIMLINEULER_properties,   only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: primlineuler_extrapolate_t

    contains

        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type primlineuler_extrapolate_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(primlineuler_extrapolate_t),    intent(inout)   :: self

        !
        ! Set name
        ! 
        call self%set_name('primlineuler_extrapolate')


    end subroutine add_options
    !*******************************************************************************************
















    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(primlineuler_extrapolate_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho_r, iu_r, iv_r, iw_r, ip_r
        integer(ik)     :: irho_i, iu_i, iv_i, iw_i, ip_i


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r,  u_r, v_r, w_r, p_r,             &
            rho_i,  u_i, v_i, w_i, p_i,             &
            flux_x, flux_y, flux_z, integrand

        real(rK),   allocatable, dimension(:)   ::  &
            normx, normy, normz




        !
        ! Get equation indices
        !
        irho_r  = prop%get_eqn_index("rho_r")
        iu_r    = prop%get_eqn_index("u_r")
        iv_r    = prop%get_eqn_index("v_r")
        iw_r    = prop%get_eqn_index("w_r")
        ip_r    = prop%get_eqn_index("p_r")


        irho_i  = prop%get_eqn_index("rho_i")
        iu_i    = prop%get_eqn_index("u_i")
        iv_i    = prop%get_eqn_index("v_i")
        iw_i    = prop%get_eqn_index("w_i")
        ip_i    = prop%get_eqn_index("p_i")



        !
        ! Interpolate interior solution to quadrature nodes
        !
        rho_r = worker%interpolate(irho_r, 'value', ME)
        u_r   = worker%interpolate(iu_r,   'value', ME)
        v_r   = worker%interpolate(iv_r,   'value', ME)
        w_r   = worker%interpolate(iw_r,   'value', ME)
        p_r   = worker%interpolate(ip_r,   'value', ME)


        rho_i = worker%interpolate(irho_i, 'value', ME)
        u_i   = worker%interpolate(iu_i,   'value', ME)
        v_i   = worker%interpolate(iv_i,   'value', ME)
        w_i   = worker%interpolate(iw_i,   'value', ME)
        p_i   = worker%interpolate(ip_i,   'value', ME)






        !=================================================
        ! Mass flux
        !=================================================
        flux_x = rho_x_rho  * rho_r  + &
                 rho_x_u    * u_r    + &
                 rho_x_v    * v_r    + &
                 rho_x_w    * w_r    + &
                 rho_x_p    * p_r
        flux_y = rho_y_rho  * rho_r  + &
                 rho_y_u    * u_r    + &
                 rho_y_v    * v_r    + &
                 rho_y_w    * w_r    + &
                 rho_y_p    * p_r
        flux_z = rho_z_rho  * rho_r  + &
                 rho_z_u    * u_r    + &
                 rho_z_v    * v_r    + &
                 rho_z_w    * w_r    + &
                 rho_z_p    * p_r

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irho_r, integrand)





        flux_x = rho_x_rho  * rho_i  + &
                 rho_x_u    * u_i    + &
                 rho_x_v    * v_i    + &
                 rho_x_w    * w_i    + &
                 rho_x_p    * p_i
        flux_y = rho_y_rho  * rho_i  + &
                 rho_y_u    * u_i    + &
                 rho_y_v    * v_i    + &
                 rho_y_w    * w_i    + &
                 rho_y_p    * p_i
        flux_z = rho_z_rho  * rho_i  + &
                 rho_z_u    * u_i    + &
                 rho_z_v    * v_i    + &
                 rho_z_w    * w_i    + &
                 rho_z_p    * p_i

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irho_i, integrand)









        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = u_x_rho  * rho_r  + &
                 u_x_u    * u_r    + &
                 u_x_v    * v_r    + &
                 u_x_w    * w_r    + &
                 u_x_p    * p_r
        flux_y = u_y_rho  * rho_r  + &
                 u_y_u    * u_r    + &
                 u_y_v    * v_r    + &
                 u_y_w    * w_r    + &
                 u_y_p    * p_r
        flux_z = u_z_rho  * rho_r  + &
                 u_z_u    * u_r    + &
                 u_z_v    * v_r    + &
                 u_z_w    * w_r    + &
                 u_z_p    * p_r

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iu_r, integrand)


        flux_x = u_x_rho  * rho_i  + &
                 u_x_u    * u_i    + &
                 u_x_v    * v_i    + &
                 u_x_w    * w_i    + &
                 u_x_p    * p_i
        flux_y = u_y_rho  * rho_i  + &
                 u_y_u    * u_i    + &
                 u_y_v    * v_i    + &
                 u_y_w    * w_i    + &
                 u_y_p    * p_i
        flux_z = u_z_rho  * rho_i  + &
                 u_z_u    * u_i    + &
                 u_z_v    * v_i    + &
                 u_z_w    * w_i    + &
                 u_z_p    * p_i

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iu_i, integrand)







        !=================================================
        ! y-momentum flux
        !=================================================

        flux_x = v_x_rho  * rho_r  + &
                 v_x_u    * u_r    + &
                 v_x_v    * v_r    + &
                 v_x_w    * w_r    + &
                 v_x_p    * p_r
        flux_y = v_y_rho  * rho_r  + &
                 v_y_u    * u_r    + &
                 v_y_v    * v_r    + &
                 v_y_w    * w_r    + &
                 v_y_p    * p_r
        flux_z = v_z_rho  * rho_r  + &
                 v_z_u    * u_r    + &
                 v_z_v    * v_r    + &
                 v_z_w    * w_r    + &
                 v_z_p    * p_r

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iv_r, integrand)


        flux_x = v_x_rho  * rho_i  + &
                 v_x_u    * u_i    + &
                 v_x_v    * v_i    + &
                 v_x_w    * w_i    + &
                 v_x_p    * p_i
        flux_y = v_y_rho  * rho_i  + &
                 v_y_u    * u_i    + &
                 v_y_v    * v_i    + &
                 v_y_w    * w_i    + &
                 v_y_p    * p_i
        flux_z = v_z_rho  * rho_i  + &
                 v_z_u    * u_i    + &
                 v_z_v    * v_i    + &
                 v_z_w    * w_i    + &
                 v_z_p    * p_i

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iv_i, integrand)











        !=================================================
        ! z-momentum flux
        !=================================================

        flux_x = w_x_rho  * rho_r  + &
                 w_x_u    * u_r    + &
                 w_x_v    * v_r    + &
                 w_x_w    * w_r    + &
                 w_x_p    * p_r
        flux_y = w_y_rho  * rho_r  + &
                 w_y_u    * u_r    + &
                 w_y_v    * v_r    + &
                 w_y_w    * w_r    + &
                 w_y_p    * p_r
        flux_z = w_z_rho  * rho_r  + &
                 w_z_u    * u_r    + &
                 w_z_v    * v_r    + &
                 w_z_w    * w_r    + &
                 w_z_p    * p_r

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iw_r, integrand)


        flux_x = w_x_rho  * rho_i  + &
                 w_x_u    * u_i    + &
                 w_x_v    * v_i    + &
                 w_x_w    * w_i    + &
                 w_x_p    * p_i
        flux_y = w_y_rho  * rho_i  + &
                 w_y_u    * u_i    + &
                 w_y_v    * v_i    + &
                 w_y_w    * w_i    + &
                 w_y_p    * p_i
        flux_z = w_z_rho  * rho_i  + &
                 w_z_u    * u_i    + &
                 w_z_v    * v_i    + &
                 w_z_w    * w_i    + &
                 w_z_p    * p_i

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iw_i, integrand)







        !=================================================
        ! Energy flux
        !=================================================


        flux_x = p_x_rho  * rho_r  + &
                 p_x_u    * u_r    + &
                 p_x_v    * v_r    + &
                 p_x_w    * w_r    + &
                 p_x_p    * p_r
        flux_y = p_y_rho  * rho_r  + &
                 p_y_u    * u_r    + &
                 p_y_v    * v_r    + &
                 p_y_w    * w_r    + &
                 p_y_p    * p_r
        flux_z = p_z_rho  * rho_r  + &
                 p_z_u    * u_r    + &
                 p_z_v    * v_r    + &
                 p_z_w    * w_r    + &
                 p_z_p    * p_r

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(ip_r, integrand)


        flux_x = p_x_rho  * rho_i  + &
                 p_x_u    * u_i    + &
                 p_x_v    * v_i    + &
                 p_x_w    * w_i    + &
                 p_x_p    * p_i
        flux_y = p_y_rho  * rho_i  + &
                 p_y_u    * u_i    + &
                 p_y_v    * v_i    + &
                 p_y_w    * w_i    + &
                 p_y_p    * p_i
        flux_z = p_z_rho  * rho_i  + &
                 p_z_u    * u_i    + &
                 p_z_v    * v_i    + &
                 p_z_w    * w_i    + &
                 p_z_p    * p_i

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(ip_i, integrand)





    end subroutine compute
    !**********************************************************************************************************






end module bc_primlineuler_extrapolate
