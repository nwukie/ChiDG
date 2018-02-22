module graddemo_P_bc_operator
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
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: graddemo_P_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_P_bc_operator_t
    !***************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_P_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Graddemo P BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !**************************************************************************************





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
    !---------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(graddemo_P_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_pbc, grad2_pbc, grad3_pbc,        &
            grad1_p,   grad2_p,   grad3_p,          &
            lift1,     lift2,     lift3,            &
            term1,     term2,     term3,            &
            flux_1,    flux_2,    flux_3

        type(AD_D), allocatable, dimension(:)   ::          &
            density, mom1, mom2, mom3, u, v, w, invdensity, &
            grad1_density,  grad2_density,  grad3_density,  &
            grad1_mom1,     grad2_mom1,     grad3_mom1,     &
            grad1_mom2,     grad2_mom2,     grad3_mom2,     &
            grad1_mom3,     grad2_mom3,     grad3_mom3,     &
            du_ddensity, dv_ddensity, dw_ddensity,          &
            du_dmom1, dv_dmom2, dw_dmom3,                   &
            grad1_u, grad2_v, grad1_w, grad2_w, grad3_w, dpdz,                &
            p_element, p_face, dp, P, T, grad3_T, c, lambda1, lambda5, L5, L1, T1, K, pressure

        real(rk),   allocatable, dimension(:)   :: z_element, z_face, dz

!        if ( (worker%iface == 1) .or. &
!             (worker%iface == 2) .or. &
!             (worker%iface == 3) .or. &
!             (worker%iface == 4) ) then



!        !
!        ! Interpolate boundary condition state to face quadrature nodes
!        !
!        grad1_pbc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary', only_lift=.true.)
!        grad2_pbc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary', only_lift=.true.)
!        grad3_pbc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary', only_lift=.true.)
!
!        grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
!        grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
!        grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')
!
!
!        lift1 = grad1_pbc
!        lift2 = grad2_pbc
!        lift3 = grad3_pbc
!
!        !=================================================
!        ! Mass flux
!        !=================================================
!        !flux_1 = grad1_pbc - 132000._rk*sin(132._rk * worker%x('boundary'))
!        !flux_1 = grad1_pbc - grad1_p
!        !flux_2 = grad2_pbc - grad2_p
!        !flux_3 = grad3_pbc 
!        flux_1 = (grad1_p+lift1) - grad1_p
!        flux_2 = (grad2_p+lift2) - grad2_p
!        flux_3 = (grad3_p+lift3) - grad3_p
!        !flux_2 = (grad2_p+lift2)




        !if ( (worker%element_info%idomain_g == 5) .and. (worker%element_info%ielement_g == 3) .and. (worker%iface == 2)) then
        !if ( (worker%element_info%ielement_g == 1) .and. (worker%iface == 5) ) then
        if ( (worker%element_info%ielement_g == 40) .and. (worker%iface == 6) ) then
        !if ( (worker%element_info%ielement_g == 160) .and. (worker%iface == 6) ) then
        !if ( (worker%element_info%ielement_g == 640) .and. (worker%iface == 6) ) then
        !if ( (worker%element_info%ielement_g == 800) .and. (worker%iface == 6) ) then
        !if ( (worker%element_info%ielement_g == 1000) .and. (worker%iface == 6) ) then

            !
            ! Interpolate boundary condition state to face quadrature nodes
            !
            pressure  = worker%get_field('Pressure', 'value', 'boundary')
            grad1_pbc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary')
            grad2_pbc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary')
            grad3_pbc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary')
            lift1 = worker%get_field('Pressure_TEMP', 'grad1', 'boundary', only_lift=.true.)
            lift2 = worker%get_field('Pressure_TEMP', 'grad2', 'boundary', only_lift=.true.)
            lift3 = worker%get_field('Pressure_TEMP', 'grad3', 'boundary', only_lift=.true.)

            grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
            grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
            grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')


            !
            ! Compute axial pressure gradient
            !
            !----------------------------------------------------------------
            density = worker%get_field('Density',    'value', 'boundary')
            mom1    = worker%get_field('Momentum-1', 'value', 'boundary')
            mom2    = worker%get_field('Momentum-2', 'value', 'boundary')
            mom3    = worker%get_field('Momentum-3', 'value', 'boundary')


            grad1_density = worker%get_field('Density', 'grad1', 'boundary')
            grad2_density = worker%get_field('Density', 'grad2', 'boundary')
            grad3_density = worker%get_field('Density', 'grad3', 'boundary')

            grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'boundary')
            grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'boundary')
            grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'boundary')

            grad1_mom2 = worker%get_field('Momentum-2', 'grad1', 'boundary')
            grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'boundary')
            grad3_mom2 = worker%get_field('Momentum-2', 'grad3', 'boundary')

            grad1_mom3 = worker%get_field('Momentum-3', 'grad1', 'boundary')
            grad2_mom3 = worker%get_field('Momentum-3', 'grad2', 'boundary')
            grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'boundary')

            invdensity = ONE/density
            u = mom1/density
            v = mom2/density
            w = mom3/density


            !
            ! compute velocity jacobians
            !
            du_ddensity  = -invdensity*invdensity*mom1
            dv_ddensity  = -invdensity*invdensity*mom2
            dw_ddensity  = -invdensity*invdensity*mom3

            du_dmom1 = invdensity
            dv_dmom2 = invdensity
            dw_dmom3 = invdensity



            !
            ! compute velocity gradients via chain rule:
            !
            !   u = f(rho,rhou)
            !
            !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
            !
            grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
            grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2

            grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
            grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
            grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3


            term1 = grad1_mom3*u  + grad1_u*mom3
            term2 = grad2_mom3*v  + grad2_v*mom3
            term3 = grad3_mom3*w  + grad3_w*mom3
            dpdz = -(term1 + term2 + term3)
            !!----------------------------------------------------------------
            !!
            !!   Compute dpdz from interior
            !!
            !!----------------------------------------------------------------
            !p_element = worker%get_field('Pressure', 'value', 'element')
            !p_face    = worker%get_field('Pressure', 'value', 'boundary')
            !z_element = worker%coordinate('3', 'element')
            !z_face    = worker%coordinate('3', 'boundary')
            !! P1
            !dp = p_face - p_element(19:27)  
            !dz = z_face - z_element(19:27)  
            !! P2           
            !!dp = p_face - p_element(101:125)
            !!dz = z_face - z_element(101:125)
            !! P3           
            !!dp = p_face - p_element(181:216)
            !!dz = z_face - z_element(181:216)
            !! P4           
            !!dp = p_face - p_element(449:512)
            !!dz = z_face - z_element(449:512)
            !! P5
            !!dp = p_face - p_element(649:729) 
            !!dz = z_face - z_element(649:729) 
            !dpdz = dp/dz

            !----------------------------------------------------------------
            !
            !   Entropy gradient
            !
            !----------------------------------------------------------------
            !P       = worker%get_field('Pressure_TEMP',            'value', 'boundary')
            !T       = worker%get_field('Temperature',              'value', 'boundary')
            !grad3_T = worker%get_field('Temperature Gradient - 3', 'value', 'boundary')
            !dpdz = ((gam/(gam-ONE))*T**(ONE/(gam-ONE)))*grad3_T * (P/(T**(gam/(gam-ONE))))

            !P             = worker%get_field('Pressure', 'value', 'boundary')
            !density       = worker%get_field('Density',  'value', 'boundary')
            !grad3_density = worker%get_field('Density',  'grad3', 'boundary')
            !dpdz = gam*(P/density)*grad3_density

            !----------------------------------------------------------------
            !
            !   dpdz = constant
            !
            !----------------------------------------------------------------
            dpdz = ZERO


            
            !----------------------------------------------------------------
            !
            !   LODI
            !
            !----------------------------------------------------------------
            T  = worker%get_field('Temperature', 'value', 'boundary')
            c  = sqrt(gam*Rgas*T)
            lambda1 = w - c
            lambda5 = w + c
            K = c
            K = 10._rk
            T1 = HALF*( (u*grad1_p + v*grad2_p)  +  gam*pressure*(grad1_u + grad2_v)  -  density*c*(u*grad1_w + v*grad2_w) )
            L1 = K*(pressure - 100000._rk)  +  (0.1_rk - ONE)*T1
            L5 = lambda5*(grad3_p + density*c*grad3_w)
            dpdz = HALF*(L5/(w + c) + L1/(w - c))





            grad1_pbc = grad1_p
            grad2_pbc = grad2_p
            grad3_pbc = dpdz

            !flux_1 = grad1_pbc - grad1_p
            !flux_2 = grad2_pbc - grad2_p
            !flux_3 = grad3_pbc - grad3_p
            flux_1 = (grad1_pbc+lift1) - grad1_p
            flux_2 = (grad2_pbc+lift2) - grad2_p
            flux_3 = (grad3_pbc+lift3) - grad3_p


        !else if ( (worker%element_info%idomain_g == 5) .and. (worker%element_info%ielement_g /= 3) .and. (worker%iface == 2)) then
        !else if ( (worker%element_info%ielement_g /= 1) .and. (worker%iface == 5) ) then
        else if ( (worker%element_info%ielement_g /= 40) .and. (worker%iface == 6) ) then
        !else if ( (worker%element_info%ielement_g /= 160) .and. (worker%iface == 6) ) then
        !else if ( (worker%element_info%ielement_g /= 640) .and. (worker%iface == 6) ) then
        !else if ( (worker%element_info%ielement_g /= 800) .and. (worker%iface == 6) ) then
        !else if ( (worker%element_info%ielement_g /= 1000) .and. (worker%iface == 6) ) then


            !
            ! Interpolate boundary condition state to face quadrature nodes
            !
            pressure  = worker%get_field('Pressure', 'value', 'boundary')
            grad1_pbc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary')
            grad2_pbc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary')
            grad3_pbc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary')
            lift1 = worker%get_field('Pressure_TEMP', 'grad1', 'boundary', only_lift=.true.)
            lift2 = worker%get_field('Pressure_TEMP', 'grad2', 'boundary', only_lift=.true.)
            lift3 = worker%get_field('Pressure_TEMP', 'grad3', 'boundary', only_lift=.true.)

            grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
            grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
            grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')


            !
            ! Compute axial pressure gradient
            !
            !----------------------------------------------------------------
            density = worker%get_field('Density',    'value', 'boundary')
            mom1    = worker%get_field('Momentum-1', 'value', 'boundary')
            mom2    = worker%get_field('Momentum-2', 'value', 'boundary')
            mom3    = worker%get_field('Momentum-3', 'value', 'boundary')


            grad1_density = worker%get_field('Density', 'grad1', 'boundary')
            grad2_density = worker%get_field('Density', 'grad2', 'boundary')
            grad3_density = worker%get_field('Density', 'grad3', 'boundary')

            grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'boundary')
            grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'boundary')
            grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'boundary')

            grad1_mom2 = worker%get_field('Momentum-2', 'grad1', 'boundary')
            grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'boundary')
            grad3_mom2 = worker%get_field('Momentum-2', 'grad3', 'boundary')

            grad1_mom3 = worker%get_field('Momentum-3', 'grad1', 'boundary')
            grad2_mom3 = worker%get_field('Momentum-3', 'grad2', 'boundary')
            grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'boundary')

            invdensity = ONE/density
            u = mom1/density
            v = mom2/density
            w = mom3/density


            !
            ! compute velocity jacobians
            !
            du_ddensity  = -invdensity*invdensity*mom1
            dv_ddensity  = -invdensity*invdensity*mom2
            dw_ddensity  = -invdensity*invdensity*mom3

            du_dmom1 = invdensity
            dv_dmom2 = invdensity
            dw_dmom3 = invdensity



            !
            ! compute velocity gradients via chain rule:
            !
            !   u = f(rho,rhou)
            !
            !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
            !
            grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
            grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2

            grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
            grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
            grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3


            term1 = grad1_mom3*u  + grad1_u*mom3
            term2 = grad2_mom3*v  + grad2_v*mom3
            term3 = grad3_mom3*w  + grad3_w*mom3
            dpdz = -(term1 + term2 + term3)
            !-----------------------------------------------------------------

            !!----------------------------------------------------------------
            !!
            !!   Compute dpdz from interior
            !!
            !!----------------------------------------------------------------
            !p_element = worker%get_field('Pressure', 'value', 'element')
            !p_face    = worker%get_field('Pressure', 'value', 'boundary')
            !z_element = worker%coordinate('3', 'element')
            !z_face    = worker%coordinate('3', 'boundary')
            !! P1
            !dp = p_face - p_element(19:27)  
            !dz = z_face - z_element(19:27)  
            !! P2           
            !!dp = p_face - p_element(101:125)
            !!dz = z_face - z_element(101:125)
            !! P3           
            !!dp = p_face - p_element(181:216)
            !!dz = z_face - z_element(181:216)
            !! P4           
            !!dp = p_face - p_element(449:512)
            !!dz = z_face - z_element(449:512)
            !! P5
            !!dp = p_face - p_element(649:729)
            !!dz = z_face - z_element(649:729)
            !dpdz = dp/dz

            !!----------------------------------------------------------------
            !!
            !!   Compute dpdz from entropy
            !!
            !!----------------------------------------------------------------
            !p_element = worker%get_field('Pressure', 'value', 'element' )
            !p_face    = worker%get_field('Pressure', 'value', 'boundary')
            !z_element = worker%coordinate('3', 'element' )
            !z_face    = worker%coordinate('3', 'boundary')
            !t_element = worker%get_field('Temperature', 'value', 'element' )
            !t_face    = worker%get_field('Temperature', 'value', 'boundary')

            !! P1
            !dp = p_face - p_element(19:27)  
            !dz = z_face - z_element(19:27)  
            !! P2           
            !!dp = p_face - p_element(101:125)
            !!dz = z_face - z_element(101:125)
            !! P3           
            !!dp = p_face - p_element(181:216)
            !!dz = z_face - z_element(181:216)
            !! P4           
            !!dp = p_face - p_element(449:512)
            !!dz = z_face - z_element(449:512)
            !! P5
            !!dp = p_face - p_element(649:729)
            !!dz = z_face - z_element(649:729)
            !dpdz = dp/dz


            !----------------------------------------------------------------
            !
            !   Entropy gradient
            !
            !----------------------------------------------------------------
            !P       = worker%get_field('Pressure_TEMP',            'value', 'boundary')
            !T       = worker%get_field('Temperature',              'value', 'boundary')
            !grad3_T = worker%get_field('Temperature Gradient - 3', 'value', 'boundary')
            !dpdz = ((gam/(gam-ONE))*T**(ONE/(gam-ONE)))*grad3_T * (P/(T**(gam/(gam-ONE))))

            !P             = worker%get_field('Pressure', 'value', 'boundary')
            !density       = worker%get_field('Density',  'value', 'boundary')
            !grad3_density = worker%get_field('Density',  'grad3', 'boundary')
            !dpdz = gam*(P/density)*grad3_density


            !----------------------------------------------------------------
            !
            !   dpdz = constant
            !
            !----------------------------------------------------------------
            dpdz = ZERO


            !----------------------------------------------------------------
            !
            !   LODI
            !
            !----------------------------------------------------------------
            T  = worker%get_field('Temperature', 'value', 'boundary')
            c  = sqrt(gam*Rgas*T)
            lambda5 = w + c
            lambda1 = w - c
            K = c
            K = 10._rk
            T1 = HALF*( (u*grad1_p + v*grad2_p)  +  gam*pressure*(grad1_u + grad2_v)  -  density*c*(u*grad1_w + v*grad2_w) )
            L1 = K*(pressure - 100000._rk)  +  (0.1_rk - ONE)*T1
            L5 = lambda5*(grad3_p + density*c*grad3_w)
            dpdz = HALF*(L5/(w + c) + L1/(w - c))







            !=================================================
            ! Mass flux
            !=================================================
            !flux_1 = grad1_pbc - 132000._rk*sin(132._rk * worker%x('boundary'))
            !flux_1 = grad1_pbc - grad1_p
            !flux_2 = grad2_pbc - grad2_p
            !flux_3 = grad3_pbc 
            !flux_1 = (grad1_p+lift1) - grad1_p
            !!flux_2 = (grad2_p+lift2) - grad2_p
            !flux_2 = grad2_p
            !flux_3 = (grad3_p+lift3) - grad3_p

            !flux_1 = (grad1_p+lift1) - grad1_p
            !flux_2 = (grad2_p+lift2) - grad2_p
            !flux_3 = (grad3_p+lift3) - grad3_p

            grad1_pbc = grad1_p
            grad2_pbc = grad2_p
            grad3_pbc = dpdz


            flux_1 = grad1_pbc - grad1_p
            flux_2 = grad2_pbc - grad2_p
            flux_3 = grad3_pbc - grad3_p






        else


            !
            ! Interpolate boundary condition state to face quadrature nodes
            !
            grad1_pbc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary')
            grad2_pbc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary')
            grad3_pbc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary')

            grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
            grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
            grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')


            !=================================================
            ! Mass flux
            !=================================================
            !flux_1 = grad1_pbc - grad1_p
            !!flux_2 = grad2_pbc - grad2_p
            !flux_2 = grad2_pbc
            !flux_3 = grad3_pbc - grad3_p
            flux_1 = grad1_pbc - grad1_p
            flux_2 = grad2_pbc - grad2_p
            flux_3 = grad3_pbc - grad3_p

        end if



        call worker%integrate_boundary_condition('Pressure_TEMP','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !***************************************************************************************










end module graddemo_P_bc_operator
