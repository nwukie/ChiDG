module RANS_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use mod_fluid,              only: omega
    use mod_spalart_allmaras

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !>  Spalart-Allmaras Source Term.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: RANS_source_t

        real(rk)    :: gam = 1.4_rk
        real(rk)    :: R   = 287.15_rk
        real(rk)    :: Cp  = 1003.0_rk
        real(rk)    :: Pr  = 0.72_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type RANS_source_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(RANS_source_t),   intent(inout)      :: self

        ! Set operator name.
        call self%set_name('RANS Source')

        ! Set operator type.
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations being integrated.
        call self%add_primary_field('Density * NuTilde')

        ! Set auxiliary variables being used.
        call self%add_auxiliary_field('Wall Distance : p-Poisson')
        call self%add_model('Wall Distance : p-Poisson Normalization')

        ! Add Turbulent Egrad2 Viscosity model
        !call self%add_model('Spalart Allmaras Turbulent Model Fields')

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(RANS_source_t),           intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                            &
            rho, mom1, mom2, mom3, energy, rho_nutilde, invrho,             &
            mu, nu, lamda,                                                  &
            grad1_rho,          grad2_rho,          grad3_rho,              &
            grad1_mom1,         grad2_mom1,         grad3_mom1,             &
            grad1_mom2,         grad2_mom2,         grad3_mom2,             &
            grad1_mom3,         grad2_mom3,         grad3_mom3,             &
            grad1_u,            grad2_u,            grad3_u,                &
            grad1_v,            grad2_v,            grad3_v,                &
            grad1_w,            grad2_w,            grad3_w,                &
            grad1_rho_nutilde,  grad2_rho_nutilde,  grad3_rho_nutilde,      &
            du_ddensity, dv_ddensity, dw_ddensity,                          &
            du_dmom1,    dv_dmom2,    dw_dmom3,                             &
            vorticity2, vorticity, vorticity_bar, vorticity_mod,            &
            f_v1, f_v2, f_t2, f_w, f_n1, production, destruction,           &
            v,                                                              &
            dnutilde_drho, dnutilde_drho_nutilde, chi, r, rbar, g, mu_t,    &
            nutilde, grad1_nutilde, grad2_nutilde, grad3_nutilde,           &
            source, dwall, vorticity_1, vorticity_2, vorticity_3

        real(rk)    :: const, epsilon_vorticity, gam

        ! Sutherlands Law constants
        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]



        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_primary_field_element('Density',           'value')
        mom1        = worker%get_primary_field_element('Momentum-1',        'value')
        mom2        = worker%get_primary_field_element('Momentum-2',        'value')
        mom3        = worker%get_primary_field_element('Momentum-3',        'value')
        energy      = worker%get_primary_field_element('Energy',            'value')
        rho_nutilde = worker%get_primary_field_element('Density * NuTilde', 'value')



        !
        ! Interpolate solution gradients to quadrature nodes
        !
        grad1_rho         = worker%get_primary_field_element('Density',          'grad1+lift')
        grad2_rho         = worker%get_primary_field_element('Density',          'grad2+lift')
        grad3_rho         = worker%get_primary_field_element('Density',          'grad3+lift')

        grad1_mom1        = worker%get_primary_field_element('Momentum-1',       'grad1+lift')
        grad2_mom1        = worker%get_primary_field_element('Momentum-1',       'grad2+lift')
        grad3_mom1        = worker%get_primary_field_element('Momentum-1',       'grad3+lift')

        grad1_mom2        = worker%get_primary_field_element('Momentum-2',       'grad1+lift')
        grad2_mom2        = worker%get_primary_field_element('Momentum-2',       'grad2+lift')
        grad3_mom2        = worker%get_primary_field_element('Momentum-2',       'grad3+lift')

        grad1_mom3        = worker%get_primary_field_element('Momentum-3',       'grad1+lift')
        grad2_mom3        = worker%get_primary_field_element('Momentum-3',       'grad2+lift')
        grad3_mom3        = worker%get_primary_field_element('Momentum-3',       'grad3+lift')

        grad1_rho_nutilde = worker%get_primary_field_element('Density * NuTilde','grad1+lift')
        grad2_rho_nutilde = worker%get_primary_field_element('Density * NuTilde','grad2+lift')
        grad3_rho_nutilde = worker%get_primary_field_element('Density * NuTilde','grad3+lift')



        !
        ! Interpolate auxiliary field, Wall Distance
        !
        dwall = worker%get_model_field_element('Wall Distance', 'value')



        !
        ! Divide by density
        !
        invrho  = ONE/rho
        nutilde = rho_nutilde*invrho



        !
        ! Compute model values
        !
        gam = 1.4_rk

        !pressure    = (self%gam-ONE)*(energy - HALF*( (mom1*mom1) + (mom2*mom2) + (mom3*mom3) )*invrho )
        !temperature = pressure/(rho*self%R)

        
        mu = worker%get_model_field_element('Laminar Viscosity', 'value')
        !mu  = mu0*((temperature/T0)**(THREE/TWO))*(T0+S)/(temperature+S)
        !mu  = worker%get_model_field_element('Laminar Viscosity',   'value')
        nu  = mu*invrho



        !
        ! Compute turbulence viscosity
        !
        chi  = nutilde/nu
        f_v1 = (chi*chi*chi)/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)


        mu_t = rho_nutilde
        mu_t = ZERO
        where(nutilde >= ZERO)
            mu_t = rho_nutilde * f_v1
        end where
        


        !
        ! Compute f_n1
        !
        f_n1 = rho_nutilde
        f_n1 = ONE
        where(nutilde < ZERO)
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where



        !
        ! Compute vorticity and modified vorticity
        !

        !############### Computing Vorticity ######################

        !
        ! compute velocity jacobians
        !
        du_ddensity = -invrho*invrho*mom1
        dv_ddensity = -invrho*invrho*mom2
        dw_ddensity = -invrho*invrho*mom3

        du_dmom1 = invrho
        dv_dmom2 = invrho
        dw_dmom3 = invrho



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = du_ddensity*grad1_rho  +  du_dmom1*grad1_mom1
        grad2_u = du_ddensity*grad2_rho  +  du_dmom1*grad2_mom1
        grad3_u = du_ddensity*grad3_rho  +  du_dmom1*grad3_mom1

        grad1_v = dv_ddensity*grad1_rho  +  dv_dmom2*grad1_mom2
        grad2_v = dv_ddensity*grad2_rho  +  dv_dmom2*grad2_mom2
        grad3_v = dv_ddensity*grad3_rho  +  dv_dmom2*grad3_mom2

        grad1_w = dw_ddensity*grad1_rho  +  dw_dmom3*grad1_mom3
        grad2_w = dw_ddensity*grad2_rho  +  dw_dmom3*grad2_mom3
        grad3_w = dw_ddensity*grad3_rho  +  dw_dmom3*grad3_mom3




        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute vorticity:
            !
            !   vorticity = Curl(U)
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   curl(U) = (dwdx - dvdz)i  +  (dudz - dwdx)j  +  (dvdx - dudy)k
            !
            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u) 



        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   curl(U) = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  (1/r)(d(rv)dr - dudtheta)k
            !           = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  ( dvdr - (1/r)dudtheta  + (v/r) )k
            !           = (grad2_w - grad3_v)i  +  (grad3_u - grad1_w)j  +  ( grad1_v - grad2_u  + (v/r) )k
            !
            !   Note:
            !       grad1_ = d/dr
            !       grad2_ = (1/r)d/dtheta
            !       grad3_ = d/dz
            !
            v = mom2*invrho

            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u + (v/r)) 


            !
            ! Account for rotation, convert to relative vorticity
            !
            vorticity_3 = vorticity_3 - TWO*omega

        end if












        !vorticity_1 = worker%get_model_field_element('Vorticity-1', 'value')
        !vorticity_2 = worker%get_model_field_element('Vorticity-2', 'value')
        !vorticity_3 = worker%get_model_field_element('Vorticity-3', 'value')

        !vorticity2 =  (dw_dy - dv_dz)**TWO  +  (du_dz - dw_dx)**TWO  +  (dv_dx - du_dy)**TWO 
        vorticity2 =  vorticity_1**TWO  +  vorticity_2**TWO  +  vorticity_3**TWO 
        
        epsilon_vorticity = 1.e-6_rk
        vorticity = vorticity2
        where(vorticity2 < epsilon_vorticity)
            vorticity = HALF*(epsilon_vorticity + vorticity2/epsilon_vorticity)
        else where
            vorticity = sqrt(vorticity2)
        end where


        f_v2 = ONE - (chi/(ONE+chi*f_v1))
        vorticity_bar = (nutilde/(SA_kappa*SA_kappa*dwall*dwall))*f_v2


        vorticity_mod = vorticity
        where (vorticity_bar >= -SA_c_v2*vorticity)
            vorticity_mod = vorticity + vorticity_bar
        else where
            vorticity_mod = vorticity + vorticity*(SA_c_v2*SA_c_v2*vorticity + SA_c_v3*vorticity_bar)/( (SA_c_v3 - TWO*SA_c_v2)*vorticity - vorticity_bar ) 
        end where




        !
        ! Compute f_t2, f_w, g, r
        !
        rbar = nutilde/(vorticity_mod * SA_kappa * SA_kappa * dwall * dwall )

        !r = min(SA_rlim,rbar)
        r = SA_rlim - (SA_rlim - rbar)*(atan(SA_b * (SA_rlim - rbar))/PI  + HALF)  + atan(SA_b)/PI  -  HALF

        g = r + SA_c_w2*(r**SIX  -  r)

        f_w = g*(((ONE + SA_c_w3**SIX)/(g**SIX + SA_c_w3**SIX))**(ONE/SIX))

        f_t2 = SA_c_t3*exp(-SA_c_t4*chi*chi)



        !
        ! Compute Production, Destruction
        !
        production = vorticity_mod
        where ( nutilde >= ZERO )
            production = SA_c_b1*(ONE - f_t2)*vorticity_mod*nutilde
        else where
            production = SA_C_b1*(ONE - SA_c_t3)*vorticity*nutilde
        end where


        destruction = vorticity_mod
        where ( nutilde >= ZERO )
            destruction = (SA_c_w1*f_w - (SA_c_b1/(SA_kappa*SA_kappa))*f_t2) * (nutilde/dwall)**TWO
        else where
            destruction = -SA_c_w1 * (nutilde/dwall)**TWO
        end where


        !
        ! Compute jacobian of nutilde
        !
        dnutilde_drho         = -invrho*invrho*rho_nutilde
        dnutilde_drho_nutilde =  invrho


        !
        ! Compute gradient of nutilde
        !
        grad1_nutilde = dnutilde_drho*grad1_rho  +  dnutilde_drho_nutilde*grad1_rho_nutilde
        grad2_nutilde = dnutilde_drho*grad2_rho  +  dnutilde_drho_nutilde*grad2_rho_nutilde
        grad3_nutilde = dnutilde_drho*grad3_rho  +  dnutilde_drho_nutilde*grad3_rho_nutilde



        !========================================================================
        !                       Spalart-Allmaras Source Term
        !========================================================================
        source = -(                                 &
                    -rho*(production-destruction)   &
                    -(SA_c_b2/SA_sigma)*rho*(grad1_nutilde*grad1_nutilde + grad2_nutilde*grad2_nutilde + grad3_nutilde*grad3_nutilde)   &
                    +(ONE/SA_sigma)*(nu + f_n1*nutilde)*(grad1_rho*grad1_nutilde + grad2_rho*grad2_nutilde + grad3_rho*grad3_nutilde)   &
                  )

        call worker%integrate_volume('Density * NuTilde',source)


    end subroutine compute
    !*********************************************************************************************************






end module RANS_source
