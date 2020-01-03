module model_blended_reynolds_shear_stress
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use ieee_arithmetic,    only: ieee_is_nan
    use DNAD_D
    implicit none


    


    !>  A model for computing shear stress
    !!
    !!  Model Fields:
    !!      : shear_11, shear_12, shear_13
    !!      :           shear_22, shear_23
    !!      :                     shear_33
    !!
    !!  Lower-triangular components are not computed because the tensor is symmetric.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: blended_reynolds_shear_stress_t
        real(rk)    :: blend = 0.0_rk



    contains

        procedure   :: init
        procedure   :: compute

    end type blended_reynolds_shear_stress_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(blended_reynolds_shear_stress_t), intent(inout)   :: self
        real(rk)            :: blend 
        integer             :: unit, msg
        logical             :: file_exists

        namelist /blended_rstm/   blend 



        call self%set_name('Blended Reynolds Shear Stress')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Shear-11')
        call self%add_model_field('Shear-22')
        call self%add_model_field('Shear-33')
        call self%add_model_field('Shear-12')
        call self%add_model_field('Shear-13')
        call self%add_model_field('Shear-23')
        !
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%mu
        !   2: if not available, do nothing and mu retains default value
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=blended_rstm,iostat=msg)
            if (msg == 0) self%blend = blend 
            close(unit)
        end if


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(blended_reynolds_shear_stress_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3, energy,              &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_energy,  grad2_energy,  grad3_energy,     &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            mu,            mu_l,          mu_t,             &
            lamda,         lamda_l,       lamda_t,          &
            mu_a, lamda_a,                                 &
            shear_11,      shear_22,      shear_33,         &
            shear_12,      shear_13,      shear_23,         &
            du_ddensity,   dv_ddensity,   dw_ddensity,      &
            du_dmom1,      dv_dmom2,      dw_dmom3,         &
            reynolds_11,    reynolds_22, reynolds_33,       &
            reynolds_12,    reynolds_13, reynolds_23,       &
            invdensity, div_velocity, u, v, r

        integer(ik) :: ii, nnodes

        real(rk) :: blend

        blend = self%blend

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')
        energy  = worker%get_field('Energy',     'value')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density    = worker%get_field('Density'   , 'grad1')
        grad2_density    = worker%get_field('Density'   , 'grad2')
        grad3_density    = worker%get_field('Density'   , 'grad3')

        grad1_mom1       = worker%get_field('Momentum-1', 'grad1')
        grad2_mom1       = worker%get_field('Momentum-1', 'grad2')
        grad3_mom1       = worker%get_field('Momentum-1', 'grad3')

        grad1_mom2       = worker%get_field('Momentum-2', 'grad1')
        grad2_mom2       = worker%get_field('Momentum-2', 'grad2')
        grad3_mom2       = worker%get_field('Momentum-2', 'grad3')

        grad1_mom3       = worker%get_field('Momentum-3', 'grad1')
        grad2_mom3       = worker%get_field('Momentum-3', 'grad2')
        grad3_mom3       = worker%get_field('Momentum-3', 'grad3')

        grad1_energy     = worker%get_field('Energy    ', 'grad1')
        grad2_energy     = worker%get_field('Energy    ', 'grad2')
        grad3_energy     = worker%get_field('Energy    ', 'grad3')
        


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if



        !
        ! Get previously computed viscosities
        !
        mu_l    = worker%get_field('Laminar Viscosity',   'value')
        mu_t    = worker%get_field('Turbulent Viscosity', 'value')
        !mu_a    = worker%get_field('Smoothed Artificial Shear Viscosity', 'value')

        lamda_l = worker%get_field('Second Coefficient of Laminar Viscosity',   'value')
        lamda_t = worker%get_field('Second Coefficient of Turbulent Viscosity', 'value')
        !lamda_a = worker%get_field('Smoothed Artificial Bulk Viscosity', 'value')

        reynolds_11 = worker%get_field('Reynolds-11',   'value')
        reynolds_22 = worker%get_field('Reynolds-22',   'value')
        reynolds_33 = worker%get_field('Reynolds-33',   'value')
        reynolds_12 = worker%get_field('Reynolds-12',   'value')
        reynolds_13 = worker%get_field('Reynolds-13',   'value')
        reynolds_23 = worker%get_field('Reynolds-23',   'value')

        nnodes = size(reynolds_11)
        do ii=1, nnodes
            if (ieee_is_nan(reynolds_11(ii)%x_ad_)) print *, 'shear stress r_11 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_22(ii)%x_ad_)) print *, 'shear stress r_22 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_33(ii)%x_ad_)) print *, 'shear stress r_33 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_12(ii)%x_ad_)) print *, 'shear stress r_12 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_13(ii)%x_ad_)) print *, 'shear stress r_13 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_23(ii)%x_ad_)) print *, 'shear stress r_23 is nan', ' interp source : ', worker%interpolation_source
            if (ieee_is_nan(reynolds_23(ii)%x_ad_)) stop
                
            if (any(ieee_is_nan(reynolds_11(ii)%xp_ad_))) print *, 'shear stress r_11 deriv is nan'
            if (any(ieee_is_nan(reynolds_22(ii)%xp_ad_))) print *, 'shear stress r_22 deriv is nan'
            if (any(ieee_is_nan(reynolds_33(ii)%xp_ad_))) print *, 'shear stress r_33 deriv is nan'
            if (any(ieee_is_nan(reynolds_12(ii)%xp_ad_))) print *, 'shear stress r_12 deriv is nan'
            if (any(ieee_is_nan(reynolds_13(ii)%xp_ad_))) print *, 'shear stress r_13 deriv is nan'
            if (any(ieee_is_nan(reynolds_23(ii)%xp_ad_))) print *, 'shear stress r_23 deriv is nan'


        end do

        invdensity = ONE/density


        !
        ! Compute effective viscosities. Laminar + Turbulent
        !
        mu    = mu_l    + (1.0_rk-blend)*mu_t
        lamda = lamda_l + (1.0_rk-blend)*lamda_t
        !mu    = mu_l    + mu_t + mu_a
        !lamda = lamda_l + lamda_t + lamda_a

        !mu    = mu_l    !+ mu_t 
        !lamda = lamda_l !+ lamda_t 


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
        grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1
        grad3_u = du_ddensity*grad3_density  +  du_dmom1*grad3_mom1

        grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
        grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2
        grad3_v = dv_ddensity*grad3_density  +  dv_dmom2*grad3_mom2

        grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
        grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
        grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3








        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   div(V) = dudx + dvdy + dwdz
            !
            div_velocity = grad1_u + grad2_v + grad3_w

            !
            ! Compute shear stress components
            !
            shear_11 = TWO*mu*grad1_u  +  lamda*(div_velocity)  -   blend*density*reynolds_11
            shear_22 = TWO*mu*grad2_v  +  lamda*(div_velocity)  -   blend*density*reynolds_22
            shear_33 = TWO*mu*grad3_w  +  lamda*(div_velocity)  -   blend*density*reynolds_33

            shear_12 = mu*(grad2_u + grad1_v)   -   blend*density*reynolds_12
            shear_13 = mu*(grad3_u + grad1_w)   -   blend*density*reynolds_13
            shear_23 = mu*(grad2_w + grad3_v)   -   blend*density*reynolds_23




        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   div(U) = (1/r)d(r*u)dr   + (1/r)dv/dtheta + dw/dz
            !          = (du/dr + u/r) + (1/r)dv/dtheta + dw/dz
            !
            !   Note:
            !       grad1_u = du/dr
            !       grad2_v = (1/r)dv/dtheta
            !       grad3_w = dw/dz
            !
            u = mom1/density
            v = mom2/density
            div_velocity = (grad1_u + u/r) + grad2_v + grad3_w

            !
            ! Compute shear stress components
            !
            shear_11 = TWO*mu*(grad1_u        )  +  lamda*(div_velocity)
            shear_22 = TWO*mu*(grad2_v + (u/r))  +  lamda*(div_velocity)
            shear_33 = TWO*mu*(grad3_w        )  +  lamda*(div_velocity)

            shear_12 = mu*(grad2_u + grad1_v  - (v/r))
            shear_13 = mu*(grad3_u + grad1_w         )
            shear_23 = mu*(grad2_w + grad3_v         )

        end if



        call worker%store_model_field('Shear-11', 'value', shear_11)
        call worker%store_model_field('Shear-22', 'value', shear_22)
        call worker%store_model_field('Shear-33', 'value', shear_33)
        call worker%store_model_field('Shear-12', 'value', shear_12)
        call worker%store_model_field('Shear-13', 'value', shear_13)
        call worker%store_model_field('Shear-23', 'value', shear_23)


    end subroutine compute
    !***************************************************************************************




end module model_blended_reynolds_shear_stress
