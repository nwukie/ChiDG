module type_spalart_allmaras_turbulent_model_fields
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D

    use mod_spalart_allmaras,   only: SA_c_v1, SA_Pr_t
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Pressure
    !!      - Temperature
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: spalart_allmaras_turbulent_model_fields_t

        real(rk)    :: Cp = 1003._rk

    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_turbulent_model_fields_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(spalart_allmaras_turbulent_model_fields_t), intent(inout)   :: self

        call self%set_name('Spalart Allmaras Turbulent Model Fields')

        call self%add_model_field('Turbulent Viscosity')
        call self%add_model_field('Second Coefficient of Turbulent Viscosity')
        call self%add_model_field('Turbulent Thermal Conductivity')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the contributions to:
    !!      'Viscosity'
    !!      'Second Coefficient of Viscosity'
    !!      'Thermal Conductivity'
    !!
    !!  from the turbulent eddy viscosity.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(spalart_allmaras_turbulent_model_fields_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            rho, rhou, rhov, rhow, rhoE, mu, nu, rho_nutilde, mu_t, nutilde, lamda_t, &
            chi, f_v1, k_t


        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_primary_field_general('Density',           'value')
        rhou        = worker%get_primary_field_general('X-Momentum',        'value')
        rhov        = worker%get_primary_field_general('Y-Momentum',        'value')
        rhow        = worker%get_primary_field_general('Z-Momentum',        'value')
        rhoE        = worker%get_primary_field_general('Energy',            'value')
        rho_nutilde = worker%get_primary_field_general('Density * NuTilde', 'value')


        ! Get viscosity and second coefficient of viscosity
        mu = worker%get_model_field_general('Laminar Viscosity', 'value')


        ! Compute kinematic viscosity
        nu      = mu/rho
        nutilde = rho_nutilde/rho


        ! Compute chi
        chi = nutilde/nu


        ! Compute f_v1
        f_v1 = chi*chi*chi/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)


        ! Initialize derivatives, compute mu_t
        mu_t = rho
        where (nutilde >= 0)
            mu_t = rho * nutilde * f_v1
        else where
            mu_t = ZERO
        end where


        ! Compute Second Coefficient of Turbulent Viscosity, Stokes' Hypothesis.
        lamda_t = (-TWO/THREE)*mu_t

        ! Compute Turbulent Thermal Conductivity, Reynolds' analogy.
        k_t = self%Cp*mu_t/SA_Pr_t


        call worker%store_model_field('Turbulent Viscosity',                       'value', mu_t   )
        call worker%store_model_field('Second Coefficient of Turbulent Viscosity', 'value', lamda_t)
        call worker%store_model_field('Turbulent Thermal Conductivity',            'value', k_t    )


    end subroutine compute
    !***************************************************************************************




end module type_spalart_allmaras_turbulent_model_fields
