module rans_volume_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use mod_spalart_allmaras,   only: SA_c_n1, SA_c_w1, SA_c_b1, SA_c_b2, SA_c_t3,  &
                                      SA_kappa, SA_sigma, SA_b, SA_c_t4, SA_c_v1,   &
                                      SA_c_v2, SA_c_w2, SA_c_v3, SA_c_w3, SA_rlim
    use mod_rans_efficient

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private

    
    !>  Volume source terms from viscous fluxes in cylindrical coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/22/2017
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: rans_volume_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rans_volume_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/22/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_volume_source_t),   intent(inout)      :: self

        integer         :: unit, msg
        logical         :: file_exists

!        namelist /fluid/ turbulence_model




        ! Set operator name
        call self%set_name('RANS Volume Source')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )







        ! Check if input from 'models.nml' is available.
        !   1: if available, read 
        !   2: if not available, do nothing and turbulence_model retains default value
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=fluid,iostat=msg)
            close(unit)
        end if



        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                call self%add_primary_field('Density * NuTilde')
                call self%add_auxiliary_field('Wall Distance : p-Poisson')
                call self%add_model('Wall Distance : p-Poisson Normalization')
            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select


    end subroutine init
    !********************************************************************************



    !> Volume flux routine for viscous equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rans_volume_source_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        ! Inviscid/Viscous
        type(AD_D), allocatable, dimension(:) :: tau_22, source 
        real(rk),   allocatable               :: ale_g(:), ale_Dinv(:,:,:)

        ! Spalart Allmaras
        type(AD_D), allocatable, dimension(:)   ::                  &
            density, grad1_density, grad2_density, grad3_density,   &
            mom1,    grad1_mom1,    grad2_mom1,    grad3_mom1,      &
            mom2,    grad1_mom2,    grad2_mom2,    grad3_mom2,      &
            mom3,    grad1_mom3,    grad2_mom3,    grad3_mom3,      &
            energy,  grad1_energy,  grad2_energy,  grad3_energy,    &
            density_nutilde,    grad1_density_nutilde,    grad2_density_nutilde,    grad3_density_nutilde,      &
            vorticity_1, vorticity_2, vorticity_3,                  &
            mu_l, mu_t, lambda_l, lambda_t, nu, nutilde,            &
            dwall, chi, f_v1, f_n1, vorticity2, vorticity, vorticity_bar, vorticity_mod, r_sa, rbar, g, f_w, f_t2, &
            destruction, production, dnutilde_ddensity, dnutilde_ddensity_nutilde, grad1_nutilde, grad2_nutilde, grad3_nutilde, &
            f_v2, invdensity, T, p, r

        real(rk)    :: eps, epsilon_vorticity



        !=================================================
        ! mass flux
        !=================================================


        !=================================================
        ! momentum-1 flux
        !=================================================

        ! Source term due to cylindrical coordinate system
        if (worker%coordinate_system() == 'Cylindrical') then

            ! Get grid data
            r = worker%coordinate('1','volume')
            ale_g    = worker%get_det_jacobian_grid_element('value')
            ale_Dinv = worker%get_inv_jacobian_grid_element()


            ! get shear stress
            tau_22 = worker%get_field('Shear-22','value','element')


            ! Compute/integrate source term
            !source = -tau_22 / r
            source = -ale_g*ale_Dinv(2,2,:)*tau_22 / r

            call worker%integrate_volume_source('Momentum-1',source)

        end if


        !=================================================
        ! momentum-2 flux
        !=================================================


        !=================================================
        ! momentum-3 flux
        !=================================================


        !=================================================
        ! energy flux
        !=================================================



        !=================================================
        ! Spalart Allmaras Source
        !=================================================


        select case (trim(turbulence_model))
            case('Spalart Allmaras')

                ! Interpolate solution gradients to quadrature nodes
                density       = worker%get_field('Density', 'value', 'element')
                grad1_density = worker%get_field('Density', 'grad1', 'element')
                grad2_density = worker%get_field('Density', 'grad2', 'element')
                grad3_density = worker%get_field('Density', 'grad3', 'element')

                mom1       = worker%get_field('Momentum-1', 'value', 'element')
                grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'element')
                grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'element')
                grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'element')

                mom2       = worker%get_field('Momentum-2', 'value', 'element')
                grad1_mom2 = worker%get_field('Momentum-2', 'grad1', 'element')
                grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'element')
                grad3_mom2 = worker%get_field('Momentum-2', 'grad3', 'element')

                mom3       = worker%get_field('Momentum-3', 'value', 'element')
                grad1_mom3 = worker%get_field('Momentum-3', 'grad1', 'element')
                grad2_mom3 = worker%get_field('Momentum-3', 'grad2', 'element')
                grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'element')

                energy       = worker%get_field('Energy', 'value', 'element')
                grad1_energy = worker%get_field('Energy', 'grad1', 'element')
                grad2_energy = worker%get_field('Energy', 'grad2', 'element')
                grad3_energy = worker%get_field('Energy', 'grad3', 'element')

                density_nutilde       = worker%get_field('Density * NuTilde', 'value', 'element')
                grad1_density_nutilde = worker%get_field('Density * NuTilde', 'grad1', 'element')
                grad2_density_nutilde = worker%get_field('Density * NuTilde', 'grad2', 'element')
                grad3_density_nutilde = worker%get_field('Density * NuTilde', 'grad3', 'element')


                ! Compute pressure and temperature
                call compute_pressure_temperature(density,mom1,mom2,mom3,energy,p,T)



                ! Interpolate auxiliary field, Wall Distance
                eps = 1.e-11_rk
                dwall = worker%get_field('Wall Distance', 'value', 'element')
                if (any(ieee_is_nan(dwall(:)%x_ad_))) call write_line('dwall is nan',io_proc=GLOBAL_MASTER)


                ! Divide by density
                invdensity  = ONE/density
                nutilde = density_nutilde*invdensity


                ! Compute model values
                call compute_viscosity(density,T,density_nutilde,mu_l, mu_t, lambda_l, lambda_t)
                nu = mu_l*invdensity


                ! Compute turbulence viscosity
                chi  = nutilde/nu
                f_v1 = (chi*chi*chi)/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)
                

                ! Compute f_n1
                f_n1 = density_nutilde
                f_n1 = ONE
                where(nutilde < ZERO)
                    f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
                end where


                ! Compute vorticity and modified vorticity
                call compute_vorticity(worker%coordinate_system(),worker%coordinate('1','volume'),                                &
                                       density,grad1_density,grad2_density,grad3_density,  &
                                       mom1,   grad1_mom1,   grad2_mom1,   grad3_mom1,     &
                                       mom2,   grad1_mom2,   grad2_mom2,   grad3_mom2,     &
                                       mom3,   grad1_mom3,   grad2_mom3,   grad3_mom3,     &
                                       vorticity_1, vorticity_2, vorticity_3)



                vorticity2 =  vorticity_1**TWO  +  vorticity_2**TWO  +  vorticity_3**TWO 
                epsilon_vorticity = 1.e-6_rk
                vorticity = vorticity2
                where(vorticity2 < epsilon_vorticity)
                    vorticity = HALF*(epsilon_vorticity + vorticity2/epsilon_vorticity)
                else where
                    vorticity = sqrt(vorticity2)
                end where



                f_v2 = ONE - (chi/(ONE+chi*f_v1))
                vorticity_bar = (nutilde/(SA_kappa*SA_kappa*(dwall*dwall + eps)))*f_v2


                vorticity_mod = vorticity
                where (vorticity_bar >= -SA_c_v2*vorticity)
                    vorticity_mod = vorticity + vorticity_bar
                else where
                    vorticity_mod = vorticity + vorticity*(SA_c_v2*SA_c_v2*vorticity + SA_c_v3*vorticity_bar)/( (SA_c_v3 - TWO*SA_c_v2)*vorticity - vorticity_bar ) 
                end where


                !
                ! Compute f_t2, f_w, g, r
                !
                rbar = nutilde/(vorticity_mod * SA_kappa * SA_kappa * (dwall*dwall + eps) )
                r_sa = SA_rlim - (SA_rlim - rbar)*(atan(SA_b * (SA_rlim - rbar))/PI  + HALF)  + atan(SA_b)/PI  -  HALF
                g = r_sa + SA_c_w2*(r_sa**SIX  -  r_sa)
                f_w = g*(((ONE + SA_c_w3**SIX)/(g**SIX + SA_c_w3**SIX))**(ONE/SIX))
                f_t2 = SA_c_t3*exp(-SA_c_t4*chi*chi)



                !
                ! Compute Production, Destruction
                !
                production = vorticity_mod
                where ( nutilde >= ZERO )
                    production = SA_c_b1*(ONE - f_t2)*vorticity_mod*nutilde
                else where
                    production = SA_c_b1*(ONE - SA_c_t3)*vorticity*nutilde
                end where


                destruction = vorticity_mod
                where ( nutilde >= ZERO )
                    destruction = (SA_c_w1*f_w - (SA_c_b1/(SA_kappa*SA_kappa))*f_t2) * (nutilde*nutilde/(dwall*dwall + eps))
                else where
                    destruction = -SA_c_w1 * (nutilde*nutilde/(dwall*dwall + eps))
                end where


                !
                ! Compute jacobian of nutilde
                !
                dnutilde_ddensity         = -invdensity*invdensity*density_nutilde
                dnutilde_ddensity_nutilde =  invdensity


                !
                ! Compute gradient of nutilde
                !
                grad1_nutilde = dnutilde_ddensity*grad1_density  +  dnutilde_ddensity_nutilde*grad1_density_nutilde
                grad2_nutilde = dnutilde_ddensity*grad2_density  +  dnutilde_ddensity_nutilde*grad2_density_nutilde
                grad3_nutilde = dnutilde_ddensity*grad3_density  +  dnutilde_ddensity_nutilde*grad3_density_nutilde


                !========================================================================
                !                       Spalart-Allmaras Source Term
                !========================================================================
                source = -(                                     &
                            -density*(production-destruction)   &
                            -(SA_c_b2/SA_sigma)*density*(grad1_nutilde*grad1_nutilde + grad2_nutilde*grad2_nutilde + grad3_nutilde*grad3_nutilde)   &
                            +(ONE/SA_sigma)*(nu + f_n1*nutilde)*(grad1_density*grad1_nutilde + grad2_density*grad2_nutilde + grad3_density*grad3_nutilde)   &
                          )

                call worker%integrate_volume_source('Density * NuTilde',source)


            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select







    end subroutine compute
    !*********************************************************************************************************






end module rans_volume_source
