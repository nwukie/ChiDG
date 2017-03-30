module type_model_vorticity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use mod_fluid,          only: omega
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing vorticity
    !!
    !!  Model Fields:
    !!      : Vorticity-1  
    !!      : Vorticity-2
    !!      : Vorticity-3
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/23/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: model_vorticity_t

        real(rk)    :: gam = 1.4_rk     ! ratio of specific heats
        real(rk)    :: R   = 287.15_rk  ! ideal gas constant [J/(kg*K)]

    contains

        procedure   :: init
        procedure   :: compute

    end type model_vorticity_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/23/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(model_vorticity_t), intent(inout)   :: self

        call self%set_name('Vorticity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Vorticity-1')
        call self%add_model_field('Vorticity-2')
        call self%add_model_field('Vorticity-3')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(model_vorticity_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3,                      &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            vorticity_1,   vorticity_2,   vorticity_3,      &
            du_ddensity,   dv_ddensity,   dw_ddensity,      &
            du_dmom1,      dv_dmom2,      dw_dmom3,         &
            invdensity, v

        real(rk),   allocatable,    dimension(:) :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_primary_field_general('Density',    'value')
        mom1    = worker%get_primary_field_general('Momentum-1', 'value')
        mom2    = worker%get_primary_field_general('Momentum-2', 'value')
        mom3    = worker%get_primary_field_general('Momentum-3', 'value')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density = worker%get_primary_field_general('Density'   , 'grad1+lift')
        grad2_density = worker%get_primary_field_general('Density'   , 'grad2+lift')
        grad3_density = worker%get_primary_field_general('Density'   , 'grad3+lift')

        grad1_mom1    = worker%get_primary_field_general('Momentum-1', 'grad1+lift')
        grad2_mom1    = worker%get_primary_field_general('Momentum-1', 'grad2+lift')
        grad3_mom1    = worker%get_primary_field_general('Momentum-1', 'grad3+lift')

        grad1_mom2    = worker%get_primary_field_general('Momentum-2', 'grad1+lift')
        grad2_mom2    = worker%get_primary_field_general('Momentum-2', 'grad2+lift')
        grad3_mom2    = worker%get_primary_field_general('Momentum-2', 'grad3+lift')

        grad1_mom3    = worker%get_primary_field_general('Momentum-3', 'grad1+lift')
        grad2_mom3    = worker%get_primary_field_general('Momentum-3', 'grad2+lift')
        grad3_mom3    = worker%get_primary_field_general('Momentum-3', 'grad3+lift')


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
            mom2       = mom2/r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if



        !
        ! compute velocity jacobians
        !
        invdensity  = ONE/density
        du_ddensity = -invdensity*invdensity*mom1
        dv_ddensity = -invdensity*invdensity*mom2
        dw_ddensity = -invdensity*invdensity*mom3

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
            v = mom2/density

            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u + (v/r)) 


            !
            ! Account for rotation, convert to relative vorticity
            !
            vorticity_3 = vorticity_3 - TWO*omega

        end if



        call worker%store_model_field('Vorticity-1', 'value', vorticity_1)
        call worker%store_model_field('Vorticity-2', 'value', vorticity_2)
        call worker%store_model_field('Vorticity-3', 'value', vorticity_3)


    end subroutine compute
    !***************************************************************************************




end module type_model_vorticity
