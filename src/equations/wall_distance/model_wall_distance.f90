module model_wall_distance
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO, RKTOL
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use ieee_arithmetic
!    use eqn_wall_distance,  only: get_p_poisson_parameter
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Wall Distance
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: wall_distance_m

    contains

        procedure   :: init
        procedure   :: compute

    end type wall_distance_m
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(wall_distance_m), intent(inout)   :: self

        call self%set_name('Wall Distance : p-Poisson Normalization')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Wall Distance')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Wall Distance using a normalization of the output
    !!  from a p-Poisson equation based method.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(wall_distance_m),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            d, grad1_d, grad2_d, grad3_d, d_normalization, sumsqr, rho

        real(rk) :: p
        integer(ik) :: igq


        ! Get primary field to initialize derivatives
        rho     = worker%get_field('Density', 'value')
        d       = rho
        grad1_d = rho
        grad2_d = rho
        grad3_d = rho


        !
        ! Interpolate solution to quadrature nodes
        !
        d       = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'value')
        grad1_d = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad1')
        grad2_d = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad2')
        grad3_d = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad3')


        !
        ! Compute wall distance normalization
        !
        !p = get_p_poisson_parameter()
        p = 6._rk
        sumsqr = grad1_d*grad1_d + grad2_d*grad2_d + grad3_d*grad3_d


        ! Beware of sumsqr==0, produces NaN, so don't normalize in this case. 
        ! Might happen if running P0 with a P1 wall distance?
        d_normalization = ZERO*d !allocate
        do igq = 1,size(d)
            ! Don't allow negative. Can't take fractional powers of negative
            ! numbers
            if (d_normalization(igq) < RKTOL) then
                d_normalization(igq) = ZERO
                
            else if (sumsqr(igq) < 1.e-8_rk) then
                d_normalization(igq) = d(igq)
                if (ieee_is_nan(d_normalization(igq)%x_ad_)) then
                    call write_line('d is nan', io_proc=GLOBAL_MASTER)
                    call write_line(d(igq)%x_ad_, grad1_d(igq)%x_ad_,grad2_d(igq)%x_ad_, grad3_d(igq)%x_ad_, io_proc=GLOBAL_MASTER)
                end if
            else
                d_normalization(igq) = (((p/(p-ONE))*d(igq)) + sumsqr(igq)**(p/TWO))**((p-ONE)/p) - sumsqr(igq)**((p-ONE)/TWO)
                if (ieee_is_nan(d_normalization(igq)%x_ad_)) then
                    call write_line('d_normalization is nan', io_proc=GLOBAL_MASTER)
                    call write_line(d(igq)%x_ad_, grad1_d(igq)%x_ad_,grad2_d(igq)%x_ad_, grad3_d(igq)%x_ad_, io_proc=GLOBAL_MASTER)
                end if
            end if
        end do






        !
        ! Store Wall Distance to model field
        !
        call worker%store_model_field('Wall Distance','value', d_normalization)


    end subroutine compute
    !***************************************************************************************




end module model_wall_distance
