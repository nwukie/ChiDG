!>
!! Description: This model computes the Reynolds stress tensor for the mean flow equations
!!              from the computed values of density*Rij.
!!
!! @author Eric M. Wolf
!! @date   01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_rstm_ssglrrw_realizable_stress
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use mod_rstm_ssglrrw

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_realizable_stress_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_realizable_stress_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric Wolf 
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_realizable_stress_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Realizable Reynolds Stress')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Reynolds-11')
        call self%add_model_field('Reynolds-22')
        call self%add_model_field('Reynolds-33')
        call self%add_model_field('Reynolds-12')
        call self%add_model_field('Reynolds-13')
        call self%add_model_field('Reynolds-23')



    end subroutine init
    !***************************************************************************************




    !>
    !! Description: Computes the Reynolds stress tensor
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_realizable_stress_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density,   &
            real_11, real_22, real_33, real_12, real_13, real_23, &
            rbar_12, rbar_13, rbar_23, det_R, &
            reynolds_11, reynolds_22, reynolds_33, &
            reynolds_12, reynolds_13, reynolds_23

        real(rk) :: det_R_pt
        real(rk), allocatable :: det_val(:)
        integer(ik) :: inode, maxiter, niter
        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        
        
        reynolds_11 = worker%get_field('Density * Reynolds-11',    'value')/density
        reynolds_22 = worker%get_field('Density * Reynolds-22',    'value')/density
        reynolds_33 = worker%get_field('Density * Reynolds-33',    'value')/density
        reynolds_12 = worker%get_field('Density * Reynolds-12',    'value')/density
        reynolds_13 = worker%get_field('Density * Reynolds-13',    'value')/density
        reynolds_23 = worker%get_field('Density * Reynolds-23',    'value')/density

        
        real_11 = reynolds_11
        real_22 = reynolds_22
        real_33 = reynolds_33
        rbar_12 = reynolds_12
        rbar_13 = reynolds_13
        rbar_23 = reynolds_23
        det_R   = ZERO*reynolds_11

        real_11 = reynolds_11*sin_ramp(reynolds_11, 0.0_rk, 0.1_rk*rstm_ssglrrw_R_infty) 
        real_22 = reynolds_22*sin_ramp(reynolds_22, 0.0_rk, 0.1_rk*rstm_ssglrrw_R_infty) 
        real_33 = reynolds_33*sin_ramp(reynolds_33, 0.0_rk, 0.1_rk*rstm_ssglrrw_R_infty) 

        rbar_12 = reynolds_12*sin_ramp(real_11*real_22-reynolds_12**TWO, 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO)
        rbar_13 = reynolds_13*sin_ramp(real_11*real_33-reynolds_13**TWO, 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO)
        rbar_23 = reynolds_23*sin_ramp(real_22*real_33-reynolds_23**TWO, 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO)

        !
        !do inode = 1, size(density)
        !    
        !    det_R_pt = real_11(inode)%x_ad_*(real_22(inode)%x_ad_*real_33(inode)%x_ad_-rbar_23(inode)%x_ad_*rbar_23(inode)%x_ad_)&
        !    - rbar_12(inode)%x_ad_*(rbar_12(inode)%x_ad_*real_33(inode)%x_ad_-rbar_23(inode)%x_ad_*rbar_13(inode)%x_ad_) &
        !    + rbar_13(inode)%x_ad_*(rbar_12(inode)%x_ad_*rbar_23(inode)%x_ad_-real_22(inode)%x_ad_*rbar_13(inode)%x_ad_)

        !    do while (det_R_pt < -1.0e-15_rk)
        !        rbar_12(inode) = 0.95_rk*rbar_12(inode)
        !        rbar_13(inode) = 0.95_rk*rbar_13(inode)
        !        rbar_23(inode) = 0.95_rk*rbar_23(inode)

        !        det_R_pt = real_11(inode)%x_ad_*(real_22(inode)%x_ad_*real_33(inode)%x_ad_-rbar_23(inode)%x_ad_*rbar_23(inode)%x_ad_)&
        !        - rbar_12(inode)%x_ad_*(rbar_12(inode)%x_ad_*real_33(inode)%x_ad_-rbar_23(inode)%x_ad_*rbar_13(inode)%x_ad_) &
        !        + rbar_13(inode)%x_ad_*(rbar_12(inode)%x_ad_*rbar_23(inode)%x_ad_-real_22(inode)%x_ad_*rbar_13(inode)%x_ad_)

        !    end do 

        !end do
        det_R = real_11*(real_22*real_33-rbar_23*rbar_23) - rbar_12*(rbar_12*real_33-rbar_23*rbar_13) + rbar_13*(rbar_12*rbar_23-real_22*rbar_13)
        !maxiter = 10_ik
        !niter = 0_ik
        !det_val = det_R%x_ad_
        !do while (any(det_val<-1.0e-15_rk))
        !    niter = niter + 1_ik
        !    rbar_12 = 0.95_rk*rbar_12
        !    rbar_13 = 0.95_rk*rbar_13
        !    rbar_23 = 0.95_rk*rbar_23
        !    if (niter .eq. maxiter) then
        !        rbar_12 = ZERO
        !        rbar_13 = ZERO
        !        rbar_23 = ZERO
        !    end if

        !    det_R = real_11*(real_22*real_33-rbar_23*rbar_23) - rbar_12*(rbar_12*real_33-rbar_23*rbar_13) + rbar_13*(rbar_12*rbar_23-real_22*rbar_13)
        !    det_val = det_R%x_ad_
        !end do

        !real_12 = rbar_12*sin_ramp(det_R, 0.0_rk, rstm_ssglrrw_k_infty)
        !real_13 = rbar_13*sin_ramp(det_R, 0.0_rk, rstm_ssglrrw_k_infty)
        !real_23 = rbar_23*sin_ramp(det_R, 0.0_rk, rstm_ssglrrw_k_infty)

        !real_11 = reynolds_11
        !real_22 = reynolds_22
        !real_33 = reynolds_33

        real_12 = reynolds_12
        real_13 = reynolds_13
        real_23 = reynolds_23

        !real_12 = HALF*(rbar_12+reynolds_12)
        !real_13 = HALF*(rbar_13+reynolds_13)
        !real_23 = HALF*(rbar_23+reynolds_23)

        call worker%store_model_field('Reynolds-11', 'value', real_11)
        call worker%store_model_field('Reynolds-22', 'value', real_22)
        call worker%store_model_field('Reynolds-33', 'value', real_33)
        call worker%store_model_field('Reynolds-12', 'value', real_12)
        call worker%store_model_field('Reynolds-13', 'value', real_13)
        call worker%store_model_field('Reynolds-23', 'value', real_23)




    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_realizable_stress
