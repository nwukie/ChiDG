module DLA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none



    !> This equation set exists really just to test equationsets with more than one equation. 
    !! The idea is just to compute the linear advection solution twice at the same time. 
    !! The equations are independent of each other. So, we can verify, for example, the volume 
    !! flux jacobians for each equation. They should be the same as for the single 
    !! LinearAdvection equation set.
    !!
    !!  @author Nathan A. Wukie
    !!
    !---------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: DLA_boundary_average_advective_flux_t

    contains
        
        procedure   :: init
        procedure   :: compute

    end type DLA_boundary_average_advective_flux_t
    !*********************************************************************************************

contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(DLA_boundary_average_advective_flux_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("DLA Boundary Average Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field("u_a")
        call self%add_primary_field("u_b")

    end subroutine init
    !********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(DLA_boundary_average_advective_flux_t),   intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop

        type(AD_D), dimension(:), allocatable   ::  &
            ua_l, ua_r, ub_l, ub_r,                 &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p

        real(rk) :: c1, c2, c3


        !
        ! Get equation set properties
        !
        c1 = ONE
        c2 = ZERO
        c3 = ZERO


        !
        ! Interpolate solution to quadrature nodes
        !
        ua_r = worker%get_field('u_a', 'value', 'face interior')
        ub_r = worker%get_field('u_b', 'value', 'face interior')

        ua_l = worker%get_field('u_a', 'value', 'face exterior')
        ub_l = worker%get_field('u_b', 'value', 'face exterior')


        !
        ! Compute boundary average flux for u_a
        !
        flux_1_m = c1*ua_r
        flux_2_m = c2*ua_r
        flux_3_m = c3*ua_r

        flux_1_p = c1*ua_l
        flux_2_p = c2*ua_l
        flux_3_p = c3*ua_l

        call worker%integrate_boundary_average('u_a','Advection',               &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)



        !
        ! Compute boundary average flux for u_b
        !
        flux_1_m = c1*ub_r
        flux_2_m = c2*ub_r
        flux_3_m = c3*ub_r

        flux_1_p = c1*ub_l
        flux_2_p = c2*ub_l
        flux_3_p = c3*ub_l

        call worker%integrate_boundary_average('u_b','Advection',               &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !************************************************************************************







end module DLA_boundary_average_advective_flux
