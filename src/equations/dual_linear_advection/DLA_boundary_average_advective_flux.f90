module DLA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                      ME, NEIGHBOR

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
        call self%set_equation("u_a")
        call self%set_equation("u_b")

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

        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu_a, iu_b

        type(AD_D), dimension(:), allocatable   ::  &
            ua_l, ua_r, ub_l, ub_r,                 &
            flux_x, flux_y, flux_z, integrand

        real(rk),   dimension(:), allocatable   ::  &
            normx, normy, normz


        !
        ! Get integer data
        !
        iu_a = prop%get_equation_index("u_a")
        iu_b = prop%get_equation_index("u_b")



!        !
!        ! Get equation set properties
!        !
!        select type(prop)
!            type is (DLA_properties_t)
!                cx = prop%c(1)
!                cy = prop%c(2)
!                cz = prop%c(3)
!        end select
        cx = 1._rk
        cy = 0._rk
        cz = 0._rk



        !
        ! Interpolate solution to quadrature nodes
        !
        ua_r = worker%interpolate(iu_a, 'value', ME)
        ua_l = worker%interpolate(iu_a, 'value', NEIGHBOR)

        ub_r = worker%interpolate(iu_b, 'value', ME)
        ub_l = worker%interpolate(iu_b, 'value', NEIGHBOR)


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)

        !
        ! Compute boundary average flux for u_a
        !
        flux_x = ((cx*ua_l + cx*ua_r)/TWO )  *  normx
        flux_y = ((cy*ua_l + cy*ua_r)/TWO )  *  normy
        flux_z = ((cz*ua_l + cz*ua_r)/TWO )  *  normz

        integrand = flux_x + flux_y + flux_z

        call worker%integrate_boundary(iu_a, integrand)



        !
        ! Compute boundary average flux for u_b
        !
        flux_x = ((cx*ub_l + cx*ub_r)/TWO )  *  normx
        flux_y = ((cy*ub_l + cy*ub_r)/TWO )  *  normy
        flux_z = ((cz*ub_l + cz*ub_r)/TWO )  *  normz

        integrand = flux_x + flux_y + flux_z

        call worker%integrate_boundary(iu_b, integrand)


    end subroutine compute
    !************************************************************************************







end module DLA_boundary_average_advective_flux
