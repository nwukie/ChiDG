module DLA_volume_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_volume_flux,       only: volume_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D

    use type_properties,        only: properties_t
    use DLA_properties,         only: DLA_properties_t
    implicit none

    private



    !> This equation set exists really just to test equationsets with more than one equation. 
    !! The idea is just to compute the linear advecdtion solution twice at the same time. 
    !! The equations are independent of each other. So, we can verify, for example, the volume 
    !! flux jacobians for each equation. They should be the same as for the single 
    !! LinearAdvection equation set.
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: DLA_volume_advective_flux_t

    contains

        procedure   :: compute

    end type DLA_volume_advective_flux_t
    !***********************************************************************************************

contains

    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(DLA_volume_advective_flux_t),     intent(in)      :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        real(rk)                :: cx, cy, cz
        integer(ik)             :: iu_a, iu_b

        type(AD_D), allocatable, dimension(:)   ::  &
            ua, ub, flux_x, flux_y, flux_z



        !
        ! Get variable index from equation set
        !
        iu_a = prop%get_eqn_index('u_a')
        iu_b = prop%get_eqn_index('u_b')


        !
        ! Get equation set properties
        !
        select type(prop)
            type is (DLA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select




        !
        ! Interpolate solution to quadrature nodes
        !
        ua = worker%interpolate(iu_a, 'value')
        ub = worker%interpolate(iu_b, 'value')



        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = cx  *  ua
        flux_y = cy  *  ua
        flux_z = cz  *  ua

        call worker%integrate_volume(iu_a, flux_x, flux_y, flux_z)



        ! Compute volume flux at quadrature nodes
        flux_x = cx  *  ub
        flux_y = cy  *  ub
        flux_z = cz  *  ub

        call worker%integrate_volume(iu_b, flux_x, flux_y, flux_z)




    end subroutine compute
    !***************************************************************************************************




end module DLA_volume_advective_flux
