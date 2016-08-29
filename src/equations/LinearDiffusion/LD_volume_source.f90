module LD_volume_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,FOUR,PI

    use type_volume_flux,       only: volume_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use LD_properties,          only: LD_properties_t
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: LD_volume_source_t


    contains

        procedure   :: compute

    end type LD_volume_source_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/19/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(LD_volume_source_t),          intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        integer(ik)                             :: iu
        type(AD_D), allocatable, dimension(:)   :: source
        real(rk),   allocatable, dimension(:)   :: x

        !
        ! Get variable index from equation set
        !
        iu = prop%get_eqn_index("u")


        !
        ! Interpolate solution to quadrature nodes
        !
        source = worker%interpolate(iu, 'ddx')

        x = worker%x('volume')

        source = FOUR*PI*PI*dsin(TWO*PI*x)



        !
        ! Integrate volume flux
        !
        call worker%integrate_volume(iu, source)


    end subroutine compute
    !****************************************************************************************************






end module LD_volume_source
