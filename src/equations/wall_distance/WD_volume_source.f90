module WD_volume_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,FOUR,PI

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
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
    type, extends(operator_t), public :: WD_volume_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type WD_volume_source_t
    !*************************************************************************

contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(WD_volume_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Wall Distance Volume Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Source")

        ! Set operator equations
        call self%add_primary_field("u")

    end subroutine init
    !********************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(WD_volume_source_t),          intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        integer(ik)                             :: iu
        type(AD_D), allocatable, dimension(:)   :: source
        real(rk),   allocatable, dimension(:)   :: x, y

        !
        ! Get variable index from equation set
        !
        iu = prop%get_primary_field_index("u")


        !
        ! Interpolate solution to quadrature nodes to initialize derivatives
        !
        source = worker%get_primary_field_element('u',iu, "value")


        source = ONE


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('u',iu, source)


    end subroutine compute
    !****************************************************************************************************






end module WD_volume_source
