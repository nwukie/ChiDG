module type_recv
    use mod_kinds,  only: ik
    implicit none


    !>  A container for a specifying a location within a parallel recv container.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !----------------------------------------------------------------------------
    type, public :: recv_t

        integer(ik) :: comm
        integer(ik) :: domain
        integer(ik) :: element

    end type recv_t
    !****************************************************************************






end module type_recv
