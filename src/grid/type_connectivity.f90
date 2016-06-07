module type_connectivity
    use mod_kinds,  only: ik
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: connectivity_t

        integer(ik),    allocatable :: data(:,:)

    contains

    end type connectivity_t
    !*******************************************************************************


contains







end module type_connectivity
