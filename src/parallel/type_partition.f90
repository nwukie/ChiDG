module type_partition
    use mod_kinds,  only: ik
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !-------------------------------------------------------------------
    type, public :: partition_t

        integer(ik)                 :: ipartition
        integer(ik),    allocatable :: data(:,:)

    contains

    end type  partition_t
    !*******************************************************************




contains







end module type_partition
