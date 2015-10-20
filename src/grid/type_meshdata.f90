module type_meshdata
    use mod_kinds,  only: ik
    use type_point, only: point_t



    !> Data type for returning mesh-data from a file-read routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !----------------------------------------------------------------------
    type, public :: meshdata_t

        character(len=100)          :: name
        type(point_t), allocatable  :: points(:,:,:)    !< Rank-3 array containing mesh points
        integer(ik)                 :: nterms_c         !< Integer specifying the number of terms in the coordinate expansion

    end type meshdata_t







end module type_meshdata
