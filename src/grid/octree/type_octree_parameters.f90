module type_octree_parameters
    use mod_kinds
    use mod_constants
    implicit none


    type, public :: octree_parameters_t
        integer(ik) :: bucket_size      = 32            ! Maximum number of nodes per leaf box
        real(rk)    :: min_extent       = ZERO          ! Minium extent of a leaf box
        integer(ik) :: refine_dir(3)    = 1             ! entry = 1 if octree is to be refined in the same direction
        logical     :: copy_points      = .false.       ! Copy array of pointst to local storage?


    contains

        procedure   :: init

    end type

contains

    subroutine init(self, bucket_size, min_extent, refine_dir, copy_points)
        class(octree_parameters_t),         intent(inout)                   ::  self
        integer(ik),                        intent(in), optional            ::  bucket_size
        real(rk),                           intent(in), optional            ::  min_extent
        integer(ik),                           intent(in), optional            ::  refine_dir(3)
        logical,                            intent(in), optional            :: copy_points                                                         


        if (present(bucket_size))   self%bucket_size    = bucket_size
        if (present(min_extent))    self%min_extent     = min_extent
        if (present(refine_dir))    self%refine_dir     = refine_dir
        if (present(copy_points))   self%copy_points    = copy_points 


    end subroutine init




end module type_octree_parameters
