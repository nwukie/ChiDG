module type_prescribed_mesh_motion_wrapper
    use type_prescribed_mesh_motion,   only: prescribed_mesh_motion_t
    implicit none
    private

    !>  Wrapper for storing a polymorphic function type prescribed_mesh_motion_t
    !!      - This allows one to store an array of prescribed_mesh_motion_t. A work around for storing an array
    !!        of polymorphic entities
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !-------------------------------------------------------------
    type, public :: prescribed_mesh_motion_wrapper_t

        class(prescribed_mesh_motion_t), allocatable    :: pmm

    end type prescribed_mesh_motion_wrapper_t
    !*************************************************************


end module type_prescribed_mesh_motion_wrapper
