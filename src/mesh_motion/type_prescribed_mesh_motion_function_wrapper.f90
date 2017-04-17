module type_prescribed_mesh_motion_function_wrapper
    use type_prescribed_mesh_motion_function,   only: prescribed_mesh_motion_function_t
    implicit none
    private

    !>  Wrapper for storing a polymorphic function type prescribed_mesh_motion_function_t
    !!      - This allows one to store an array of prescribed_mesh_motion_function_t. A work around for storing an array
    !!        of polymorphic entities
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !-------------------------------------------------------------
    type, public :: prescribed_mesh_motion_function_wrapper_t

        class(prescribed_mesh_motion_function_t), allocatable    :: pmmf

    end type prescribed_mesh_motion_function_wrapper_t
    !*************************************************************


end module type_prescribed_mesh_motion_function_wrapper
