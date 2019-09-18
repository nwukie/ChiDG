module type_rbf_mm_driver_wrapper
    use type_rbf_mm_driver,  only: rbf_mm_driver_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !-----------------------------------------------------------------------
    type, public :: rbf_mm_driver_wrapper_t

        class(rbf_mm_driver_t),  allocatable :: driver

    end type rbf_mm_driver_wrapper_t
    !************************************************************************






end module type_rbf_mm_driver_wrapper
