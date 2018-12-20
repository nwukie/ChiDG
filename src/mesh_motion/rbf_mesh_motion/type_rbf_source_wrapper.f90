module type_rbf_source_wrapper
    use type_rbf_source,  only: rbf_source_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !-----------------------------------------------------------------------
    type, public :: rbf_source_wrapper_t

        class(rbf_source_t),  allocatable :: source

    end type rbf_source_wrapper_t
    !************************************************************************






end module type_rbf_source_wrapper
