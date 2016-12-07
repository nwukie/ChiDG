module type_model_wrapper
    use type_model, only: model_t
    implicit none



    !>  A wrapper enabling arrays of abstract model_t classes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !--------------------------------------------------------------------
    type, public :: model_wrapper_t

        class(model_t), allocatable :: model

    end type model_wrapper_t
    !********************************************************************


end module type_model_wrapper
