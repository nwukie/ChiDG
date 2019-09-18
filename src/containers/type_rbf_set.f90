module type_rbf_set
    implicit none

    type, public :: rbf_set_t

        type(ivector_t) :: registered_rbf_indices

    end type rbf_set_t

end module type_rbf_set
