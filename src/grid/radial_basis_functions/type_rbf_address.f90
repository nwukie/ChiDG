module type_rbf_address
    use mod_kinds
    use type_ivector, only: ivector_t
    implicit none

    type :: rbf_address_t
        
        integer(ik)     :: rbf_set_ID
        type(ivector_t) :: registered_rbf_indices


    end type rbf_address_t


end module type_rbf_address
