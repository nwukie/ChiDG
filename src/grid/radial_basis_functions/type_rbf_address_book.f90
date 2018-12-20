module type_rbf_address_book
    use type_rbf_address, only: rbf_address_t
    implicit none

    type :: rbf_address_book_t

        class(rbf_address_t), allocatable :: rbf_addresses(:)

    end type rbf_address_book_t


end module type_rbf_address_book
