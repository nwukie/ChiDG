module type_chimera_receiver
    use mod_kinds,                      only: rk, ik
    use type_chimera_receiver_data,     only: chimera_receiver_data_t



    !> A receiver for Chimera interface data.
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------------
    type, public :: chimera_receiver_t
        integer(ik)                                 :: nfaces

        type(chimera_receiver_data_t),  allocatable :: data(:)


!        integer(ik),            allocatable     :: map(:,:)     !< element index map.   'element index in donor domain => local chimera_receiver_t index in q(:)'
!        type(blockvector_t)                     :: q            !< blockvector_t for holding incoming solution data from external domains


    contains


    end type chimera_receiver_t
    !------------------------------------------------------------------------------------------------




contains





end module type_chimera_receiver
