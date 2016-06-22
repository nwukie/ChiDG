module type_chimera_receiver
    use mod_kinds,                      only: rk, ik
    use type_chimera_receiver_data,     only: chimera_receiver_data_t



    !>  A receiver for Chimera interface data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------------
    type, public :: chimera_receiver_t
        integer(ik)                                 :: nfaces

        type(chimera_receiver_data_t),  allocatable :: data(:)

    contains


    end type chimera_receiver_t
    !************************************************************************************************




contains





end module type_chimera_receiver
