module type_chimera_receiver_data
    use type_ivector,   only: ivector_t




    !> Chimera receiver data container. One instance per Chimera face receiving data.
    !! Holds donor domain indices and donor element indices corresponding to their
    !! location in the donor domain.
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------
    type, public :: chimera_receiver_data_t

        type(ivector_t)                 :: donor_domain     !< Vector of domain indices
        type(ivector_t)                 :: donor_element    !< Vector of element indices for the location in the corresponding domain
    contains

    end type chimera_receiver_data_t










end module type_chimera_receiver_data
