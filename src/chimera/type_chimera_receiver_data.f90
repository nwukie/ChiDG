module type_chimera_receiver_data
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    use type_mvector,   only: mvector_t
    use type_pvector,   only: pvector_t




    !> Chimera receiver data container. One instance per Chimera face receiving data.
    !! Holds donor domain indices and donor element indices corresponding to their
    !! location in the donor domain.
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------
    type, public :: chimera_receiver_data_t

        integer(ik)                     :: ndonors = 0

        integer(ik)                     :: receiver_domain      !< Domain index of receiver
        integer(ik)                     :: receiver_element     !< Element index of receiver
        integer(ik)                     :: receiver_face        !< Face index of receiver



        ! Access via data%donor_domain%at(idonor)
        type(ivector_t)                 :: donor_domain         !< Vector of domain indices
        type(ivector_t)                 :: donor_element        !< Vector of element indices for the location in the corresponding domain
        type(mvector_t)                 :: donor_interpolator   !< Vector of matrices defining the Chimera interpolation

        ! The access for this component is slightly different than the above components
        ! Access via data%donor_gq_indices(idonor)%data()
        type(ivector_t), allocatable    :: donor_gq_indices(:)  !< Array of integer vectors defining the GQ node indices associated with a given donor
        type(pvector_t), allocatable    :: donor_coords(:)      !< Array of points definind the local coordinates of the GQ nodes
    contains

    end type chimera_receiver_data_t










end module type_chimera_receiver_data
