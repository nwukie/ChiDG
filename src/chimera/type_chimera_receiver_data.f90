module type_chimera_receiver_data
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    use type_mvector,   only: mvector_t
    use type_pvector,   only: pvector_t
    use type_rvector,   only: rvector_t




    !> Chimera receiver data container. One instance per Chimera face receiving data.
    !! Holds donor domain indices and donor element indices corresponding to their
    !! location in the donor domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------------
    type, public :: chimera_receiver_data_t

        integer(ik)                     :: receiver_proc        !< Processor rank of receiver
        integer(ik)                     :: receiver_domain_g    !< ChiDG-global domain index of receiver
        integer(ik)                     :: receiver_domain_l    !< Proc-local domain index of receiver
        integer(ik)                     :: receiver_element_g   !< Domain-global element index of receiver
        integer(ik)                     :: receiver_element_l   !< Proc-local element index of receiver
        integer(ik)                     :: receiver_face        !< Face index of receiver


        ! Access via data%donor_domain%at(idonor)
        type(ivector_t)                 :: donor_neqns
        type(ivector_t)                 :: donor_nterms_s
        type(ivector_t)                 :: donor_proc               !< Vector of processor ranks
        type(ivector_t)                 :: donor_domain_g           !< Vector of domain indices
        type(ivector_t)                 :: donor_domain_l           !< Vector of domain indices
        type(ivector_t)                 :: donor_element_g          !< Vector of element indices for the location in the corresponding domain
        type(ivector_t)                 :: donor_element_l          !< Vector of element indices for the location in the corresponding domain
        type(mvector_t)                 :: donor_interpolator       !< Vector of matrices defining the Chimera interpolation
        type(mvector_t)                 :: donor_interpolator_ddx
        type(mvector_t)                 :: donor_interpolator_ddy
        type(mvector_t)                 :: donor_interpolator_ddz

        ! For donor elements being comminucated from off-processor, their location in the recv container for accessing components
        type(ivector_t)                 :: donor_recv_comm          !< Indices of comm container
        type(ivector_t)                 :: donor_recv_domain        !< Domain index within comm
        type(ivector_t)                 :: donor_recv_element       !< Element index within domain

        ! The access for this component is slightly different than the above components
        ! Access via data%donor_gq_indices(idonor)%data()
        type(ivector_t), allocatable    :: donor_gq_indices(:)      !< Array of integer vectors defining the GQ node indices associated with a given donor
        type(pvector_t), allocatable    :: donor_coords(:)          !< Array of points defining the local coordinates of the GQ nodes
        type(mvector_t), allocatable    :: donor_metrics(:)         !< For each donor, matrices of metric terms for each donor GQ node
        type(rvector_t), allocatable    :: donor_jinv(:)            !< For each donor, inverse element jacobian term for each donor GQ node

    contains

        procedure   :: ndonors

    end type chimera_receiver_data_t
    !***********************************************************************************************




contains






    !>  Return the number of donor elements providing the receiver with data.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/13/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    function ndonors(self) result(res)
        class(chimera_receiver_data_t), intent(in)  :: self

        integer(ik) :: res

        res = self%donor_neqns%size()

    end function ndonors
    !**************************************************************************************************






end module type_chimera_receiver_data
