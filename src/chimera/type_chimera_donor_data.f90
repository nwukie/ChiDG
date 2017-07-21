module type_chimera_donor_data
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    implicit none



    !>  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !!
    !---------------------------------------------------------------
    type, public :: chimera_donor_data_t

        integer(ik)                 :: donor_ID
        integer(ik)                 :: donor_proc

        ! For donor elements being communicated from off-processor, 
        ! their location in the recv container for accessing 
        ! solution data.
        integer(ik)                 :: donor_recv_comm
        integer(ik)                 :: donor_recv_domain
        integer(ik)                 :: donor_recv_element

        ! Interpolation node information
        integer(ik),    allocatable :: donor_gq_indices(:)    ! index in a node set the donor is providing data for
        real(rk),       allocatable :: donor_coords(:,:)      ! donor-local node coordinate(xi,eta,zeta)
        real(rk),       allocatable :: donor_metrics(:,:,:,:) ! 
        real(rk),       allocatable :: donor_jinv(:)

        ! Interpolators
        real(rk),       allocatable :: value(:,:)
        real(rk),       allocatable :: grad1(:,:)
        real(rk),       allocatable :: grad2(:,:)
        real(rk),       allocatable :: grad3(:,:)

    
    end type chimera_donor_data_t
    !***************************************************************









end module type_chimera_donor_data
