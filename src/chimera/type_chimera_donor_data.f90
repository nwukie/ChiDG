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

        ! Donor information
        integer(ik)                 :: donor_ID             ! ID in chimera%donors(donor_ID)
        integer(ik)                 :: donor_proc           ! host processor rank
        integer(ik)                 :: donor_recv_comm      ! location of donor solution
        integer(ik)                 :: donor_recv_domain    ! location of donor solution
        integer(ik)                 :: donor_recv_element   ! location of donor solution

        ! Node information
        integer(ik),    allocatable :: donor_gq_indices(:)    ! index in a node set the donor is providing data for
        real(rk),       allocatable :: donor_coords(:,:)      ! donor-local node coordinate(xi,eta,zeta)
        real(rk),       allocatable :: donor_metrics(:,:,:,:) ! 
        real(rk),       allocatable :: donor_jinv(:)

        ! Interpolators
        real(rk),       allocatable :: value(:,:)
        real(rk),       allocatable :: grad1(:,:)
        real(rk),       allocatable :: grad2(:,:)
        real(rk),       allocatable :: grad3(:,:)


    contains

        procedure   :: add_node
        procedure   :: nnodes

    
    end type chimera_donor_data_t
    !***************************************************************


contains




    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !----------------------------------------------------------------
    subroutine add_node(self,inode,coord)
        class(chimera_donor_data_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: inode
        real(rk),                       intent(in)      :: coord(3)







    end subroutine add_node
    !****************************************************************


    !>
    !!
    !!
    !!
    !----------------------------------------------------------------
    function nnodes(self) result(n)
        class(chimera_donor_data_t),    intent(in)  :: self

        integer(ik) :: n



    end function nnodes
    !****************************************************************




end module type_chimera_donor_data
