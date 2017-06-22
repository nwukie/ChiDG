module type_partition
#include <messenger.h>
    use mod_kinds,                  only: ik
    use type_domain_connectivity,   only: domain_connectivity_t
    implicit none



    !>  A partition contains an array of connectivitites. Each is potentially
    !!  a subset of an entire connectivity for a domain, or it could be the
    !!  entire connectivity.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, public :: partition_t

        integer(ik)                                 :: nconn
        type(domain_connectivity_t),    allocatable :: connectivities(:)

    contains

        procedure   :: init

    end type  partition_t
    !*****************************************************************************




contains






    !>  Initialize the partition to have a connectivity instance for each domain
    !!  included in the partition.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!  @param[in]  nconn   Number of domain connectivities included in the partition
    !!
    !!
    !----------------------------------------------------------------------------
    subroutine init(self,nconn)
        class(partition_t), intent(inout)   :: self
        integer(ik),        intent(in)      :: nconn

        integer(ik) :: ierr

        self%nconn = nconn

        if (allocated(self%connectivities)) deallocate(self%connectivities)
        allocate(self%connectivities(nconn), stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !****************************************************************************







end module type_partition
