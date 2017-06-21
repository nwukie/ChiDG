module type_boundary_connectivity
#include <messenger.h>
    use mod_kinds,              only: ik
    use type_face_connectivity, only: face_connectivity_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: boundary_connectivity_t

        type(face_connectivity_t),  allocatable :: data(:)

    contains
        
        procedure   :: init
        procedure   :: nfaces            ! Return number of face connectivities contained

    end type boundary_connectivity_t
    !*******************************************************************************


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!  @param[in]  nfaces  Number of faces in the connectivity
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self,nfaces)
        class(boundary_connectivity_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: nfaces
        
        integer(ik) :: ierr

        allocate(self%data(nfaces), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !********************************************************************************








    !>  Return the number of faces in the connectivity structure.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function nfaces(self) result(nfaces_)
        class(boundary_connectivity_t),  intent(in)  :: self

        integer(ik) :: nfaces_

        nfaces_ = size(self%data)

    end function nfaces
    !********************************************************************************









end module type_boundary_connectivity
