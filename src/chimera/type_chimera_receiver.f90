module type_chimera_receiver
#include <messenger.h>
    use mod_kinds,                      only: rk, ik
    use type_chimera_receiver_data,     only: chimera_receiver_data_t
    implicit none



    !>  A domain-level container for Chimera interface data.
    !!
    !!  Stores a chimera_receiver_data_t for each Chimera face in a given domain. This way,
    !!  a domain knows everything about the chimera communication that is going on. How many
    !!  Chimera faces, etc.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------------
    type, public :: chimera_receiver_t

        type(chimera_receiver_data_t),  allocatable :: data(:)

    contains

        procedure   :: init     ! Initialize the chimera receiver to allocated storage for nchimera_faces
        procedure   :: clear    ! Clear the Chimera data for each face
        procedure   :: nfaces   ! Return the number of chimera faces that currently have a 
                                ! chimera_receiver_data_t allocation

    end type chimera_receiver_t
    !************************************************************************************************




contains




    !>  Allocate a chimera_receiver_data_t instance for each Chimera face in a domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/11/2016
    !!
    !!  @param[in]  nchimera_faces  Integer indicating the number of chimera faces in a 
    !!                              domain. A chimera_receiver_data_t will be allocated for each.
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self,nchimera_faces)
        class(chimera_receiver_t),  intent(inout)   :: self
        integer(ik),                intent(in)      :: nchimera_faces

        integer(ik) :: ierr


        if (allocated(self%data)) deallocate(self%data)

        allocate(self%data(nchimera_faces), stat=ierr)
        if (ierr /= 0) call AllocationError
        

    end subroutine init
    !**********************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/11/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function nfaces(self) result(res)
        class(chimera_receiver_t),  intent(in)  :: self

        integer(ik)  :: res

        if (allocated(self%data)) then
            res = size(self%data)
        else
            res = 0
        end if

    end function nfaces
    !*********************************************************************************************







    !>  Call clear for each chimera_receiver_data_t allocation. This clears the data
    !!  for each chimera face, but the allocation for each face is still there. So 
    !!  those containers can be refilled with information.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/11/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chimera_receiver_t),  intent(inout)   :: self

        integer(ik) :: iface

        do iface = 1,self%nfaces()
            call self%data(iface)%clear()
        end do

    end subroutine clear
    !**************************************************************************************


end module type_chimera_receiver
