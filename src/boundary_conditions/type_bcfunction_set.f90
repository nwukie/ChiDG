module type_bcfunction_set
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_bcfunction,    only: bcfunction_t
    implicit none







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !--------------------------------------------------------------------------------
    type, public :: bcfunction_set_t

        integer(ik)                         :: nfunctions_ = 0  !< Number of functions
        type(bcfunction_t),     allocatable :: bcfcn(:)         !< Bounday function list

    contains

        procedure   :: add      !< Procedure for adding functions to the list

    end type bcfunction_set_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine add(self,fname,ftype)
        class(bcfunction_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: fname
        character(*),               intent(in)      :: ftype

        type(bcfunction_t), allocatable :: temp_bcfcn(:)
        integer(ik)                     :: ierr, ifcn

        !
        ! Increment nfunctions
        !
        self%nfunctions_ = self%nfunctions_ + 1
        ifcn = self%nfunctions_


        !
        ! Resize array storage
        !
        allocate(temp_bcfcn(self%nfunctions_), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy previously initialized instances to new array.
        !
        if (self%nfunctions_ > 1) then
            temp_bcfcn(1:size(self%bcfcn)) = self%bcfcn(1:size(self%bcfcn))
        end if


        !
        ! Initialize new bcfunction
        !
        temp_bcfcn(ifcn)%name_      = fname
        temp_bcfcn(ifcn)%type_      = ftype
        !temp_bcfcn(ifcn)%status_    = 'Empty'


        !
        ! Move resized temp allocation back to self%bcfcn(:)
        !
        call move_alloc(temp_bcfcn,self%bcfcn)


    end subroutine add
    !********************************************************************************





end module type_bcfunction_set
