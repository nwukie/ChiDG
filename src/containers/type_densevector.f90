!> Data type for storing the dense block matrices for the linearization of each element
!!  @author Nathan A. Wukie
module type_densevector
#include <messenger.h>
    use mod_kinds,  only: rk,ik
    implicit none

    type, public :: densevector_t
        !> Element Associativity
        integer(ik), private    :: parent_ = 0 !> parent element

        !> Vector storage
        real(rk),  dimension(:), allocatable :: vec     !> Vector storage
        real(rk),  dimension(:,:), pointer   :: mat     !> Matrix-view alias of vec  (nterms, neq)

    contains
        !> Initializers
        generic, public :: init => init_vector
        procedure, private :: init_vector       !> Initialize vector storage

        procedure :: parent     !> return parent element
        procedure :: nentries   !> return number of matrix entries
        procedure :: resize     !> resize vector storage
        procedure :: reparent   !> reassign parent

        final :: destructor
    end type densevector_t



    private
contains

    !> Subroutine for initializing square dense-vector storage
    !-----------------------------------------------------------
    subroutine init_vector(self,vsize,parent)
        class(densevector_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: vsize, parent

        integer(ik) :: ierr

        ! Block parents
        self%parent_ = parent

        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        if (allocated(self%vec)) then
            deallocate(self%vec)
            allocate(self%vec(vsize),stat=ierr)
        else
            allocate(self%vec(vsize),stat=ierr)
        end if
        if (ierr /= 0) call AllocationError

        ! Initialize to zero
        self%vec = 0._rk
    end subroutine




    !> return number of entries in block storage
    !------------------------------------------------------------
    function nentries(self) result(n)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: n

        n = size(self%vec)
    end function

    !> return index of block parent
    !------------------------------------------------------------
    function parent(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%parent_
    end function


    !> Resize dense-block storage
    !------------------------------------------------------------
    subroutine resize(self,vsize)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: vsize

        integer(ik) :: ierr

        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        if (allocated(self%vec)) then
            deallocate(self%vec)
            allocate(self%vec(vsize),stat=ierr)
        else
            allocate(self%vec(vsize),stat=ierr)
        end if
        if (ierr /= 0) call AllocationError

    end subroutine


    !> set index of parent
    !------------------------------------------------------------
    subroutine reparent(self,par)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: par

        self%parent_ = par
    end subroutine


    subroutine destructor(self)
        type(densevector_t),    intent(inout)   :: self

    end subroutine

end module type_densevector
