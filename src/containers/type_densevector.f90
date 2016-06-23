!> Data type for storing the dense block matrices for the linearization of each element
!!  @author Nathan A. Wukie
module type_densevector
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    implicit none






    !> Container for dense floating-point vector storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    type, public :: densevector_t

        ! Element Associativity
        integer(ik), private    :: parent_ = 0                  !< Associated parent element

        ! Storage size and equation information
        integer(ik), private    :: nterms_                      !< Number of terms in an expansion
        integer(ik), private    :: nvars_                       !< Number of equations included
    
        ! Vector storage
        real(rk),  dimension(:), allocatable :: vec             !< Vector storage

    contains

        ! Initializers
        procedure, public :: init           !< Initialize vector storage

        procedure :: parent                 !< return parent element
        procedure :: nentries               !< return number of vector entries
        procedure :: reparent               !< reassign parent
        procedure :: nterms                 !< return nterms_
        procedure :: nvars                  !< return nvars_


        procedure, public   :: setvar
        procedure, public   :: getvar

        procedure, public   :: getterm
        procedure, public   :: setterm


        procedure, public   :: clear


    end type densevector_t
    !*************************************************************************************************************











    !-------------------    OPERATORS   ---------------------
    public operator (*)
    interface operator (*)
        module procedure mult_real_dv   ! real * densevector,   ELEMENTAL
        module procedure mult_dv_real   ! densevector * real,   ELEMENTAL
    end interface



    public operator (/)
    interface operator (/)
        module procedure div_real_dv    ! real / densevector,   ELEMENTAL
        module procedure div_dv_real    ! densevector / real,   ELEMENTAL
    end interface



    public operator (+)
    interface operator (+)
        module procedure add_dv_dv      ! densevector + densevector,    ELEMENTAL
    end interface



    public operator (-)
    interface operator (-)
        module procedure sub_dv_dv      ! densevector - densevector,    ELEMENTAL
    end interface











    private
contains




    !> Subroutine for initializing dense-vector storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  nterms  Number of terms in an expansion
    !!  @param[in]  nvars   Number of equations being represented
    !!  @param[in]  parent  Index of associated parent element
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine init(self,nterms,nvars,parent)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms
        integer(ik),            intent(in)      :: nvars
        integer(ik),            intent(in)      :: parent

        integer(ik) :: ierr, vsize

        !
        ! Set dense-vector integer data
        !
        self%parent_ = parent
        self%nterms_ = nterms
        self%nvars_  = nvars


        !
        ! Compute total number of elements for densevector storage
        !
        vsize = nterms * nvars


        !
        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        !
        if (allocated(self%vec)) then
            deallocate(self%vec)
            allocate(self%vec(vsize), stat=ierr)
        else
            allocate(self%vec(vsize), stat=ierr)
        end if
        if (ierr /= 0) call AllocationError


        !
        ! Initialize to zero
        !
        self%vec = 0._rk

    end subroutine init
    !******************************************************************************************************












    !> Function returns the stored vector data associated with variable index ivar
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar        Integer index of the variable, for which modes will be returned
    !!  @result     modes_out   Array of modes from the variable, ivar
    !!
    !-------------------------------------------------------------------------------------------------------
    function getvar(self,ivar) result(modes_out)
        class(densevector_t),   intent(in)      :: self
        integer(ik),            intent(in)      :: ivar

        real(rk)                                :: modes_out(self%nterms_)
        integer(ik)                             :: istart, iend


        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + 1
        iend   = istart + (self%nterms_-1)


        !
        ! Return modes
        !
        modes_out = self%vec(istart:iend)

        ! ifort has occasional and inconsistent trouble remapping vec to mat so this access is sometimes wrong.
        ! Pretty difficult to diagnose. Should test as much as possible before reenabling the map
        !modes_out = self%mat(:,ivar)
    end function getvar
    !******************************************************************************************************














    !>  Set the modes for a particular variable
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar    Integer index of the variable being set
    !!  @param[in]  vals    Array of mode values that will be set
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine setvar(self,ivar,vals)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: ivar
        real(rk),               intent(in)      :: vals(:)

        integer(ik) :: istart, iend

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + 1
        iend   = istart + (self%nterms_-1)


        !
        ! Set modes
        !
        self%vec(istart:iend) = vals


    end subroutine setvar
    !******************************************************************************************************







    !>  Return an individual mode from the polynomial expansion of a variable
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar    Integer index of the variable
    !!  @param[in]  iterm   Integer index of the mode in the expansion to be returned
    !!
    !-------------------------------------------------------------------------------------------------------
    function getterm(self,ivar,iterm) result(mode_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: ivar
        integer(ik),            intent(in)  :: iterm

        real(rk)    :: mode_out
        integer(ik) :: istart, iterm_g

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + 1
        iterm_g = istart + (iterm-1)
    

        !
        ! Get mode
        !
        mode_out = self%vec(iterm_g)

    end function getterm
    !******************************************************************************************************








    !>  Set an individual mode in the polynomial expansion of a variable
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar        Integer index of the variable
    !!  @param[in]  iterm       Integer index of the mode in the expansion to be set
    !!  @param[in]  mode_in     Floating point value, which is the mode amplitude to be set
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine setterm(self,ivar,iterm,mode_in)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: ivar, iterm
        real(rk),               intent(in)      :: mode_in

        integer(ik) :: istart, iterm_g

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + 1
        iterm_g = istart + (iterm-1)


        !
        ! Set mode
        !
        self%vec(iterm_g) = mode_in

    end subroutine setterm
    !******************************************************************************************************











    !> Function that returns number of entries in block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------
    function nentries(self) result(n)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: n

        n = size(self%vec)

    end function nentries
    !******************************************************************************************************












    !> Function that returns nterms_ private component
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------
    pure function nterms(self) result(nterms_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik)                         :: nterms_out

        nterms_out = self%nterms_

    end function nterms
    !******************************************************************************************************










    !> Function that returns nvars_ private component
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------
    pure function nvars(self) result(nvars_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik)                         :: nvars_out

        nvars_out = self%nvars_

    end function nvars
    !******************************************************************************************************














    !> Function that returns index of block parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------
    function parent(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%parent_

    end function parent
    !******************************************************************************************************












    !> reset index of parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  par     Index of new parent element
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine reparent(self,par)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: par


        !
        ! Set element parent index
        !
        self%parent_ = par

    end subroutine reparent
    !******************************************************************************************************












    !> Zero vector storage, self%vec
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(densevector_t),   intent(inout)   :: self


        self%vec = ZERO

    end subroutine clear
    !******************************************************************************************************

















    !-------------------------------------------------------------------
    !-------------      OPERATOR IMPLEMENTATIONS    --------------------
    !-------------------------------------------------------------------
    elemental function mult_real_dv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res

        res%parent_ = right%parent_
        res%nvars_  = right%nvars_
        res%nterms_ = right%nterms_

        res%vec     = left * right%vec

    end function

   
    elemental function mult_dv_real(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right
        type(densevector_t) :: res


        res%parent_ = left%parent_
        res%nvars_  = left%nvars_
        res%nterms_ = left%nterms_

        res%vec     = left%vec * right

    end function








    elemental function div_real_dv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%parent_ = right%parent_
        res%nvars_  = right%nvars_
        res%nterms_ = right%nterms_
        
        res%vec     = left / right%vec

    end function


    elemental function div_dv_real(left,right) result(res)
        type(densevector_t),        intent(in)  :: left
        real(rk),                   intent(in)  :: right
        type(densevector_t) :: res


        res%parent_ = left%parent_
        res%nvars_  = left%nvars_
        res%nterms_ = left%nterms_

        res%vec     = left%vec / right

    end function







    elemental function add_dv_dv(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%parent_ = left%parent_
        res%nvars_  = left%nvars_
        res%nterms_ = left%nterms_

        res%vec     = left%vec + right%vec

    end function




    elemental function sub_dv_dv(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%parent_ = left%parent_
        res%nvars_  = left%nvars_
        res%nterms_ = left%nterms_

        res%vec     = left%vec - right%vec

    end function


    












end module type_densevector
