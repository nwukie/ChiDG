module type_densevector
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, TWO
    implicit none






    !>  Container for dense floating-point vector storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: densevector_t

        ! Element Associativity
        integer(ik)             :: dparent_g_
        integer(ik)             :: dparent_l_
        integer(ik)             :: eparent_g_
        integer(ik)             :: eparent_l_

        ! Storage size and equation information
        integer(ik), private    :: nterms_              ! Number of terms in an expansion
        integer(ik), private    :: nvars_               ! Number of equations included
        integer(ik), private    :: ntime_               ! Number of time instances stored

        ! Vector storage
        real(rk),  dimension(:), allocatable :: vec     ! Vector storage

    contains

        procedure, public :: init           ! Initialize vector storage

        procedure, public :: dparent_g
        procedure, public :: dparent_l
        procedure, public :: eparent_g
        procedure, public :: eparent_l
        procedure, public :: nentries       ! return number of vector entries
        procedure, public :: nterms         ! return nterms_
        procedure, public :: nvars          ! return nvars_
        procedure, public :: ntime          ! return ntime_
        procedure, public :: set_ntime      ! procedure for changing ntime 

        procedure, public :: settime
        procedure, public :: gettime
        procedure, public :: setvar
        procedure, public :: getvar
        procedure, public :: setterm
        procedure, public :: getterm

        procedure, public :: get_time_start
        procedure, public :: get_time_end

        procedure, public :: sumsqr
        procedure, public :: sumsqr_fields

        procedure, public :: clear

        procedure, public :: restrict
        procedure, public :: prolong

    end type densevector_t
    !******************************************************************************************











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




    !>  Subroutine for initializing dense-vector storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  nterms  Number of terms in an expansion
    !!  @param[in]  nvars   Number of equations being represented
    !!  @param[in]  parent  Index of associated parent element
    !!  
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!` @date   11/15/2016
    !!
    !!  @param[in]  ntime   Number of time levels in solution
    !! 
    !------------------------------------------------------------------------------------------
    subroutine init(self,nterms,nvars,ntime,dparent_g,dparent_l,eparent_g,eparent_l)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms
        integer(ik),            intent(in)      :: nvars
        integer(ik),            intent(in)      :: ntime
        integer(ik),            intent(in)      :: dparent_g
        integer(ik),            intent(in)      :: dparent_l
        integer(ik),            intent(in)      :: eparent_g
        integer(ik),            intent(in)      :: eparent_l

        integer(ik) :: ierr, vsize
        logical     :: reallocate

        !
        ! Set dense-vector integer data
        !
        self%dparent_g_ = dparent_g
        self%dparent_l_ = dparent_l
        self%eparent_g_ = eparent_g
        self%eparent_l_ = eparent_l

        self%nterms_ = nterms
        self%nvars_  = nvars
        self%ntime_  = ntime


        !
        ! Compute total number of elements for densevector storage
        !
        vsize = ntime * nterms * nvars


        !
        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        !
        if (allocated(self%vec)) then

            reallocate = (vsize /= size(self%vec))
            if (reallocate) then
                deallocate(self%vec)
                allocate(self%vec(vsize), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        else
            allocate(self%vec(vsize), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if


        !
        ! Initialize to zero
        !
        self%vec = ZERO

    end subroutine init
    !******************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_time_start(self,itime) result(istart)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: itime

        integer(ik) :: istart

        istart = (self%nvars_ * self%nterms_)*(itime - 1) + 1

    end function get_time_start
    !*****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_time_end(self,itime) result(iend)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: itime

        integer(ik) :: iend

        iend = self%nvars_ * self%nterms_ * itime

    end function get_time_end
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !----------------------------------------------------------------------------------------
    subroutine settime(self,itime,vals)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: itime
        real(rk),               intent(in)      :: vals(:)

        integer(ik) :: istart, iend

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (self%nvars_ * self%nterms_)*(itime - 1) + 1
        iend = istart + (self%nvars_ * self%nterms_ - 1)

        !
        ! Set modes
        !
        self%vec(istart:iend) = vals

    end subroutine settime
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !----------------------------------------------------------------------------------------
    function gettime(self,itime) result(vals_out)
        class(densevector_t),   intent(in)      :: self
        integer(ik),            intent(in)      :: itime

        integer(ik)             :: istart, iend
        real(rk),   allocatable :: vals_out(:)

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (self%nvars_ * self%nterms_)*(itime - 1) + 1
        iend = istart + (self%nvars_ * self%nterms_ - 1)

        !
        ! Return values at particular time
        !
        vals_out = self%vec(istart:iend)

    end function gettime
    !****************************************************************************************









    !> Function returns the stored vector data associated with variable index ivar
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar        Integer index of the variable, for which modes will be returned
    !!  @result     modes_out   Array of modes from the variable, ivar
    !!
    !!
    !!  @author Matteo Ugolotti + Sharma Mayank
    !!  @date   11/14/2016
    !!
    !!  @param[in]  itime       Integer index of the time level, for which modes will be returned
    !!
    !-----------------------------------------------------------------------------------------
    function getvar(self,ivar,itime) result(modes_out)
        class(densevector_t),   intent(in)      :: self
        integer(ik),            intent(in)      :: ivar, itime

        real(rk)                                :: modes_out(self%nterms_)
        integer(ik)                             :: istart, iend
       
        
        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + (self%nvars_*self%nterms_)*(itime - 1) + 1
        iend = istart + (self%nterms_ - 1)

        !
        ! Return modes
        !
        modes_out = self%vec(istart:iend)

    end function getvar
    !****************************************************************************************













    !>  Set the modes for a particular variable at a particular time
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar    Integer index of the variable being set
    !!  @param[in]  vals    Array of mode values that will be set
    !!
    !!  @author Matteo Ugolotti + Sharma Mayank
    !!  @date   11/14/2016
    !!
    !!  @param[in]  itime   Integer index of the variable time level being set 
    !!
    !-----------------------------------------------------------------------------------------
    subroutine setvar(self,ivar,itime,vals)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: ivar, itime
        real(rk),               intent(in)      :: vals(:)

        integer(ik) :: istart, iend    


        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + (self%nvars_*self%nterms_)*(itime - 1) + 1
        iend = istart + (self%nterms_ - 1) 

        !
        ! Set modes
        !
        self%vec(istart:iend) = vals

    end subroutine setvar
    !****************************************************************************************







    !>  Return an individual mode from the polynomial expansion of a variable
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar    Integer index of the variable
    !!  @param[in]  iterm   Integer index of the mode in the expansion to be returned
    !!
    !!  @author Matteo Ugolotti + Sharma Mayank
    !!  @date   11/14/2016
    !!
    !!  @param[in]  itime   Integer index of the time level, for which modes will be returned
    !!
    !-----------------------------------------------------------------------------------------
    function getterm(self,ivar,iterm,itime) result(mode_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: ivar
        integer(ik),            intent(in)  :: iterm
        integer(ik),            intent(in)  :: itime

        real(rk)    :: mode_out
        integer(ik) :: istart, iterm_g
 
        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + (self%nvars_*self%nterms_)*(itime - 1) + 1
        iterm_g = istart + iterm
        
        !
        ! Get mode
        !
        mode_out = self%vec(iterm_g)

    end function getterm
    !****************************************************************************************








    !>  Set an individual mode in the polynomial expansion of a variable
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  ivar        Integer index of the variable
    !!  @param[in]  iterm       Integer index of the mode in the expansion to be set
    !!  @param[in]  mode_in     Floating point value, which is the mode amplitude to be set
    !!
    !!  @author Matteo Ugolotti + Sharma Mayank
    !!  @date   11/14/2016
    !!
    !!  @param[in]  itime       Integer index of time level, for which modes will be set
    !!
    !--------------------------------------------------------------------------------
    subroutine setterm(self,ivar,iterm,itime,mode_in)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: ivar, iterm, itime
        real(rk),               intent(in)      :: mode_in

        integer(ik) :: istart, iterm_g

        !
        ! Compute start and end indices for accessing modes of a variable
        !
        istart = (ivar-1) * self%nterms_ + (self%nvars_*self%nterms_)*(itime - 1) + 1
        iterm_g = istart + iterm


        !
        ! Set mode
        !
        self%vec(iterm_g) = mode_in

    end subroutine setterm
    !********************************************************************************











    !> Function that returns number of entries in block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------
    function nentries(self) result(n)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: n

        n = size(self%vec)

    end function nentries
    !********************************************************************************









    !> Function that returns nterms_ private component
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------
    pure function nterms(self) result(nterms_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik)                         :: nterms_out

        nterms_out = self%nterms_

    end function nterms
    !********************************************************************************








    !> Function that returns nvars_ private component
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------
    pure function nvars(self) result(nvars_out)
        class(densevector_t),   intent(in)  :: self
        integer(ik)                         :: nvars_out

        nvars_out = self%nvars_

    end function nvars
    !********************************************************************************








    !> Function that returns ntime_ private component
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/03/2016
    !!
    !--------------------------------------------------------------------------------
    pure function ntime(self) result(ntime_out)
        class(densevector_t),  intent(in)  :: self

        integer(ik)                        :: ntime_out
        
        ntime_out = self%ntime_
        
    end function ntime
    !********************************************************************************
  
    






    !>  Subroutine to change ntime in the densevector to a value passed in 
    !!  Called in type_chidg_vector
    !!
    !!  @author Mayank Sharma
    !!  @date   3/9/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine set_ntime(self,ntime_in)
        class(densevector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: ntime_in

        !
        ! Set ntime
        !
        self%ntime_ = ntime_in


        !
        ! Reinitialize the densevector
        ! TODO: Move outside this subroutine?
        !
        call init(self,self%nterms_,self%nvars_,self%ntime_,self%dparent_g_, &
                  self%dparent_l_, self%eparent_g_, self%eparent_l_)


    end subroutine set_ntime
    !********************************************************************************
  
    







    !> Function that returns index of block parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/30/2016
    !!
    !--------------------------------------------------------------------------------
    function dparent_g(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%dparent_g_

    end function dparent_g
    !********************************************************************************



    !> Function that returns index of block parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/30/2016
    !!
    !--------------------------------------------------------------------------------
    function dparent_l(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%dparent_l_

    end function dparent_l
    !********************************************************************************



    !> Function that returns index of block parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/30/2016
    !!
    !--------------------------------------------------------------------------------
    function eparent_g(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%eparent_g_

    end function eparent_g
    !********************************************************************************


    !> Function that returns index of block parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/30/2016
    !!
    !--------------------------------------------------------------------------------
    function eparent_l(self) result(par)
        class(densevector_t),   intent(in)      :: self
        integer(ik)                             :: par

        par = self%eparent_l_

    end function eparent_l
    !********************************************************************************










    !> Zero vector storage, self%vec
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine clear(self)
        class(densevector_t),   intent(inout)   :: self

        self%vec = ZERO

    end subroutine clear
    !********************************************************************************





    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !!  @param[in]  nterms_r    Number of terms in the restricted basis.
    !!
    !--------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        integer(ik)             :: nvars, ntime, dparent_g, dparent_l, eparent_g, eparent_l, &
                                   ivar, itime
        real(rk),   allocatable :: field(:)
        type(densevector_t)     :: restricted


        ! Get parameters from initial vector
        nvars     = self%nvars()
        ntime     = self%ntime()
        dparent_g = self%dparent_g()
        dparent_l = self%dparent_l()
        eparent_g = self%eparent_g()
        eparent_l = self%eparent_l()


        ! Size to restricted basis
        call restricted%init(nterms_r,nvars,ntime,dparent_g,dparent_l,eparent_g,eparent_l)


        ! Copy restricted fields
        do ivar = 1,nvars
            do itime = 1,ntime
                field   = self%getvar(ivar,itime)
                call restricted%setvar(ivar,itime,field(1:nterms_r) )
            end do
        end do

    end function restrict
    !********************************************************************************







    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !!  @param[in]  nterms_p    Number of terms in the restricted basis.
    !!
    !--------------------------------------------------------------------------------
    function prolong(self,nterms_p) result(prolonged)
        class(densevector_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_p

        integer(ik)             :: nterms, nvars, ntime, dparent_g, dparent_l,  &
                                   eparent_g, eparent_l, ivar, itime, ierr
        real(rk),   allocatable :: field(:)
        type(densevector_t)     :: prolonged


        ! Get parameters from initial vector
        nterms    = self%nterms()
        nvars     = self%nvars()
        ntime     = self%ntime()
        dparent_g = self%dparent_g()
        dparent_l = self%dparent_l()
        eparent_g = self%eparent_g()
        eparent_l = self%eparent_l()


        ! Size to restricted basis
        call prolonged%init(nterms_p,nvars,ntime,dparent_g,dparent_l,eparent_g,eparent_l)


        ! Copy restricted fields
        allocate(field(nterms_p), stat=ierr)
        if (ierr /= 0) call AllocationError
        field = ZERO

        do ivar = 1,nvars
            do itime = 1,ntime
                field(1:nterms) = self%getvar(ivar,itime)
                call prolonged%setvar(ivar,itime,field)
            end do
        end do

    end function prolong
    !********************************************************************************





























    !>  Return the sum of squared terms in the dense vector.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2017
    !!
    !!
    !-------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(densevector_t),   intent(in)  :: self

        real(rk)    :: res

        ! Square vector values and sum
        res = sum( self%vec**TWO )

    end function sumsqr
    !*************************************************************************************






    !>  Return the sum of squared terms in the dense vector for each field independently.
    !!
    !!  Returns
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    function sumsqr_fields(self) result(res)
        class(densevector_t),   intent(in)  :: self

        real(rk)    :: res(self%nvars_)
        real(rk)    :: modes(self%nterms_)
        integer(ik) :: ifield, itime


        do ifield = 1,self%nvars_
            res(ifield) = 0

            do itime = 1,self%ntime_
                modes = self%getvar(ifield,itime)
                res(ifield) = res(ifield) + sum( modes**TWO )
            end do

        end do

    end function sumsqr_fields
    !************************************************************************************






    !-------------------------------------------------------------------
    !-------------      OPERATOR IMPLEMENTATIONS    --------------------
    !-------------------------------------------------------------------
    !
    ! @author Mayank Sharma + Matteo Ugolotti
    ! @date   11/15/2016
    ! 
    ! Added ntime_ attribute to all operators
    !
    elemental function mult_real_dv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res

        res%dparent_g_ = right%dparent_g_
        res%dparent_l_ = right%dparent_l_
        res%eparent_g_ = right%eparent_g_
        res%eparent_l_ = right%eparent_l_
        res%nvars_     = right%nvars_
        res%nterms_    = right%nterms_
        res%ntime_     = right%ntime_

        res%vec        = left * right%vec

    end function

   
    elemental function mult_dv_real(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right
        type(densevector_t) :: res


        res%dparent_g_ = left%dparent_g_
        res%dparent_l_ = left%dparent_l_
        res%eparent_g_ = left%eparent_g_
        res%eparent_l_ = left%eparent_l_
        res%nvars_     = left%nvars_
        res%nterms_    = left%nterms_
        res%ntime_     = left%ntime_

        res%vec        = left%vec * right

    end function








    elemental function div_real_dv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%dparent_g_ = right%dparent_g_
        res%dparent_l_ = right%dparent_l_
        res%eparent_g_ = right%eparent_g_
        res%eparent_l_ = right%eparent_l_
        res%nvars_     = right%nvars_
        res%nterms_    = right%nterms_
        res%ntime_     = right%ntime_

        res%vec        = left / right%vec

    end function


    elemental function div_dv_real(left,right) result(res)
        type(densevector_t),        intent(in)  :: left
        real(rk),                   intent(in)  :: right
        type(densevector_t) :: res


        res%dparent_g_ = left%dparent_g_
        res%dparent_l_ = left%dparent_l_
        res%eparent_g_ = left%eparent_g_
        res%eparent_l_ = left%eparent_l_
        res%nvars_     = left%nvars_
        res%nterms_    = left%nterms_
        res%ntime_     = left%ntime_

        res%vec        = left%vec / right

    end function







    elemental function add_dv_dv(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%dparent_g_ = left%dparent_g_
        res%dparent_l_ = left%dparent_l_
        res%eparent_g_ = left%eparent_g_
        res%eparent_l_ = left%eparent_l_
        res%nvars_     = left%nvars_
        res%nterms_    = left%nterms_
        res%ntime_     = left%ntime_

        res%vec        = left%vec + right%vec

    end function




    elemental function sub_dv_dv(left,right) result(res)
        type(densevector_t),    intent(in)  :: left
        type(densevector_t),    intent(in)  :: right
        type(densevector_t) :: res


        res%dparent_g_ = left%dparent_g_
        res%dparent_l_ = left%dparent_l_
        res%eparent_g_ = left%eparent_g_
        res%eparent_l_ = left%eparent_l_
        res%nvars_     = left%nvars_
        res%nterms_    = left%nterms_
        res%ntime_     = left%ntime_

        res%vec        = left%vec - right%vec

    end function


    












end module type_densevector
