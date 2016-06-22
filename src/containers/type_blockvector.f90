module type_blockvector
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, TWO
    use type_mesh,          only: mesh_t
    use type_densevector
    use DNAD_D
    implicit none


    !> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------------
    type, public :: blockvector_t

        type(densevector_t), allocatable :: lvecs(:)        !< Local element vectors
        integer(ik),         allocatable :: ldata(:,:)      !< Local block data     (nvars, nterms)


    contains
        ! Initializers
        generic,   public   :: init => init_vector          !< Initialize local vector
        procedure, private  :: init_vector

        procedure, public   :: distribute                   !< Given a full-vector representation, distribute it to the denseblock format
        procedure, public   :: clear                        !< Zero all vector storage elements
        
        procedure, public   :: norm
        procedure, public   :: nentries
        procedure, public   :: dump

        final :: destructor

    end type blockvector_t
    !*************************************************************************************************************************











    !-----------------      OPERATORS       --------------------------


    public operator (*)
    interface operator (*)
        module procedure mult_real_bv   ! real * blockvector
        module procedure mult_bv_real   ! blockvector * real
    end interface



    public operator (/)
    interface operator (/)
        module procedure div_real_bv    ! real / blockvector
        module procedure div_bv_real    ! blockvector / real
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_bv_bv      ! blockvector - blockvector
    end interface


    public operator (+)
    interface operator (+)
        module procedure add_bv_bv     ! blockvector + blockvector
    end interface












    private
contains



    !> Subroutine for initializing blockvector storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    mesh_t instance containing initialized elements and faces
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_vector(self,mesh)
        class(blockvector_t), intent(inout), target  :: self
        class(mesh_t),        intent(in)             :: mesh

        integer(ik)                     :: nelem, nblk, ierr, ielem, iblk, size1d, parent, nterms, neqns
        logical                         :: new_elements
        type(densevector_t), pointer    :: temp(:)

        nelem = mesh%nelem  ! Number of elements in the local block


        !
        ! ALLOCATE SIZE FOR 'lvecs'
        ! If vector was already allocated, deallocate and then reallocate vector size
        ! Reallocation would take place if the number of elements were changed
        !
        if (allocated(self%lvecs)) then
            !
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            !
            new_elements = (mesh%nelem /= size(self%lvecs))
            if (new_elements) then
                deallocate(self%lvecs, self%ldata)
                allocate(self%lvecs(nelem), self%ldata(nelem,2), stat=ierr)
            end if

        else

            allocate(self%lvecs(nelem), self%ldata(nelem,2), stat=ierr)

        end if
        if (ierr /= 0) call AllocationError




        !
        ! Loop through elements and call initialization for densevectors
        !
        do ielem = 1,mesh%nelem
            parent = mesh%elems(ielem)%ielement_l
            nterms = mesh%elems(ielem)%nterms_s
            neqns  = mesh%elems(ielem)%neqns

            ! Call densevector initialization routine
            call self%lvecs(ielem)%init(nterms,neqns,parent)

            ! Store data about number of equations and number of terms in solution expansion
            self%ldata(ielem,1) = mesh%elems(ielem)%neqns
            self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
        end do




    end subroutine init_vector
    !*******************************************************************************************************************














    !> Given a full-vector, distribute it's values to the blockvector format, 
    !! assuming contiguously stored data for each element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  fullvec     Full-vector
    !!
    !--------------------------------------------------------------------------------------------------------------------
    subroutine distribute(self,fullvec)
        class(blockvector_t),    intent(inout)   :: self
        real(rk),                intent(in)      :: fullvec(:) 

        integer(ik)     :: ndof, ielem, nvars, nterms, fstart, fend, ndof_l

        !
        ! Compute total entries allocated in the blockvector container. 
        !
        ndof = 0
        do ielem = 1,size(self%lvecs)
            nvars  = self%ldata(ielem,1)
            nterms = self%ldata(ielem,2)
            
            ndof = ndof  +  (nvars * nterms)
        end do

        !
        ! Test that the number of dof's match between the full and block format's
        !
        if (ndof /= size(fullvec) ) call chidg_signal(FATAL,"blockvector_t%distribute: Storage sizes of full-vector and block-vector are not equal.")



        !
        ! Loop through elements and store data from full-vector
        !
        fstart = 1
        do ielem = 1,size(self%lvecs)
            ! Get number of entries for current block
            nvars  = self%ldata(ielem,1)
            nterms = self%ldata(ielem,2)
            ndof_l = nvars * nterms

            fend = fstart + (ndof_l-1)


            !& DEBUG
            if (fend > size(fullvec) ) call chidg_signal(FATAL,"blockvector%distribute: array bounds exceeded")


            self%lvecs(ielem)%vec = fullvec(fstart:fend)

            fstart = fend + 1
        end do




    end subroutine distribute
    !***************************************************************************************************************************************














    !> Zero all vector storage elements in blockvector%lvecs
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(blockvector_t),   intent(inout)   :: self

        integer(ik) :: iblk


        !
        ! Call clear for each densevector component
        !
        do iblk = 1,size(self%lvecs)
            call self%lvecs(iblk)%clear()
        end do


    end subroutine clear
    !****************************************************************************************************************************************














    !> Compute the L2-norm of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------------
    function norm(self) result(res)
        class(blockvector_t),   intent(inout)   :: self

        real(rk)    :: res
        integer(ik) :: ielem



        res = ZERO

        !
        ! Loop through block vectors and compute contribution to vector L2-Norm
        !
        do ielem = 1,size(self%lvecs)


            !
            ! Square vector values and sum
            !
            res = res + sum( self%lvecs(ielem)%vec ** TWO )

        end do


        !
        ! Take the square root of the result
        !
        res = sqrt(res)

    end function
    !*****************************************************************************************************************************************












    !> Compute number of entries in the block vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------
    function nentries(self) result(res)
        class(blockvector_t),   intent(in)  :: self

        integer(ik) :: res, ielem

        res = ZERO
        !
        ! Loop through block vectors and compute contribution to number of entries
        !
        do ielem = 1,size(self%lvecs)

            res = res + self%lvecs(ielem)%nentries()

        end do

    end function nentries
    !**********************************************************************************************************************************







    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------
    subroutine dump(self)
        class(blockvector_t),   intent(in)  :: self
        integer(ik) :: ielem, ientry

        do ielem = 1,size(self%lvecs)
            print*, ielem
            do ientry = 1,size(self%lvecs(ielem)%vec)
                print*, self%lvecs(ielem)%vec(ientry)
            end do
        end do


    end subroutine
    !************************************************









    !------------------------------------------------------------------------
    !
    !                       OPERATOR IMPLEMENTATIONS
    !
    !------------------------------------------------------------------------


    !> Multiply real * blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function mult_real_bv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(blockvector_t),    intent(in)  :: right

        type(blockvector_t), target     :: res



        res%ldata = right%ldata

        res%lvecs = left * right%lvecs

    end function mult_real_bv
    !************************************************************************




    !> Multiply blockvector * real
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function mult_bv_real(left,right) result(res)
        type(blockvector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(blockvector_t), target     :: res




        res%ldata = left%ldata

        res%lvecs = left%lvecs * right

    end function mult_bv_real
    !************************************************************************





    !> Divide real / blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function div_real_bv(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(blockvector_t),    intent(in)  :: right

        type(blockvector_t), target     :: res




        res%ldata = right%ldata

        res%lvecs = left / right%lvecs

    end function div_real_bv
    !************************************************************************




    !> Divide blockvector / real
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function div_bv_real(left,right) result(res)
        type(blockvector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(blockvector_t), target     :: res




        res%ldata = left%ldata

        res%lvecs = left%lvecs / right

    end function div_bv_real
    !*************************************************************************






    !> Add blockvector + blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function add_bv_bv(left,right) result(res)
        type(blockvector_t),  intent(in)  :: left
        type(blockvector_t),  intent(in)  :: right

        type(blockvector_t), target     :: res


        res%ldata = right%ldata

        res%lvecs = left%lvecs + right%lvecs

    end function add_bv_bv
    !*************************************************************************





    !> Add blockvector - blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function sub_bv_bv(left,right) result(res)
        type(blockvector_t),  intent(in)  :: left
        type(blockvector_t),  intent(in)  :: right

        type(blockvector_t), target     :: res



        res%ldata = right%ldata

        res%lvecs = left%lvecs - right%lvecs

    end function sub_bv_bv
    !*************************************************************************





    subroutine destructor(self)
        type(blockvector_t), intent(inout) :: self

    end subroutine

end module type_blockvector
