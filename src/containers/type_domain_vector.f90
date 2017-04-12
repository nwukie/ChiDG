module type_domain_vector
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, TWO
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_Recv, MPI_ANY_TAG, MPI_STATUS_IGNORE, MPI_INTEGER4
    use type_domain,        only: domain_t
    use type_ivector,       only: ivector_t
    use type_densevector
    use DNAD_D
    implicit none



    !> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------
    type, public :: domain_vector_t

        type(densevector_t), allocatable :: vecs(:)     !< Local element vectors

    contains

        generic,    public  :: init => init_local, init_recv    !< Initialize vector
        procedure, private  :: init_local                       !< Initialize vector to store data for local domain
        procedure, private  :: init_recv                        !< Initialize vector to store data for domains being received from another processor

        procedure,  public  :: distribute               !< Given a full-vector representation, distribute it to the denseblock format
        procedure,  public  :: clear                    !< Zero all vector storage elements
        
        procedure,  public  :: norm             !< Return the L2 vector norm of the block-vector
        procedure,  public  :: sumsqr           !< Return the sum of the squared block-vector entries
        procedure,  public  :: sumsqr_fields    !< Return the sum of squared entries for fields independently
        procedure,  public  :: nentries
        procedure,  public  :: dump

        final :: destructor

    end type domain_vector_t
    !*********************************************************************************************











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



    !>  Subroutine for initializing blockvector storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    mesh_t instance containing initialized elements and faces
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_local(self,domain)
        class(domain_vector_t),   intent(inout) :: self
        type(domain_t),         intent(in)    :: domain

        integer(ik) :: nelem, ierr, ielem, nterms, neqns, ntime
        integer(ik) :: dparent_g, dparent_l, eparent_g, eparent_l
        logical     :: new_elements


        nelem = domain%nelem  ! Number of elements in the local block

        !
        ! ALLOCATE SIZE FOR 'vecs'
        ! If vector was already allocated, deallocate and then reallocate vector size
        ! Reallocation would take place if the number of elements were changed
        !
        if (allocated(self%vecs)) then
            !
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            !
            new_elements = (domain%nelem /= size(self%vecs))
            if (new_elements) then
                deallocate(self%vecs)
                allocate(self%vecs(nelem), stat=ierr)
            end if

        else

            allocate(self%vecs(nelem), stat=ierr)

        end if
        if (ierr /= 0) call AllocationError




        !
        ! Loop through elements and call initialization for densevectors
        !
        do ielem = 1,domain%nelem
            dparent_g = domain%elems(ielem)%idomain_g
            dparent_l = domain%elems(ielem)%idomain_l
            eparent_g = domain%elems(ielem)%ielement_g
            eparent_l = domain%elems(ielem)%ielement_l
            nterms    = domain%elems(ielem)%nterms_s
            neqns     = domain%elems(ielem)%neqns
            ntime     = domain%elems(ielem)%ntime

            ! Call densevector initialization routine
            call self%vecs(ielem)%init(nterms,neqns,ntime,dparent_g,dparent_l,eparent_g,eparent_l)

        end do




    end subroutine init_local
    !*****************************************************************************************













    !>  Subroutine for initializing blockvector storage for elements being received from 
    !!  another processor.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!  @param[in]  domain    domain_t instance containing initialized elements and faces
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_recv(self,iproc)
        class(domain_vector_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: iproc

        type(ivector_t) :: recv_elems
        integer(ik)     :: nelem_recv, ielem_recv, ierr, ielem, iface, nterms, neqns, loc, recv_element
        integer(ik)     :: idomain_g, idomain_l, ielement_g, ielement_l, ntime
        logical         :: new_elements, proc_element, already_added, comm_element



        !
        ! Get the domain index we are receiving
        !
        call MPI_Recv(idomain_g, 1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
        call MPI_Recv(idomain_l, 1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)


        !
        ! Get the number of elements being received from domain
        !
        call MPI_Recv(nelem_recv, 1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
        


        

        !
        ! Allocate densevector size. If vector was already allocated, deallocate and then reallocate vector size
        ! Reallocation would take place if the number of elements were changed
        if (allocated(self%vecs)) then
            !
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            !
            new_elements = (nelem_recv /= size(self%vecs))
            if (new_elements) then
                deallocate(self%vecs)
                allocate(self%vecs(nelem_recv), stat=ierr)
            end if

        else

            allocate(self%vecs(nelem_recv), stat=ierr)

        end if
        if (ierr /= 0) call AllocationError




        !
        ! Loop through and recv element information from sending proc and call initialization on densevector storage
        !
        do ielem_recv = 1,nelem_recv


            call MPI_Recv(ielement_g, 1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(ielement_l, 1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(nterms,     1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(neqns,      1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(ntime,      1, MPI_INTEGER4, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

            !
            ! Call densevector initialization routine
            !
            call self%vecs(ielem_recv)%init(nterms,neqns,ntime,idomain_g,idomain_l,ielement_g,ielement_l)


        end do




    end subroutine init_recv
    !************************************************************************************












    !> Given a full-vector, distribute it's values to the blockvector format, 
    !! assuming contiguously stored data for each element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  fullvec     Full-vector
    !!
    !------------------------------------------------------------------------------------
    subroutine distribute(self,fullvec)
        class(domain_vector_t),    intent(inout)   :: self
        real(rk),                intent(in)      :: fullvec(:) 

        integer(ik)     :: ndof, ielem, nvars, nterms, fstart, fend, ndof_l

        !
        ! Compute total entries allocated in the blockvector container. 
        !
        ndof = 0
        do ielem = 1,size(self%vecs)
            nvars = self%vecs(ielem)%nvars()
            nterms = self%vecs(ielem)%nterms()
            
            ndof = ndof  +  (nvars * nterms)
        end do

        ! Test that the number of dof's match between the full and block format's
        if (ndof /= size(fullvec) ) call chidg_signal(FATAL,"domain_vector_t%distribute: Storage sizes of full-vector and block-vector are not equal.")



        !
        ! Loop through elements and store data from full-vector
        !
        fstart = 1
        do ielem = 1,size(self%vecs)
            ! Get number of entries for current block
            nvars = self%vecs(ielem)%nvars()
            nterms = self%vecs(ielem)%nterms()
            ndof_l = nvars * nterms

            fend = fstart + (ndof_l-1)


            !& DEBUG
            if (fend > size(fullvec) ) call chidg_signal(FATAL,"blockvector%distribute: array bounds exceeded")


            self%vecs(ielem)%vec = fullvec(fstart:fend)

            fstart = fend + 1
        end do




    end subroutine distribute
    !***********************************************************************************











    !> Zero all vector storage elements in blockvector%vecs
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine clear(self)
        class(domain_vector_t),   intent(inout)   :: self

        integer(ik) :: iblk

        ! Call clear for each densevector component
        do iblk = 1,size(self%vecs)
            call self%vecs(iblk)%clear()
        end do

    end subroutine clear
    !****************************************************************************************












    !> Compute the L2-norm of the block vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function norm(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        real(rk)    :: res
        integer(ik) :: ielem

!        res = ZERO
!
!        ! Loop through block vectors and compute contribution to vector L2-Norm
!        do ielem = 1,size(self%vecs)
!            ! Square vector values and sum
!            res = res + sum( self%vecs(ielem)%vec ** TWO )
!        end do
        res = self%sumsqr()


        ! Take the square root of the result
        res = sqrt(res)

    end function norm
    !*****************************************************************************************









    !>  Sum of the squared block-vector entries
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        real(rk)    :: res
        integer(ik) :: ielem

        res = ZERO

        ! Loop through block vectors and compute contribution to sum of squared entries
        do ielem = 1,size(self%vecs)
            ! Square vector values and sum
            res = res + sum( self%vecs(ielem)%vec ** TWO )
        end do


    end function sumsqr
    !*****************************************************************************************










    !>  Sum of the squared block-vector entries, separated by field.
    !!
    !!  Returns an array of values. The sum of the squared values from each field independently.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function sumsqr_fields(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        real(rk),   allocatable :: res(:)
        integer(ik)             :: ielem

        ! Size number of field residuals being computed
        res = self%vecs(1)%sumsqr_fields()
        res = ZERO

        ! Loop through densevectors and compute contribution to sum of squared entries
        do ielem = 1,size(self%vecs)
            res = res + self%vecs(ielem)%sumsqr_fields()
        end do


    end function sumsqr_fields
    !*****************************************************************************************









    !> Compute number of entries in the block vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function nentries(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        integer(ik) :: res, ielem

        res = ZERO
        !
        ! Loop through block vectors and compute contribution to number of entries
        !
        do ielem = 1,size(self%vecs)

            res = res + self%vecs(ielem)%nentries()

        end do

    end function nentries
    !*****************************************************************************************







    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------
    subroutine dump(self)
        class(domain_vector_t),   intent(in)  :: self
        integer(ik) :: ielem, ientry

        do ielem = 1,size(self%vecs)
            print*, ielem
            do ientry = 1,size(self%vecs(ielem)%vec)
                print*, self%vecs(ielem)%vec(ientry)
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
        type(domain_vector_t),    intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left * right%vecs

    end function mult_real_bv
    !************************************************************************




    !> Multiply blockvector * real
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function mult_bv_real(left,right) result(res)
        type(domain_vector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left%vecs * right

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
        type(domain_vector_t),    intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left / right%vecs

    end function div_real_bv
    !************************************************************************




    !> Divide blockvector / real
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function div_bv_real(left,right) result(res)
        type(domain_vector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left%vecs / right

    end function div_bv_real
    !*************************************************************************






    !> Add blockvector + blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function add_bv_bv(left,right) result(res)
        type(domain_vector_t),  intent(in)  :: left
        type(domain_vector_t),  intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left%vecs + right%vecs

    end function add_bv_bv
    !*************************************************************************





    !> Add blockvector - blockvector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function sub_bv_bv(left,right) result(res)
        type(domain_vector_t),  intent(in)  :: left
        type(domain_vector_t),  intent(in)  :: right

        type(domain_vector_t) :: res

        res%vecs = left%vecs - right%vecs

    end function sub_bv_bv
    !*************************************************************************





    subroutine destructor(self)
        type(domain_vector_t), intent(inout) :: self

    end subroutine

end module type_domain_vector
