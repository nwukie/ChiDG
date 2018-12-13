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



    !>  Domain-level vector container.
    !!
    !!  Contains vector information for a single domain. For each element on the domain,
    !!  a densevector_t instance is created that contains data for the element.
    !!
    !!  Use in local problem:
    !!  ---------------------
    !!  For each domain in the processor-local data%mesh, an instance of domain_vector_t
    !!  is allocated on chidg_vector%dom(:) 
    !!
    !!  Use in global problem:
    !!  ----------------------
    !!  This is also used to receive domain data from other processors. The container
    !!  access chidg_vector%recv%comm(:)%dom(:) designates that for any processor
    !!  being communicated with, for each domain that is being received, an instance
    !!  of domain_vector_t is allocated to receive that data.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: domain_vector_t

        type(densevector_t), allocatable :: vecs(:)     ! Local element vectors

    contains

        generic,    public  :: init => init_local, init_recv    ! Initialize vector.
        procedure, private  :: init_local                       ! Initialize vector for local domain.
        procedure, private  :: init_recv                        ! Initialize vector for domain received from another proc.

        procedure,  public  :: distribute       ! Given a full-vector representation, distribute it to a dense array.
        procedure,  public  :: clear            ! Zero all vector storage elements.
        
        procedure,  public  :: norm             ! Return the L2 vector norm of the block-vector.
        procedure,  public  :: sumsqr           ! Return the sum of the squared block-vector entries.
        procedure,  public  :: sumsqr_fields    ! Return the sum of squared entries for fields independently.
        procedure,  public  :: nentries
        procedure,  public  :: nelements
        procedure,  public  :: dump

        procedure,  public  :: restrict
        procedure,  public  :: prolong

        final :: destructor

    end type domain_vector_t
    !**************************************************************************************











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
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%vecs(nelem), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if




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
    subroutine init_recv(self,proc)
        class(domain_vector_t), intent(inout)   :: self
        integer(ik),            intent(in)      :: proc

        type(ivector_t) :: recv_elems
        integer(ik)     :: nelem_recv, ielem_recv, ierr, ielem, iface, nterms,  &
                           neqns, loc, recv_element, idomain_g, idomain_l,      &
                           ielement_g, ielement_l, ntime, element_location(5),  &
                           element_data(8)
        logical         :: new_elements, proc_element, already_added, comm_element



!        !
!        ! Get the domain index we are receiving
!        !
!        call MPI_Recv(idomain_g, 1, MPI_INTEGER4, proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
!        call MPI_Recv(idomain_l, 1, MPI_INTEGER4, proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)


        !
        ! Get the number of elements being received from domain
        !
        call MPI_Recv(nelem_recv, 1, MPI_INTEGER4, proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
        


        

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
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%vecs(nelem_recv), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if




        !
        ! Loop through and recv element information from sending proc and call initialization on densevector storage
        !
        do ielem_recv = 1,nelem_recv


            call MPI_Recv(element_location, 5, MPI_INTEGER4, proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(element_data,     8, MPI_INTEGER4, proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            
            ! Interpret information that has been received
            idomain_g  = element_location(1)
            idomain_l  = element_location(2)
            ielement_g = element_location(3)
            ielement_l = element_location(4)
            nterms = element_data(5)
            neqns  = element_data(4)
            ntime  = element_data(7)

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
    !! !$OMP PARALLEL num_threads(2)
    !! !$OMP  PARALLEL DO 
    !! !$OMP& DEFAULT(SHARED) PRIVATE(ielem) 
    !! !$OMP& SCHEDULE(STATIC,CHUNK) 
    !! !$OMP& REDUCTION(+:res)
    !!
    !----------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        real(rk)    :: res
        integer(ik) :: ielem

        res = ZERO

! !$OMP PARALLEL num_threads(2)
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ielem) SCHEDULE(STATIC,CHUNK) REDUCTION(+:res)
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ielem) REDUCTION(+:res)

        ! Loop through block vectors and compute contribution to sum of squared entries
        do ielem = 1,size(self%vecs)
            ! Square vector values and sum
            res = res + sum( self%vecs(ielem)%vec ** TWO )
        end do

! !$OMP END PARALLEL DO

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









    !>  Compute number of entries in the block vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------
    function nentries(self) result(res)
        class(domain_vector_t),   intent(in)  :: self

        integer(ik) :: res, ielem

        ! Loop through block vectors and compute contribution to number of entries
        res = ZERO
        do ielem = 1,size(self%vecs)
            res = res + self%vecs(ielem)%nentries()
        end do

    end function nentries
    !******************************************************************************





    !>  Return number of densevector_t instances, one for each finite-element.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------
    function nelements(self) result(nelem)
        class(domain_vector_t), intent(in)  :: self

        integer(ik) :: nelem

        if (allocated(self%vecs)) then
            nelem = size(self%vecs)
        else
            nelem = 0
        end if

    end function nelements
    !******************************************************************************










    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------
    subroutine dump(self)
        class(domain_vector_t),   intent(in)  :: self
        integer(ik) :: ielem, ientry

        do ielem = 1,size(self%vecs)
            print*, ielem
            do ientry = 1,size(self%vecs(ielem)%vec)
                print*, self%vecs(ielem)%vec(ientry)
            end do
        end do


    end subroutine dump
    !******************************************************************************







    !>  Return a domain_vector_t restricted to the basis defined by nterms_r.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(domain_vector_t), intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(domain_vector_t)   :: restricted
        integer(ik)             :: ierr, ielem


        ! Alocate storage on restricted domain_vector
        allocate(restricted%vecs(self%nelements()), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Copy restricted densevector for each element
        do ielem = 1,self%nelements()
            restricted%vecs(ielem) = self%vecs(ielem)%restrict(nterms_r)
        end do !ielem


    end function restrict
    !******************************************************************************






    !>  Return a domain_vector_t prolonged to the basis defined by nterms_p.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2017
    !!
    !-----------------------------------------------------------------------------
    function prolong(self,nterms_p) result(prolonged)
        class(domain_vector_t), intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_p

        type(domain_vector_t)   :: prolonged
        integer(ik)             :: ierr, ielem


        ! Alocate storage on restricted domain_vector
        allocate(prolonged%vecs(self%nelements()), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Copy restricted densevector for each element
        do ielem = 1,self%nelements()
            prolonged%vecs(ielem) = self%vecs(ielem)%prolong(nterms_p)
        end do !ielem


    end function prolong
    !******************************************************************************






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
