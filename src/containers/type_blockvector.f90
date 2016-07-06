module type_blockvector
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, TWO
    use type_mesh,          only: mesh_t
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
    !------------------------------------------------------------------------------------------------------------------------
    type, public :: blockvector_t

        type(densevector_t), allocatable :: vecs(:)     !< Local element vectors

    contains

        generic,    public  :: init => init_local, init_recv    !< Initialize vector
        procedure, private  :: init_local                       !< Initialize vector to store data for local domain
        procedure, private  :: init_recv                        !< Initialize vector to store data for domains being received from another processor

        procedure,  public  :: distribute               !< Given a full-vector representation, distribute it to the denseblock format
        procedure,  public  :: clear                    !< Zero all vector storage elements
        
        procedure,  public  :: norm                     !< Return the L2 vector norm of the block-vector
        procedure,  public  :: sumsqr                   !< Return the sum of the squared block-vector entries
        procedure,  public  :: nentries
        procedure,  public  :: dump

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



    !>  Subroutine for initializing blockvector storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    mesh_t instance containing initialized elements and faces
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_local(self,mesh)
        class(blockvector_t),   intent(inout) :: self
        type(mesh_t),           intent(in)    :: mesh

        integer(ik) :: nelem, ierr, ielem, nterms, neqns
        integer(ik) :: dparent_g, dparent_l, eparent_g, eparent_l
        logical     :: new_elements


        nelem = mesh%nelem  ! Number of elements in the local block

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
            new_elements = (mesh%nelem /= size(self%vecs))
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
        do ielem = 1,mesh%nelem
            dparent_g = mesh%elems(ielem)%idomain_g
            dparent_l = mesh%elems(ielem)%idomain_l
            eparent_g = mesh%elems(ielem)%ielement_g
            eparent_l = mesh%elems(ielem)%ielement_l
            nterms = mesh%elems(ielem)%nterms_s
            neqns  = mesh%elems(ielem)%neqns

            ! Call densevector initialization routine
            call self%vecs(ielem)%init(nterms,neqns,dparent_g,dparent_l,eparent_g,eparent_l)

        end do




    end subroutine init_local
    !*******************************************************************************************************************













    !>  Subroutine for initializing blockvector storage for elements being received from another processor
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!  @param[in]  mesh    mesh_t instance containing initialized elements and faces
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_recv(self,mesh,iproc,icomm,idom_recv)
        class(blockvector_t),   intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: iproc
        integer(ik),            intent(in)      :: icomm
        integer(ik),            intent(in)      :: idom_recv

        type(ivector_t) :: recv_elems
        integer(ik)     :: nelem_recv, ielem_recv, ierr, ielem, iface, nterms, neqns, loc, recv_element
        integer(ik)     :: dparent_g, dparent_l, eparent_g, eparent_l
        logical         :: new_elements, proc_element, already_added, comm_element


        !
        ! Compute number of densevectors being received from proc. One for each comm neighbor element.
        !
        do ielem = 1,mesh%nelem
            do iface = 1,size(mesh%faces,2)

                proc_element = ( iproc == mesh%faces(ielem,iface)%ineighbor_proc )
            

                ! Test if element has been added to recv list
                if (proc_element) then
                    ! Get index of recv element and test if it was already added
                    recv_element = mesh%faces(ielem,iface)%ineighbor_element_g

                    loc = recv_elems%loc(recv_element)
                    already_added = (loc /= 0)

                    if (.not. already_added) then
                        call recv_elems%push_back(recv_element)
                    end if

                end if

            end do ! iface
        end do ! ielem


        nelem_recv = recv_elems%size()
        if (nelem_recv == 0) call chidg_signal(FATAL,"blockvector%init_recv: no elements detected in mesh that are being received")



        !
        ! Allocate densevector size. If vector was already allocated, deallocate and then reallocate vector size
        ! Reallocation would take place if the number of elements were changed
        !
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
        ! Loop through elements and call initialization for densevectors being received. Loop through all elements so
        ! recv_comm, recv_domain, and recv_element can be initialized since multiple local elements could be receiving 
        ! the same element and need the recv indices initialized
        !
        ! ASSUME: nterms_s of neighbor being received is same as ielem
        !
        do ielem_recv = 1,recv_elems%size()

            recv_element = recv_elems%at(ielem_recv)

            do ielem = 1,mesh%nelem
                do iface = 1,size(mesh%faces,2)

                    proc_element = ( iproc == mesh%faces(ielem,iface)%ineighbor_proc )
                
                    comm_element = ( recv_element == mesh%faces(ielem,iface)%ineighbor_element_g )
                    

                    if (proc_element .and. comm_element) then
                        dparent_g = mesh%faces(ielem,iface)%ineighbor_domain_g
                        dparent_l = mesh%faces(ielem,iface)%ineighbor_domain_l
                        eparent_g = mesh%faces(ielem,iface)%ineighbor_element_g
                        eparent_l = mesh%faces(ielem,iface)%ineighbor_element_l

                        nterms = mesh%elems(ielem)%nterms_s ! assume same throughout domain
                        neqns  = mesh%elems(ielem)%neqns    ! assume same throughout domain

                        !
                        ! Call densevector initialization routine
                        !
                        call self%vecs(ielem_recv)%init(nterms,neqns,dparent_g,dparent_l,eparent_g,eparent_l)


                        mesh%faces(ielem,iface)%recv_comm    = icomm
                        mesh%faces(ielem,iface)%recv_domain  = idom_recv
                        mesh%faces(ielem,iface)%recv_element = ielem_recv

                    end if


                end do !iface
            end do ! ielem


        end do !ielem_recv



    end subroutine init_recv
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
        do ielem = 1,size(self%vecs)
            nvars = self%vecs(ielem)%nvars()
            nterms = self%vecs(ielem)%nterms()
            
            ndof = ndof  +  (nvars * nterms)
        end do

        ! Test that the number of dof's match between the full and block format's
        if (ndof /= size(fullvec) ) call chidg_signal(FATAL,"blockvector_t%distribute: Storage sizes of full-vector and block-vector are not equal.")



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
    !***************************************************************************************************************************************














    !> Zero all vector storage elements in blockvector%vecs
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(blockvector_t),   intent(inout)   :: self

        integer(ik) :: iblk

        ! Call clear for each densevector component
        do iblk = 1,size(self%vecs)
            call self%vecs(iblk)%clear()
        end do

    end subroutine clear
    !****************************************************************************************************************************************














    !> Compute the L2-norm of the block vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------------
    function norm(self) result(res)
        class(blockvector_t),   intent(in)  :: self

        real(rk)    :: res
        integer(ik) :: ielem

        res = ZERO

        ! Loop through block vectors and compute contribution to vector L2-Norm
        do ielem = 1,size(self%vecs)
            ! Square vector values and sum
            res = res + sum( self%vecs(ielem)%vec ** TWO )
        end do


        ! Take the square root of the result
        res = sqrt(res)

    end function norm
    !*****************************************************************************************************************************************









    !>  Sum of the squared block-vector entries
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(blockvector_t),   intent(in)  :: self

        real(rk)    :: res
        integer(ik) :: ielem

        res = ZERO

        ! Loop through block vectors and compute contribution to sum of squared entries
        do ielem = 1,size(self%vecs)
            ! Square vector values and sum
            res = res + sum( self%vecs(ielem)%vec ** TWO )
        end do


    end function sumsqr
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
        do ielem = 1,size(self%vecs)

            res = res + self%vecs(ielem)%nentries()

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
        type(blockvector_t),    intent(in)  :: right

        type(blockvector_t) :: res

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
        type(blockvector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(blockvector_t) :: res

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
        type(blockvector_t),    intent(in)  :: right

        type(blockvector_t) :: res

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
        type(blockvector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(blockvector_t) :: res

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
        type(blockvector_t),  intent(in)  :: left
        type(blockvector_t),  intent(in)  :: right

        type(blockvector_t) :: res

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
        type(blockvector_t),  intent(in)  :: left
        type(blockvector_t),  intent(in)  :: right

        type(blockvector_t) :: res

        res%vecs = left%vecs - right%vecs

    end function sub_bv_bv
    !*************************************************************************





    subroutine destructor(self)
        type(blockvector_t), intent(inout) :: self

    end subroutine

end module type_blockvector
