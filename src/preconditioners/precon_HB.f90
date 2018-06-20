module precon_HB
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, ONE, ZERO
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK
    use mod_io,                 only: verbosity

    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_chidg_vector,      only: chidg_vector_t
    use type_element,           only: element_t
    use type_mesh,              only: mesh_t
    implicit none


    !----------------------------------------------------------------------------------
    type, public :: HB_dense_matrix_t
        real(rk), allocatable :: mat(:,:)
    contains
        procedure :: init => init_hb_dense_matrix
    end type HB_dense_matrix_t

    type, public :: HB_matrix_t
        type(HB_dense_matrix_t), allocatable :: elems(:)
    contains
        procedure :: init => init_hb_matrix
        procedure :: clear => clear_hb_matrix
    end type HB_matrix_t
    !**********************************************************************************


    !>  HB preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !----------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_HB_t
        type(HB_matrix_t)       :: HB_matrix
        real(rk),   allocatable :: D(:,:)
    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply
    end type precon_HB_t
    !**********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine init_hb_dense_matrix(self,element)
        class(HB_dense_matrix_t), intent(inout)   :: self
        type(element_t),    intent(in)      :: element

        integer(ik) :: ierr, mat_size

        mat_size = element%nterms_s*element%neqns*element%ntime 

        allocate(self%mat(mat_size,mat_size), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init_hb_dense_matrix
    !**********************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine init_hb_matrix(self,mesh)
        class(HB_matrix_t), intent(inout)   :: self
        type(mesh_t),       intent(in)      :: mesh

        integer(ik) :: ierr, ielem

        allocate(self%elems(mesh%domain(1)%nelements()), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ielem = 1,mesh%domain(1)%nelements()
            call self%elems(ielem)%init(mesh%domain(1)%elems(ielem))
        end do

    end subroutine init_hb_matrix
    !**********************************************************************************

    
    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine clear_hb_matrix(self)
        class(HB_matrix_t), intent(inout)   :: self

        integer(ik) :: ielem

        do ielem = 1,size(self%elems)
            self%elems(ielem)%mat = ZERO
        end do

    end subroutine clear_hb_matrix
    !**********************************************************************************




    !>  Initialize the HB preconditioner. This is for allocating storage. In this case, 
    !!  we allocate a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !!  @param[inout]   domain  domain_t instance containing a mesh component used to 
    !!                          initialize the block matrix
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_HB_t),     intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        call self%HB_matrix%init(data%mesh)
        call self%HB_matrix%clear()
        self%D = data%time_manager%D
        self%initialized = .true.

    end subroutine init
    !*******************************************************************************************






    !>  Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !-------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_HB_t),     intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: b

        integer(ik) :: ielem, itime, sz, irow_start, irow_end, icol_start, icol_end, &
                       itime_a, itime_b, idiag, ierr, ifield, nfields, nterms
        real(rk), allocatable   :: hb_mat(:,:)

        do ielem = 1,size(self%hb_matrix%elems)
            nterms  = A%dom(1)%lblks(ielem,1)%data_(1)%nterms
            nfields = A%dom(1)%lblks(ielem,1)%data_(1)%nfields

            ! Assemble diagonal elements
            do itime = 1,size(A%dom(1)%lblks,2)

                ! Access indices and size
                sz = size(A%dom(1)%lblks(ielem,itime)%data_(1)%mat,1)
                irow_start = 1 + (itime-1)*sz
                irow_end = irow_start + (sz-1)
                icol_start = 1 + (itime-1)*sz
                icol_end = icol_start + (sz-1)

                ! Get diagonal index and store
                idiag = A%dom(1)%lblks(ielem,itime)%get_diagonal()
                self%hb_matrix%elems(ielem)%mat(irow_start:irow_end,icol_start:icol_end) = A%dom(1)%lblks(ielem,itime)%data_(idiag)%mat

            end do !itime

            ! Assemble off-diagonal coupling 
            do itime_a = 1,size(A%dom(1)%lblks,2)
                do itime_b = 1,size(A%dom(1)%lblks,2)

                    ! LHS HB contribution for each variable
                    if (allocated(hb_mat)) deallocate(hb_mat)
                    allocate(hb_mat(nterms*nfields,nterms*nfields), stat=ierr)
                    if (ierr /= 0) call AllocationError
                    hb_mat(:,:) = ZERO

                    ! Accumulate hb coupling for each field
                    do ifield = 1,nfields
                        irow_start = 1 + (ifield-1)*nterms
                        irow_end   = irow_start + (nterms-1)
                        icol_start = 1 + (ifield-1)*nterms
                        icol_end   = icol_start + (nterms-1)
                        hb_mat(irow_start:irow_end,icol_start:icol_end) = self%D(itime_a,itime_b)*A%dom(1)%lblks(ielem,itime_a)%mass
                    end do  ! ifield

                    ! Store hb coupling for (itime_a,itime_b)
                    sz = size(A%dom(1)%lblks(ielem,itime_a)%data_(1)%mat,1)
                    irow_start = 1 + (itime_a-1)*sz
                    irow_end = irow_start + (sz-1)
                    icol_start = 1 + (itime_b-1)*sz
                    icol_end = icol_start + (sz-1)
                    self%hb_matrix%elems(ielem)%mat(irow_start:irow_end,icol_start:icol_end) = hb_mat

                end do !itime_b
            end do !itime_a

            ! Invert and store
            self%hb_matrix%elems(ielem)%mat = inv(self%hb_matrix%elems(ielem)%mat)

        end do !ielem

    end subroutine update
    !*******************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/19/2018
    !!
    !-------------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_HB_t),     intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: v

        type(chidg_vector_t)         :: z

        integer(ik)             :: ielem, nterms, nfields, ntime, ierr, idom, istart, iend, itime
        real(rk),   allocatable :: big_vec(:), new_vec(:)

        ! Copy v to z
        z = v

        ! Assume 1 Domain
        idom = 1

        do ielem = 1,size(A%dom(idom)%lblks,1)
            nterms  = A%dom(idom)%lblks(ielem,1)%data_(1)%nterms
            nfields = A%dom(idom)%lblks(ielem,1)%data_(1)%nfields
            ntime   = size(A%dom(idom)%lblks,2)

            ! Allocate multi-time vector
            if (allocated(big_vec)) deallocate(big_vec)
            allocate(big_vec(nterms*nfields*ntime), stat=ierr)
            if (ierr /= 0) call AllocationError

            ! Assemble multi-time vector
            do itime = 1,ntime
                istart = 1 + (itime-1)*nterms*nfields
                iend = istart + (nterms*nfields-1)
                big_vec(istart:iend) = z%dom(idom)%vecs(ielem)%gettime(itime)
            end do !itime

            ! Apply hb preconditioner
            new_vec = matmul(self%hb_matrix%elems(ielem)%mat,big_vec)

            ! Store multi-time vector
            ! Assemble multi-time vector
            do itime = 1,ntime
                istart = 1 + (itime-1)*nterms*nfields
                iend = istart + (nterms*nfields-1)
                call z%dom(idom)%vecs(ielem)%settime(itime,new_vec(istart:iend))
            end do !itime

        end do !ielem


    end function apply
    !-----------------------------------------------------------------------------------------








end module precon_HB