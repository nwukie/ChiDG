module precon_line
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
    use type_densematrix,       only: densematrix_t
    use type_element,           only: element_t
    use type_mesh,              only: mesh_t
    use type_ivector,           only: ivector_t
    implicit none





    !>  One row of a particular line.
    !!
    !!  [...                       ]
    !!  [   ...                    ]
    !!  [       ...                ]
    !!  [       B   A   C          ]  <-- This is one row
    !!  [               ...        ]
    !!  [                   ...    ]
    !!  [                       ...]
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2018
    !!
    !--------------------------------------------------------------------
    type, public :: row_t
        ! Matrix storage
        type(densematrix_t) :: A    ! Diagonal
        type(densematrix_t) :: B    ! Lower off-diagonal
        type(densematrix_t) :: C    ! Upper off-diagonal
        ! Block access indices in the lhs matrix densematrix_vector 
        integer(ik) :: A_index
        integer(ik) :: B_index
        integer(ik) :: C_index
    end type row_t
    !********************************************************************






    !>  A line data type that holds the line matrix data and applies the
    !!  effect of the solve on a constrained portion of the vector
    !!  via the block Thomas algorithm.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2018
    !!
    !--------------------------------------------------------------------
    type, public :: line_t
        ! List of elements coupled together
        type(ivector_t) :: idomain_g
        type(ivector_t) :: ielement_g
        type(ivector_t) :: idomain_l
        type(ivector_t) :: ielement_l

        ! One row for each element in the line
        type(row_t),    allocatable :: rows(:)
    contains
        procedure   :: init => init_line
        procedure   :: update => update_line
        procedure   :: apply => apply_line
    end type line_t
    !********************************************************************



    !>  Line preconditioner.
    !!
    !!  Create a group of lines. Each line solves its own block-tridiagonal subproblem
    !!  and applies this result to the incoming vector as a preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2018
    !!
    !----------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_line_t
        type(line_t),   allocatable :: lines(:)
    contains
        procedure   :: init   => init_preconditioner
        procedure   :: update => update_preconditioner
        procedure   :: apply  => apply_preconditioner
    end type precon_line_t
    !**********************************************************************************


contains

    
    !>  Initialize line preconditioner.
    !!
    !!  Allocating 1 line per every three elements.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2018
    !!
    !------------------------------------------------------------------------
    subroutine init_preconditioner(self,data)
        class(precon_line_t),   intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        integer(ik) :: idom, ielem, nelements, nlines, ierr, elem_count, line_ID


        ! Accumulate total number of elements. 
        nelements = 0
        do idom = 1,data%sdata%lhs%ndomains() 
            nelements = nelements + data%sdata%lhs%dom(idom)%nelements()
        end do !idom

        
        ! Compute number of lines 
        nlines = nint(real(nelements,rk)/real(3,rk))


        ! Loop through and count how many lines there will be
        nlines = 0
        elem_count = 1
        do idom = 1,data%sdata%lhs%ndomains()
            do ielem = 1,data%sdata%lhs%dom(idom)%nelements()
                ! Always put line on first element
                if ( (ielem==1) .or. (elem_count==2) ) then
                    nlines = nlines + 1
                    elem_count = 0
                end if
                elem_count = elem_count + 1
            end do !ielem
        end do !idom


        ! Allocate lines
        if (allocated(self%lines)) deallocate(self%lines)
        allocate(self%lines(nlines), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initialize lines
        line_ID = 0
        elem_count = 1
        do idom = 1,data%sdata%lhs%ndomains()
            do ielem = 1,data%sdata%lhs%dom(idom)%nelements()
                ! Always put line on first element
                if ( (ielem==1) .or. (elem_count==2) ) then
                    line_ID = line_ID + 1
                    elem_count = 0
                    ! Initialize line
                    call self%lines(line_ID)%init(data,idom,ielem)
                end if
                elem_count = elem_count + 1
            end do !ielem
        end do !idom

        
        ! Store spectral matrix and indicate initialization
        self%initialized = .true.

    end subroutine init_preconditioner
    !*****************************************************************************


    


    !>  Initialize line with seed (idomain_l, ielement_l)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/24/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine init_line(self,data,idomain_l,ielement_l)
        class(line_t),      intent(inout)   :: self
        type(chidg_data_t), intent(in)      :: data
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_l

        integer(ik) :: nrows, irow, max_rows, idom_search, ielem_search, itime, ierr, &
                       ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l
        logical :: neighbor_exists

        ! Start at the seed and travel in the XI_MAX direction to construct line path
        nrows = 1 ! At least 1 row from seed.
        idom_search = idomain_l
        ielem_search = ielement_l
        do while (nrows < max_rows) 
            ! Check neighbor exists, if exists then increment total number of rows
            neighbor_exists = (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l /= 0)
            if (neighbor_exists) nrows = nrows + 1
        end do

        ! Allocate rows
        if (allocated(self%rows)) deallocate(self%rows)
        allocate(self%rows(nrows), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Start at the seed and travel in the XI_MAX direction to construct line path
        nrows = 1 ! At least 1 row from seed.
        idom_search = idomain_l
        ielem_search = ielement_l
        itime = 1
        do irow = 1,nrows   ! row 1 was already started

            ! Store diagonal
            self%rows(irow)%A_index = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%get_diagonal()
            self%rows(irow)%A = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%at(self%rows(irow)%A_index)

            ! If not first row, find/store coupling of (search element) with (up-line neighbor element)
            if (irow /= 1) then
                self%rows(irow)%B_index = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%loc(self%rows(irow-1)%A%dparent_g(),self%rows(irow-1)%A%eparent_g(),itime)
                self%rows(irow)%B       = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%at(self%rows(irow)%C_index)
            end if

            ! Check down-line element exists, if exists then initialize row
            if (irow /= nrows) then
                neighbor_exists = (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l /= 0)
                if (neighbor_exists) then
                    ! Get neighbor indices
                    ineighbor_domain_g = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_domain_g
                    ineighbor_domain_l = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_domain_l
                    ineighbor_element_g = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_g
                    ineighbor_element_l = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l

                    ! Find/store coupling of (search element) with (down-line neighbor element)
                    self%rows(irow)%C_index = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%loc(ineighbor_domain_g,ineighbor_element_g,itime)
                    self%rows(irow)%C       = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%at(self%rows(irow)%C_index)

                    ! Update search element to be down-line neighbor
                    idom_search  = ineighbor_domain_l
                    ielem_search = ineighbor_element_l
                end if
            end if
        end do


    end subroutine init_line
    !******************************************************************************







    !>  Update the line preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine update_preconditioner(self,A,b)
        class(precon_line_t),   intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: b

        integer(ik) :: iline

        ! Each line, update its matrix entries
        do iline = 1,size(self%lines)
            call self%lines(iline)%update(A)
        end do !iline

    end subroutine update_preconditioner
    !*****************************************************************************





    !>  Update line matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/25/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine update_line(self,lhs)
        class(line_t),          intent(inout)   :: self
        class(chidg_matrix_t),  intent(in)      :: lhs

        integer(ik) :: idom, ielem, irow, itime

        itime = 1
        do irow = 1,size(self%rows)
            idom = self%rows(irow)%A%dparent_l()
            ielem = self%rows(irow)%A%eparent_l()
            self%rows(irow)%A = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%A_index)
            if (irow /= 1) self%rows(irow)%B = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%B_index)
            if (irow /= size(self%rows)) self%rows(irow)%C = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%C_index)
        end do !irow


    end subroutine update_line
    !*****************************************************************************









    !>  Apply the line preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2018
    !!
    !-----------------------------------------------------------------------------
    function apply_preconditioner(self,A,v) result(z)
        class(precon_line_t),   intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: v

        type(chidg_vector_t)    :: z

        integer(ik) :: iline

        ! Allocate z
        z = v

        ! Each line, apply its effect on v and store in z
        do iline = 1,size(self%lines)
            call self%lines(iline)%apply(A,v,z)
        end do !iline

    end function apply_preconditioner
    !*****************************************************************************







    !>  Apply the line contribution.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine apply_line(self,A,v,z)
        class(line_t),          intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: v
        type(chidg_vector_t),   intent(in)      :: z

        ! assemble rhs components into densevec
        itime = 1
        do irow = 1,size(self%rows)

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
                new_vec = matmul(self%hb_matrix(idom)%elems(ielem)%mat,big_vec)

                ! Store multi-time vector
                do itime = 1,ntime
                    istart = 1 + (itime-1)*nterms*nfields
                    iend = istart + (nterms*nfields-1)
                    call z%dom(idom)%vecs(ielem)%settime(itime,new_vec(istart:iend))
                end do !itime

        end do !idom






    end subroutine apply_line
    !*****************************************************************************










!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !----------------------------------------------------------------------------------
!    subroutine init_hb_dense_matrix(self,element)
!        class(hb_dense_matrix_t), intent(inout)   :: self
!        type(element_t),    intent(in)      :: element
!
!        integer(ik) :: ierr, mat_size
!
!        mat_size = element%nterms_s*element%neqns*element%ntime 
!
!        allocate(self%mat(mat_size,mat_size), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!    end subroutine init_hb_dense_matrix
!    !**********************************************************************************
!
!
!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !----------------------------------------------------------------------------------
!    subroutine init_hb_matrix(self,mesh)
!        class(hb_matrix_t), intent(inout)   :: self
!        type(mesh_t),       intent(in)      :: mesh
!
!        integer(ik) :: ierr, ielem
!
!        allocate(self%elems(mesh%domain(1)%nelements()), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        do ielem = 1,mesh%domain(1)%nelements()
!            call self%elems(ielem)%init(mesh%domain(1)%elems(ielem))
!        end do
!
!    end subroutine init_hb_matrix
!    !**********************************************************************************
!
!    
!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !----------------------------------------------------------------------------------
!    subroutine clear_hb_matrix(self)
!        class(hb_matrix_t), intent(inout)   :: self
!
!        integer(ik) :: ielem
!
!        do ielem = 1,size(self%elems)
!            self%elems(ielem)%mat = ZERO
!        end do
!
!    end subroutine clear_hb_matrix
!    !**********************************************************************************
!
!
!
!
!    !>  Initialize the HB preconditioner. This is for allocating storage. In this case, 
!    !!  we allocate a Lower-Diagonal block matrix for storing the LU decomposition.
!    !!  
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !!  @param[inout]   domain  domain_t instance containing a mesh component used to 
!    !!                          initialize the block matrix
!    !!
!    !-------------------------------------------------------------------------------------------
!    subroutine init(self,data)
!        class(precon_hb_t),     intent(inout)   :: self
!        type(chidg_data_t),     intent(in)      :: data
!
!        integer(ik) :: idom, ierr
!
!        ! Allocate domain storage
!        if (allocated(self%hb_matrix)) deallocate(self%hb_matrix)
!        allocate(self%hb_matrix(data%sdata%lhs%ndomains()), stat=ierr)
!        if (ierr /= 0) call AllocationError
!        
!        ! Initialize domain storage
!        do idom = 1,data%sdata%lhs%ndomains()
!            call self%hb_matrix(idom)%init(data%mesh)
!            call self%hb_matrix(idom)%clear()
!        end do !idom
!
!        ! Store spectral matrix and indicate initialization
!        self%D = data%time_manager%D
!        self%initialized = .true.
!
!    end subroutine init
!    !*******************************************************************************************
!
!
!
!
!
!
!    !>  Compute the block diagonal inversion and store so it can be applied.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !-------------------------------------------------------------------------------------
!    subroutine update(self,A,b)
!        class(precon_hb_t),     intent(inout)   :: self
!        type(chidg_matrix_t),   intent(in)      :: A
!        type(chidg_vector_t),   intent(in)      :: b
!
!        integer(ik) :: ielem, itime, sz, irow_start, irow_end, icol_start, icol_end, &
!                       itime_a, itime_b, idiag, ierr, ifield, nfields, nterms, idom
!        real(rk), allocatable   :: hb_mat(:,:)
!
!        do idom = 1,size(self%hb_matrix)
!            do ielem = 1,size(self%hb_matrix(idom)%elems)
!
!                ! Assemble diagonal elements
!                do itime = 1,size(A%dom(idom)%lblks,2)
!
!                    ! Access indices and size
!                    idiag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
!                    sz = size(A%dom(idom)%lblks(ielem,itime)%data_(idiag)%mat,1)
!                    irow_start = 1 + (itime-1)*sz
!                    irow_end = irow_start + (sz-1)
!                    icol_start = 1 + (itime-1)*sz
!                    icol_end = icol_start + (sz-1)
!
!                    ! Get diagonal index and store
!                    self%hb_matrix(idom)%elems(ielem)%mat(irow_start:irow_end,icol_start:icol_end) = A%dom(idom)%lblks(ielem,itime)%data_(idiag)%mat
!
!                end do !itime
!
!                ! Assemble off-diagonal coupling 
!                do itime_a = 1,size(A%dom(idom)%lblks,2)
!                    do itime_b = 1,size(A%dom(idom)%lblks,2)
!
!                        idiag   = A%dom(idom)%lblks(ielem,itime_a)%get_diagonal()
!                        nterms  = A%dom(idom)%lblks(ielem,itime_a)%data_(idiag)%nterms
!                        nfields = A%dom(idom)%lblks(ielem,itime_a)%data_(idiag)%nfields
!
!                        if (itime_a /= itime_b) then
!
!                            ! LHS HB contribution for each variable
!                            if (allocated(hb_mat)) deallocate(hb_mat)
!                            allocate(hb_mat(nterms*nfields,nterms*nfields), stat=ierr)
!                            if (ierr /= 0) call AllocationError
!                            hb_mat(:,:) = ZERO
!
!                            ! Accumulate hb coupling for each field
!                            do ifield = 1,nfields
!                                irow_start = 1 + (ifield-1)*nterms
!                                irow_end   = irow_start + (nterms-1)
!                                icol_start = 1 + (ifield-1)*nterms
!                                icol_end   = icol_start + (nterms-1)
!                                hb_mat(irow_start:irow_end,icol_start:icol_end) = self%D(itime_a,itime_b)*A%dom(idom)%lblks(ielem,itime_a)%mass
!                            end do  ! ifield
!
!                            ! Store hb coupling for (itime_a,itime_b)
!                            sz = size(A%dom(idom)%lblks(ielem,itime_a)%data_(1)%mat,1)
!                            irow_start = 1 + (itime_a-1)*sz
!                            irow_end = irow_start + (sz-1)
!                            icol_start = 1 + (itime_b-1)*sz
!                            icol_end = icol_start + (sz-1)
!                            self%hb_matrix(idom)%elems(ielem)%mat(irow_start:irow_end,icol_start:icol_end) = hb_mat
!
!                        end if
!
!                    end do !itime_b
!                end do !itime_a
!
!                ! Invert and store
!                self%hb_matrix(idom)%elems(ielem)%mat = inv(self%hb_matrix(idom)%elems(ielem)%mat)
!
!            end do !ielem
!        end do !idom
!
!    end subroutine update
!    !*******************************************************************************************
!
!
!
!
!
!
!
!
!    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/19/2018
!    !!
!    !-------------------------------------------------------------------------------------------
!    function apply(self,A,v) result(z)
!        class(precon_hb_t),     intent(inout)   :: self
!        type(chidg_matrix_t),   intent(in)      :: A
!        type(chidg_vector_t),   intent(in)      :: v
!
!        type(chidg_vector_t)    :: z
!
!        integer(ik)             :: ielem, nterms, nfields, ntime, ierr, idom, istart, iend, itime
!        real(rk),   allocatable :: big_vec(:), new_vec(:)
!
!        ! Copy v to z
!        z = v
!
!        ! Loop domains/elements and apply inverted HB-diagonal
!        do idom = 1,size(self%hb_matrix)
!            do ielem = 1,size(self%hb_matrix(idom)%elems)
!
!                nterms  = A%dom(idom)%lblks(ielem,1)%data_(1)%nterms
!                nfields = A%dom(idom)%lblks(ielem,1)%data_(1)%nfields
!                ntime   = size(A%dom(idom)%lblks,2)
!
!                ! Allocate multi-time vector
!                if (allocated(big_vec)) deallocate(big_vec)
!                allocate(big_vec(nterms*nfields*ntime), stat=ierr)
!                if (ierr /= 0) call AllocationError
!
!                ! Assemble multi-time vector
!                do itime = 1,ntime
!                    istart = 1 + (itime-1)*nterms*nfields
!                    iend = istart + (nterms*nfields-1)
!                    big_vec(istart:iend) = z%dom(idom)%vecs(ielem)%gettime(itime)
!                end do !itime
!
!                ! Apply hb preconditioner
!                new_vec = matmul(self%hb_matrix(idom)%elems(ielem)%mat,big_vec)
!
!                ! Store multi-time vector
!                do itime = 1,ntime
!                    istart = 1 + (itime-1)*nterms*nfields
!                    iend = istart + (nterms*nfields-1)
!                    call z%dom(idom)%vecs(ielem)%settime(itime,new_vec(istart:iend))
!                end do !itime
!
!            end do !ielem
!        end do !idom
!
!
!    end function apply
!    !-----------------------------------------------------------------------------------------








end module precon_line
