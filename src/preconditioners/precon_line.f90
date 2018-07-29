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
    !!  [       A   B   C          ]  <-- This is one row
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
        type(densematrix_t) :: B    ! Diagonal 
        type(densematrix_t) :: A    ! Lower off-diagonal
        type(densematrix_t) :: C    ! Upper off-diagonal

        real(rk), allocatable   :: gam(:,:)
        real(rk), allocatable   :: beta(:)
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
        procedure   :: init   => init_line
        procedure   :: update => update_line
        procedure   :: apply  => apply_line
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

        integer(ik) :: seed_per_nelem = 1

        !print*, 'init preconditioner - 1'

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
                if ( (ielem==1) .or. (elem_count==seed_per_nelem) ) then
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
                if ( (ielem==1) .or. (elem_count==seed_per_nelem) ) then
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

        !print*, 'init preconditioner - 2'

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
                       ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l, ineighbor_proc
        logical :: neighbor_exists

        max_rows = 1
        !print*, 'init line - 1'

        ! Start at the seed and travel in the XI_MAX direction to construct line path
        nrows = 1 ! At least 1 row from seed.
        idom_search = idomain_l
        ielem_search = ielement_l
        neighbor_exists = ( (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l /= 0) .and. &
                            (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_proc == IRANK ) )
                            
        do while ( (nrows < max_rows) .and. neighbor_exists) 
            ! Get neighbor indices
            ineighbor_domain_l  = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_domain_l
            ineighbor_element_l = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l
            ineighbor_proc      = data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_proc

            ! Check neighbor exists, if exists then increment total number of rows
            neighbor_exists = (ineighbor_element_l /= 0 .and. ineighbor_proc == IRANK)

            ! If exists, increment and move downstream
            if (neighbor_exists) then
                nrows = nrows + 1
                idom_search  = ineighbor_domain_l
                ielem_search = ineighbor_element_l
            end if
        end do


        ! Allocate rows
        if (allocated(self%rows)) deallocate(self%rows)
        allocate(self%rows(nrows), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Start at the seed and travel in the XI_MAX direction to construct line path
        !nrows = 1 ! At least 1 row from seed.
        idom_search = idomain_l
        ielem_search = ielement_l
        itime = 1
        do irow = 1,nrows   ! row 1 was already started

            ! Store diagonal
            self%rows(irow)%B_index = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%get_diagonal()
            self%rows(irow)%B = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%at(self%rows(irow)%B_index)

            ! If not first row, find/store coupling of (search element) with (up-line neighbor element)
            if (irow /= 1) then
                self%rows(irow)%A_index = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%loc(self%rows(irow-1)%B%dparent_g(),self%rows(irow-1)%B%eparent_g(),itime)
                self%rows(irow)%A       = data%sdata%lhs%dom(idom_search)%lblks(ielem_search,itime)%at(self%rows(irow)%A_index)
            end if

            ! Check down-line element exists, if exists then initialize row
            if (irow /= nrows) then
                neighbor_exists = ( (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_element_l /= 0) .and. &
                                    (data%mesh%domain(idom_search)%faces(ielem_search,XI_MAX)%ineighbor_proc == IRANK) )

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


        !print*, 'init line - 2'

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

        !print*, 'update preconditioner - 1'

        ! Each line, update its matrix entries
        do iline = 1,size(self%lines)
            call self%lines(iline)%update(A)
        end do !iline

        !print*, 'update preconditioner - 2'

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

        real(rk),   allocatable :: tmp(:,:)

        !print*, 'update line - 1'

        ! Transfer A, B, C matrices from lhs to line, row-by-row
        itime = 1
        do irow = 1,size(self%rows)
            idom = self%rows(irow)%B%dparent_l()
            ielem = self%rows(irow)%B%eparent_l()
            self%rows(irow)%B = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%B_index)
            if (irow /= 1              ) self%rows(irow)%A = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%A_index)
            if (irow /= size(self%rows)) self%rows(irow)%C = lhs%dom(idom)%lblks(ielem,itime)%at(self%rows(irow)%C_index)
        end do !irow


        ! Precompute some quantities used in the forward/reverse sweep.
        ! Precompute B_prime and overwrite diagonal. Precompute gam
        irow = 1
        self%rows(irow)%B%mat = inv(self%rows(irow)%B%mat)  ! b_1 = b_prime_1 = inv(b_1)
        self%rows(irow)%gam   = ZERO*self%rows(irow)%B%mat  ! gam_1 = 0
        do irow = 2,size(self%rows)
            self%rows(irow)%gam = matmul(self%rows(irow-1)%B%mat,-self%rows(irow-1)%C%mat)
            tmp = (matmul(self%rows(irow)%A%mat,self%rows(irow)%gam) + self%rows(irow)%B%mat)
            self%rows(irow)%B%mat = inv(tmp)
        end do

        !print*, 'update line - 2'

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

        integer(ik) :: iline, idom, ielem, diag_index, itime

        type(densematrix_t) :: D


        ! Allocate/clear z
        z = v

        ! Each line, apply its effect on v and store in z
        do iline = 1,size(self%lines)
            call self%lines(iline)%apply(A,v,z)
        end do !iline


        ! Loop through and check for elements that got missed by lines
        ! and apply Jacobi
        do idom = 1,z%ndomains()
            do ielem = 1,z%dom(idom)%nelements()
                do itime = 1,z%get_ntime()
                    if (norm2(z%dom(idom)%vecs(ielem)%vec) < 1.e-8_rk) then
                        diag_index = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                        D = A%dom(idom)%lblks(ielem,itime)%at(diag_index)
                        z%dom(idom)%vecs(ielem)%vec = matmul(inv(D%mat),v%dom(idom)%vecs(ielem)%vec)
                    end if
                end do
            end do
        end do

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
        type(chidg_vector_t),   intent(inout)   :: z

        integer(ik) :: irow, idom, ielem
        real(rk),   allocatable :: x(:), x_prev(:), y_prev(:), tmp_vec(:)

        !print*, 'apply line - 1'

        ! Forward sweep: compute beta
        irow = 1
        self%rows(irow)%A%mat = ZERO*self%rows(irow)%b%mat
        self%rows(irow)%beta  = ZERO*self%rows(irow)%B%mat(:,1)  ! beta_1 = 0
        do irow = 2,size(self%rows)
            ! Get y(irow-1)
            idom = self%rows(irow-1)%B%dparent_l()
            ielem = self%rows(irow-1)%B%eparent_l()
            y_prev = v%dom(idom)%vecs(ielem)%vec
            ! Compute beta
            tmp_vec = matmul(self%rows(irow-1)%A%mat,self%rows(irow-1)%beta)
            self%rows(irow)%beta = matmul(self%rows(irow-1)%B%mat, y_prev - tmp_vec)
        end do !irow



        ! Reverse sweep: compute x

        ! Last element, stores and initializes x_prev for recursive loop
        x_prev = self%rows(size(self%rows))%beta
        idom  = self%rows(size(self%rows))%B%dparent_l()
        ielem = self%rows(size(self%rows))%B%eparent_l()
        z%dom(idom)%vecs(ielem)%vec = x_prev

        ! Store if not already filled
        do irow = size(self%rows)-1,1,-1
            x = matmul(self%rows(irow+1)%gam,x_prev) + self%rows(irow+1)%beta

            ! Store if not already filled from another line
            idom = self%rows(irow)%B%dparent_l()
            ielem = self%rows(irow)%B%eparent_l()
            z%dom(idom)%vecs(ielem)%vec = x

            ! Store x to x_prev
            x_prev = x

        end do !irow


        !print*, 'apply line - 2'

    end subroutine apply_line
    !*****************************************************************************








end module precon_line
