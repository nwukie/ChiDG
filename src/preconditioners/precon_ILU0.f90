module precon_ILU0
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_chidgVector
    use type_blockmatrix,       only: blockmatrix_t

    use mod_inv,                only: inv

    use mod_chidg_mpi,          only: IRANK
    implicit none


    !> ILU0 preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_ILU0_t

        type(blockmatrix_t), allocatable     :: LD(:)       !< Lower-Diagonal, sparse-block matrix representation

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply

    end type precon_ILU0_t
    !*****************************************************************************************************************




contains



    !> Initialize the ILU0 preconditioner. This is for allocating storage. In this case, we allocate
    !! a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  @param[inout]   domain      domain_t instance containing a mesh component used to initialize the block matrix
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        integer(ik) :: idom, ndom, ierr

        ndom = data%ndomains()

        !
        ! Allocate a Lower-Diagonal block matrix for each domain
        !
        allocate(self%LD(ndom), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Initialize each blockmatrix
        !
        do idom = 1,ndom
            call self%LD(idom)%init(mesh=data%mesh(idom),mtype='LowerDiagonal')
            call self%LD(idom)%clear()
        end do


        self%initialized = .true.

    end subroutine init
    !*****************************************************************************************************************








    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: b


        integer(ik)             :: ielem, irow, icol, eparent_l, iblk, idom, ndom
        integer(ik)             :: lower_blocks(3), upper_blocks(3)
        real(rk), allocatable   :: pdiag(:,:)


        call write_line(' Computing ILU0 factorization', io_proc=GLOBAL_MASTER)

        !
        ! Test preconditioner initialization
        !
        if ( .not. self%initialized ) call chidg_signal(FATAL,'preconditioner::ILU0%update - preconditioner has not yet been initialized')


        !
        ! For each domain
        !
        ndom = size(A%dom)
        do idom = 1,ndom


            !
            ! Store diagonal blocks of A
            !
            do ielem = 1,size(A%dom(idom)%lblks,1)
                self%LD(idom)%lblks(ielem,DIAG)%mat = A%dom(idom)%lblks(ielem,DIAG)%mat
            end do

            !
            ! Invert first diagonal
            !
            self%LD(idom)%lblks(1,DIAG)%mat = inv(self%LD(idom)%lblks(1,DIAG)%mat)





            lower_blocks = [XI_MIN, ETA_MIN, ZETA_MIN]
            upper_blocks = [XI_MAX, ETA_MAX, ZETA_MAX]

            !
            ! Loop through all rows
            !
            do irow = 2,size(A%dom(idom)%lblks,1)


                !
                ! Operate on all the L blocks for the current row
                !
                do iblk = 1,size(lower_blocks)

                    if (allocated(self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat) .and. self%LD(idom)%lblks(irow,lower_blocks(iblk))%parent_proc() == IRANK) then
                    !if (allocated(self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat) ) then


                        ! Get parent index
                        eparent_l = self%LD(idom)%lblks(irow,lower_blocks(iblk))%eparent_l()

                        ! Compute and store the contribution to the lower-triangular part of LD
                        self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat = matmul(A%dom(idom)%lblks(irow,lower_blocks(iblk))%mat,self%LD(idom)%lblks(eparent_l,DIAG)%mat)

                        ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                        self%LD(idom)%lblks(irow,DIAG)%mat = self%LD(idom)%lblks(irow,DIAG)%mat  -  matmul(self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat,  A%dom(idom)%lblks(eparent_l,upper_blocks(iblk))%mat)

                    end if
                end do ! iblk


                !
                ! Pre-Invert current diagonal block and store
                !
                self%LD(idom)%lblks(irow,DIAG)%mat = inv(self%LD(idom)%lblks(irow,DIAG)%mat)


            end do !irow

        end do ! idom

    end subroutine update
    !*******************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: v

        type(chidgVector_t)         :: z

        integer(ik)                 :: ielem, eparent_l, irow, iblk, block_index, idom, ndom
        integer(ik), allocatable    :: lower_blocks(:), upper_blocks(:)

        !
        ! Initialize z for preconditioning
        !
        z = v


        lower_blocks = [XI_MIN, ETA_MIN, ZETA_MIN]
        upper_blocks = [XI_MAX, ETA_MAX, ZETA_MAX]


        !
        ! For each domain
        !
        ndom = size(A%dom)
        do idom = 1,ndom


            !
            ! Forward Solve
            !
            do irow = 1,size(self%LD(idom)%lblks,1)


                !
                ! Lower-Triangular blocks
                !
                do block_index = 1,size(lower_blocks)
                    iblk = lower_blocks(block_index)

                    if (allocated(self%LD(idom)%lblks(irow,iblk)%mat) .and. self%LD(idom)%lblks(irow,iblk)%parent_proc() == IRANK) then
                    !if (allocated(self%LD(idom)%lblks(irow,iblk)%mat) ) then
                        !
                        ! Get associated parent block index
                        !
                        eparent_l = self%LD(idom)%lblks(irow,iblk)%eparent_l()

                        z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(self%LD(idom)%lblks(irow,iblk)%mat, z%dom(idom)%vecs(eparent_l)%vec)

                    end if


                end do


            end do ! irow



            !
            ! Backward Solve
            !
            do irow = size(A%dom(idom)%lblks,1),1,-1

                !
                ! Upper-Triangular blocks
                !
                do block_index = 1,size(upper_blocks)
                    iblk = upper_blocks(block_index)

                    if (allocated(A%dom(idom)%lblks(irow,iblk)%mat) .and. A%dom(idom)%lblks(irow,iblk)%parent_proc() == IRANK) then
                    !if (allocated(A%dom(idom)%lblks(irow,iblk)%mat) ) then
                        !
                        ! Get associated parent block index
                        !
                        eparent_l = A%dom(idom)%lblks(irow,iblk)%eparent_l()

                        z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(A%dom(idom)%lblks(irow,iblk)%mat, z%dom(idom)%vecs(eparent_l)%vec)

                    end if

                end do


                !
                ! Diagonal block
                !
                z%dom(idom)%vecs(irow)%vec = matmul(self%LD(idom)%lblks(irow,DIAG)%mat, z%dom(idom)%vecs(irow)%vec)

            end do ! irow



        end do ! idom




    end function apply
    !-----------------------------------------------------------------------------------------

















end module precon_ILU0
