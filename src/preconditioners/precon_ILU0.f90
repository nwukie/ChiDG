module precon_ILU0
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_chidgVector
    use type_blockmatrix,       only: blockmatrix_t

    use mod_inv,    only: inv
    implicit none


    !> ILU0 preconditioner
    !!
    !!
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




contains



    !> Initialize the ILU0 preconditioner. This is for allocating storage. In this case, we allocate
    !! a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie
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
            call self%LD(idom)%init(data%mesh(idom),'LowerDiagonal')
            call self%LD(idom)%clear()
        end do


        self%initialized = .true.

    end subroutine








    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: b


        integer(ik)             :: ielem, irow, icol, eparent, iblk, idom, ndom
        integer(ik)             :: lower_blocks(3), upper_blocks(3)
        real(rk), allocatable   :: pdiag(:,:)


        call write_line(' Computing ILU0 factorization')

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
                    if (allocated(self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat)) then


                        ! Get parent index
                        eparent = self%LD(idom)%lblks(irow,lower_blocks(iblk))%eparent()

                        ! Compute and store the contribution to the lower-triangular part of LD
                        self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat = matmul(A%dom(idom)%lblks(irow,lower_blocks(iblk))%mat,self%LD(idom)%lblks(eparent,DIAG)%mat)

                        ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                        self%LD(idom)%lblks(irow,DIAG)%mat = self%LD(idom)%lblks(irow,DIAG)%mat  -  matmul(self%LD(idom)%lblks(irow,lower_blocks(iblk))%mat,  A%dom(idom)%lblks(eparent,upper_blocks(iblk))%mat)

                    end if
                end do ! iblk


                !
                ! Pre-Invert current diagonal block and store
                !
                self%LD(idom)%lblks(irow,DIAG)%mat = inv(self%LD(idom)%lblks(irow,DIAG)%mat)


            end do !irow

        end do ! idom

    end subroutine update








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: v

        type(chidgVector_t)         :: z

        integer(ik)                 :: ielem, eparent, irow, iblk, block_index, idom, ndom
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

                    if (allocated(self%LD(idom)%lblks(irow,iblk)%mat)) then
                        !
                        ! Get associated parent block index
                        !
                        eparent = self%LD(idom)%lblks(irow,iblk)%eparent()

                        z%dom(idom)%lvecs(irow)%vec = z%dom(idom)%lvecs(irow)%vec - matmul(self%LD(idom)%lblks(irow,iblk)%mat, z%dom(idom)%lvecs(eparent)%vec)

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

                    if (allocated(A%dom(idom)%lblks(irow,iblk)%mat)) then
                        !
                        ! Get associated parent block index
                        !
                        eparent = A%dom(idom)%lblks(irow,iblk)%eparent()

                        z%dom(idom)%lvecs(irow)%vec = z%dom(idom)%lvecs(irow)%vec - matmul(A%dom(idom)%lblks(irow,iblk)%mat, z%dom(idom)%lvecs(eparent)%vec)

                    end if

                end do


                !
                ! Diagonal block
                !
                !z%lvecs(irow)%vec = matmul(inv(self%LD%lblks(irow,DIAG)%mat), z%lvecs(irow)%vec)
                z%dom(idom)%lvecs(irow)%vec = matmul(self%LD(idom)%lblks(irow,DIAG)%mat, z%dom(idom)%lvecs(irow)%vec)

            end do ! irow



        end do ! idom




    end function apply
    !-----------------------------------------------------------------------------------------

















end module precon_ILU0
