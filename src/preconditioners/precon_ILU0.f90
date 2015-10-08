module precon_ILU0
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX
    use type_domain,            only: domain_t
    use type_preconditioner,    only: preconditioner_t
    use type_blockmatrix,       only: blockmatrix_t
    use type_blockvector,       only: blockvector_t

    use mod_inv,    only: inv
    implicit none


    !> ILU0 preconditioner
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_ILU0_t

        type(blockmatrix_t)                 :: LD       !< Lower-Diagonal, sparse-block matrix representation

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
    subroutine init(self,domain)
        class(precon_ILU0_t),   intent(inout)   :: self
        type(domain_t),         intent(inout)   :: domain


        call self%LD%init(domain%mesh,'LowerDiagonal')
        call self%LD%clear()


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
        type(blockmatrix_t),    intent(in)      :: A
        type(blockvector_t),    intent(in)      :: b


        integer(ik)             :: ielem, irow, icol, iparent, iblk
        integer(ik)             :: lower_blocks(3), upper_blocks(3)
        real(rk), allocatable   :: pdiag(:,:)


        print*, ' Computing ILU0 factorization'

        !
        ! Test preconditioner initialization
        !
        if ( .not. self%initialized ) call signal(FATAL,'preconditioner::ILU0%update - preconditioner has not yet been initialized')


        !
        ! Store diagonal blocks of A
        !
        do ielem = 1,size(A%lblks,1)
            self%LD%lblks(ielem,DIAG)%mat = A%lblks(ielem,DIAG)%mat
        end do




        lower_blocks = [XI_MIN, ETA_MIN, ZETA_MIN]
        upper_blocks = [XI_MAX, ETA_MAX, ZETA_MAX]

        !
        ! Loop through all rows
        !
        do irow = 2,size(A%lblks,1)


            !
            ! Operate on all the L blocks for the current row
            !
            do iblk = 1,size(lower_blocks)
                if (allocated(self%LD%lblks(irow,lower_blocks(iblk))%mat)) then


                    ! Invert parent diagonal from the preconditioner
                    iparent = self%LD%lblks(irow,lower_blocks(iblk))%parent()
                    pdiag = inv(self%LD%lblks(iparent,DIAG)%mat)

                    ! Compute and store the contribution to the lower-triangular part of LD
                    self%LD%lblks(irow,lower_blocks(iblk))%mat = matmul(A%lblks(irow,lower_blocks(iblk))%mat,pdiag)


                    ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                    self%LD%lblks(irow,DIAG)%mat = self%LD%lblks(irow,DIAG)%mat  -  matmul(self%LD%lblks(irow,lower_blocks(iblk))%mat,  A%lblks(iparent,upper_blocks(iblk))%mat)

                end if
            end do



        end do







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
        type(blockmatrix_t),    intent(in)      :: A
        type(blockvector_t),    intent(in)      :: v

        type(blockvector_t)         :: z

        integer(ik)                 :: ielem, iparent, irow, iblk, block_index
        integer(ik), allocatable    :: lower_blocks(:), upper_blocks(:)

        !
        ! Initialize z for preconditioning
        !
        z = v


        lower_blocks = [XI_MIN, ETA_MIN, ZETA_MIN]
        upper_blocks = [XI_MAX, ETA_MAX, ZETA_MAX]

        !
        ! Forward Solve
        !
        do irow = 1,size(self%LD%lblks,1)


            !
            ! Lower-Triangular blocks
            !
            do block_index = 1,size(lower_blocks)
                iblk = lower_blocks(block_index)

                if (allocated(self%LD%lblks(irow,iblk)%mat)) then
                    !
                    ! Get associated parent block index
                    !
                    iparent = self%LD%lblks(irow,iblk)%parent()

                    z%lvecs(irow)%vec = z%lvecs(irow)%vec - matmul(self%LD%lblks(irow,iblk)%mat, z%lvecs(iparent)%vec)

                end if


            end do




        end do



        !
        ! Backward Solve
        !
        do irow = size(A%lblks,1),1,-1

            !
            ! Upper-Triangular blocks
            !
            do block_index = 1,size(upper_blocks)
                iblk = upper_blocks(block_index)

                if (allocated(A%lblks(irow,iblk)%mat)) then
                    !
                    ! Get associated parent block index
                    !
                    iparent = A%lblks(irow,iblk)%parent()

                    z%lvecs(irow)%vec = z%lvecs(irow)%vec - matmul(A%lblks(irow,iblk)%mat, z%lvecs(iparent)%vec)

                end if

            end do


            !
            ! Diagonal block
            !

            z%lvecs(irow)%vec = matmul(inv(self%LD%lblks(irow,DIAG)%mat), z%lvecs(irow)%vec)

        end do








    end function apply

















end module precon_ILU0
