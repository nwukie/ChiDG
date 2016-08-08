module precon_ILU0
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, ONE
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK

    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_chidgVector
    implicit none


    !>  ILU0 preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_ILU0_t

        type(chidgMatrix_t)     :: LD

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
        class(precon_ILU0_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data



        call self%LD%init(mesh=data%mesh, mtype='LowerDiagonal')
        call self%LD%clear()


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
        class(precon_ILU0_t),    intent(inout)   :: self
        type(chidgMatrix_t),        intent(in)      :: A
        type(chidgVector_t),        intent(in)      :: b


        integer(ik)             :: ielem, irow, icol, eparent_l, idom, ndom, ilower, itranspose


        call write_line(' Computing ILU0 factorization', io_proc=GLOBAL_MASTER)




        !
        ! Test preconditioner initialization
        !
        if ( .not. self%initialized ) call chidg_signal(FATAL,'ILU0%update - preconditioner has not yet been initialized')



        !
        ! For each domain
        !
        ndom = size(A%dom)
        do idom = 1,ndom


            !
            ! Store diagonal blocks of A
            !
            do ielem = 1,size(A%dom(idom)%lblks,1)
                self%LD%dom(idom)%lblks(ielem,DIAG)%mat = A%dom(idom)%lblks(ielem,DIAG)%mat
            end do


            !
            ! Invert first diagonal block
            !
            self%LD%dom(idom)%lblks(1,DIAG)%mat = inv(self%LD%dom(idom)%lblks(1,DIAG)%mat)



            !
            ! Loop through all Proc-Local rows
            !
            do irow = 2,size(A%dom(idom)%lblks,1)


                !
                ! Operate on all the L blocks for the current row
                !
                do icol = 1,A%dom(idom)%local_lower_blocks(irow)%size()
                    ilower = A%dom(idom)%local_lower_blocks(irow)%at(icol)

                    if (allocated(A%dom(idom)%lblks(irow,ilower)%mat) .and. A%dom(idom)%lblks(irow,ilower)%parent_proc() == IRANK) then

                        ! Get parent index
                        eparent_l = self%LD%dom(idom)%lblks(irow,ilower)%eparent_l()

                        ! Compute and store the contribution to the lower-triangular part of LD
                        self%LD%dom(idom)%lblks(irow,ilower)%mat = matmul(A%dom(idom)%lblks(irow,ilower)%mat,self%LD%dom(idom)%lblks(eparent_l,DIAG)%mat)

                        ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                        itranspose = A%dom(idom)%local_transpose(irow,ilower)
                        self%LD%dom(idom)%lblks(irow,DIAG)%mat = self%LD%dom(idom)%lblks(irow,DIAG)%mat  -  matmul(self%LD%dom(idom)%lblks(irow,ilower)%mat,  A%dom(idom)%lblks(eparent_l,itranspose)%mat)

                    end if

                end do ! icol


                !
                ! Pre-Invert current diagonal block and store
                !
                self%LD%dom(idom)%lblks(irow,DIAG)%mat = inv(self%LD%dom(idom)%lblks(irow,DIAG)%mat)


            end do !irow



        end do ! idom

        call write_line(' Done Computing ILU0 factorization', io_proc=GLOBAL_MASTER)
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

        integer(ik)                 :: ielem, eparent_l, idom, ndom, irow, icol, ilower, iupper


        call self%timer%start()


        !
        ! Initialize z for preconditioning
        !
        z = v



        !
        ! For each domain
        !
        ndom = size(A%dom)
        do idom = 1,ndom


            !
            ! Forward Solve - Local
            !
            do irow = 1,size(self%LD%dom(idom)%lblks,1)

                !
                ! Lower-Triangular blocks
                !
                do icol = 1,A%dom(idom)%local_lower_blocks(irow)%size()
                    ilower = A%dom(idom)%local_lower_blocks(irow)%at(icol)

                    if (allocated(A%dom(idom)%lblks(irow,ilower)%mat) .and. A%dom(idom)%lblks(irow,ilower)%parent_proc() == IRANK) then
                        
                            ! Get associated parent block index
                            eparent_l = self%LD%dom(idom)%lblks(irow,ilower)%eparent_l()
                            z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(self%LD%dom(idom)%lblks(irow,ilower)%mat, z%dom(idom)%vecs(eparent_l)%vec)

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
                do icol = 1,A%dom(idom)%local_upper_blocks(irow)%size()
                    iupper = A%dom(idom)%local_upper_blocks(irow)%at(icol)

                    if (allocated(A%dom(idom)%lblks(irow,iupper)%mat) .and. A%dom(idom)%lblks(irow,iupper)%parent_proc() == IRANK) then

                            ! Get associated parent block index
                            eparent_l = A%dom(idom)%lblks(irow,iupper)%eparent_l()
                            z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(A%dom(idom)%lblks(irow,iupper)%mat, z%dom(idom)%vecs(eparent_l)%vec)

                    end if

                end do


                !
                ! Diagonal block
                !
                z%dom(idom)%vecs(irow)%vec = matmul(self%LD%dom(idom)%lblks(irow,DIAG)%mat, z%dom(idom)%vecs(irow)%vec)

            end do ! irow




        end do ! idom



        call self%timer%stop()

    end function apply
    !-----------------------------------------------------------------------------------------

















end module precon_ILU0
