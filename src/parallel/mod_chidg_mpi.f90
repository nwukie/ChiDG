!> Module for storing mpi-based data, such as current rank, number of global ranks
!! and any other pertinent mpi datas. Maybe also some helpful routines for managing
!! mpi data and communication.
!!
!!  @author Nathan A. Wukie
!!
!!
!!
!-----------------------------------------------------------------------------
module mod_chidg_mpi
#include <petsc/finclude/petsc.h>
    use petsc,                          only: PETSC_NULL_CHARACTER, PETSC_COMM_WORLD
    use mod_kinds,                      only: rk, ik
    use mpi_f08
    implicit none


    type(mpi_comm)      :: ChiDG_COMM          ! Communicator for ChiDG. of MPI_COMM_WORLD, 
                                               ! but could be changed at run-time.
    
    integer(ik)         :: IRANK = 0           ! Rank of current process
    integer(ik)         :: NRANK = 1           ! Number of global ranks

    integer, parameter  :: GLOBAL_MASTER = 0   ! Master rank for all global processes. This 
                                               ! shall not be modified during run-time.
    integer             :: GROUP_MASTER  = 0   ! Master rank for group of processes. This 
                                               ! could be modified during run-time for a group.

    logical :: mpi_is_initialized       = .false.
    logical :: mpi_is_finalized         = .false.
    logical :: mpi_initialized_by_chidg = .false.


contains








    !>  Initialize MPI library. Get NRANK, IRANK.
    !!
    !!      - If already initialized, skip call to MPI_Init
    !!      - Query MPI library for NRANK, IRANK
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !----------------------------------------------------------------------------------
    subroutine chidg_mpi_init(comm)
        type(mpi_comm), intent(in), optional    :: comm

        PetscErrorCode  :: ierr

        ! Check if MPI_Init has been called already or by someone else
        call MPI_Initialized(mpi_is_initialized,ierr)


        ! Initialize MPI
        if ( .not. mpi_is_initialized ) then
            call MPI_Init(ierr)
            mpi_initialized_by_chidg = .true.
        else
            mpi_initialized_by_chidg = .false.
        end if


        ! Option to use incoming communicator instead of MPI_COMM_WORLD
        if (present(comm)) then
            ChiDG_COMM = comm
        else
            ChiDG_COMM = MPI_COMM_WORLD
        end if


        ! Initialize IRANK, NRANK 
        call MPI_Comm_Size(ChiDG_COMM,NRANK,ierr)
        call MPI_Comm_Rank(ChiDG_COMM,IRANK,ierr)


        ! Initialize PETSc
        PETSC_COMM_WORLD = ChiDG_COMM%mpi_val
        call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
        if (ierr .ne. 0) then
            print*,'Unable to initialize PETSc'
            stop
        end if

    end subroutine chidg_mpi_init
    !**********************************************************************************






    !>  Finalize MPI library. 
    !!
    !!      - If already finalized, skip call to MPI_Finalize
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !----------------------------------------------------------------------------------
    subroutine chidg_mpi_finalize()

        PetscErrorCode :: ierr


        ! Finalize PETSC
        call PetscFinalize(ierr)


        ! Only finalize if mpi was initialized by ChiDG. In tests,
        ! MPI is initialized by the test framework, so we count on it
        ! to also finalize and we don't finalize here.
        if (mpi_initialized_by_chidg) then

            ! Check if MPI_Finalize has been called already or by someone else
            call MPI_Finalized(mpi_is_finalized,ierr)
            if (ierr /= 0) then
                print*, '********************* WARNING *******************'
                print*, 'chidg_mpi_finalize: MPI_Finalized returned error.'
                print*, '*************************************************'
            end if

            ! Finalize MPI, if necessary
            if ( .not. mpi_is_finalized ) then
                call MPI_Finalize(ierr)
                if (ierr /= 0) then
                    print*, '********************* WARNING ******************'
                    print*, 'chidg_mpi_finalize: MPI_Finalize returned error.'
                    print*, '************************************************'
                end if
            end if

        end if


    end subroutine chidg_mpi_finalize
    !**********************************************************************************



end module mod_chidg_mpi
