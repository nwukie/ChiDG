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
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use mpi_f08
    implicit none



    
    integer(ik) :: IRANK = 0        !< Rank of current process
    integer(ik) :: NRANK = 1        !< Number of global ranks

    integer, parameter :: GLOBAL_MASTER = 0    !< Master rank for all global processes. This shall not be modified during run-time.
    integer            :: GROUP_MASTER  = 0    !< Master rank for group of processes. This could be modified during run-time for a group.


    logical :: chidg_mpi_initialized = .false.
    logical :: chidg_mpi_finalized   = .false.

contains








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine chidg_mpi_init()

        integer :: ierr

        !
        ! Initialize MPI
        !
        if ( .not. chidg_mpi_initialized ) then
            call MPI_Init(ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"MPI_Init")

            chidg_mpi_initialized = .true.
        end if



        call MPI_Comm_Size(MPI_COMM_WORLD,NRANK,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"MPI_Comm_Size")
        call MPI_Comm_Rank(MPI_COMM_WORLD,IRANK,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"MPI_Comm_Rank")


    end subroutine chidg_mpi_init
    !**********************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine chidg_mpi_finalize()

        integer :: ierr

        !
        ! Initialize MPI
        !
        if ( .not. chidg_mpi_finalized ) then
            call MPI_Finalize(ierr)
            if (ierr /= 0) call chidg_signal(FATAL,"MPI_Comm_Rank")

            chidg_mpi_finalized = .true.
        end if


    end subroutine chidg_mpi_finalize
    !**********************************************************************************






















end module mod_chidg_mpi
