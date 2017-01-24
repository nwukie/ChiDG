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
!#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use mpi_f08
    implicit none


    type(mpi_comm)      :: ChiDG_COMM          !< Communicator for ChiDG. of MPI_COMM_WORLD, 
                                               !! but could be changed at run-time.
    
    integer(ik)         :: IRANK = 0           !< Rank of current process
    integer(ik)         :: NRANK = 1           !< Number of global ranks

    integer, parameter  :: GLOBAL_MASTER = 0   !< Master rank for all global processes. This 
                                               !! shall not be modified during run-time.
    integer             :: GROUP_MASTER  = 0   !< Master rank for group of processes. This 
                                               !! could be modified during run-time for a group.


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

            chidg_mpi_initialized = .true.
        end if



        call MPI_Comm_Size(MPI_COMM_WORLD,NRANK,ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD,IRANK,ierr)


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

            chidg_mpi_finalized = .true.
        end if


    end subroutine chidg_mpi_finalize
    !**********************************************************************************













end module mod_chidg_mpi
