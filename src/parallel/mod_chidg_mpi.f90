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
    use mod_kinds,  only: rk, ik
    implicit none



    
    integer(ik) :: IRANK = 1        !< Rank of current process
    integer(ik) :: NRANK = 2        !< Number of global ranks

    integer, parameter :: GLOBAL_MASTER = 0     !< Master rank for all global processes. This shall not be modified during run-time.
    !integer, parameter :: GROUP_MASTER  = 0    !< Master rank for group of processes. This could be modified during run-time for a group.



contains











end module mod_chidg_mpi
