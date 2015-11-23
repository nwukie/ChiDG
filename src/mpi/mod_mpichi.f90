!> Module for storing mpi-based data, such as current rank, number of global ranks
!! and any other pertinent mpi datas. Maybe also some helpful routines for managing
!! mpi data and communication.
!!
!!  @author Nathan A. Wukie
!!
!!
!!
!-----------------------------------------------------------------------------
module mod_mpichi
    use mod_kinds,  only: rk, ik
    implicit none



    
    integer(ik) :: irank        !< Rank of current process
    integer(ik) :: nrank        !< Number of global ranks

    integer, parameter :: GLOBAL_MASTER = 0  !< Master rank for all global processes
    !integer, parameter :: GROUP_MASTER  = 0  !< Master rank for group of processes



contains











end module mod_mpichi
