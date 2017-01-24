!> Contains lengths of different data types. Short, long, etc.
!!
!!  @author Nathan A. Wukie
!!  @date   2/1/2016
!!
!-------------------------------------------------------------------
module mod_kinds
    implicit none

    private

    ! generic type parameters
    public :: ibyte, ishort, ilong
    public :: rsingle, rdouble
    public :: ik, rk    ! (ik - integer kind), (rk - real kind)

    ! integer kinds
    integer, parameter :: ibyte  = selected_int_kind(1) ! byte
    integer, parameter :: ishort = selected_int_kind(4) ! short
    integer, parameter :: ilong  = selected_int_kind(8) ! long

    ! floating point kinds
    integer, parameter :: rsingle = selected_real_kind(6,  37  )
    integer, parameter :: rdouble = selected_real_kind(15, 307 )
    integer, parameter :: rquad   = selected_real_kind(33, 4931)

    ! generic/default kinds
    integer, parameter :: ik = ilong     ! default integer kind
    integer, parameter :: rk = rdouble   ! default real kind
    !integer, parameter :: rk = rquad
    !integer, parameter :: rk = rsingle

    ! TECPLOT
    integer, parameter, public :: TEC = 4

end module mod_kinds
