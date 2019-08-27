module mod_tioga
    use mod_tioga_interfaces,   only: tioga_init_f90
    use mpi_f08,                only: mpi_comm
    implicit none



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/27/2019
    !!
    !---------------------------------------------------------------
    subroutine tioga_compute_iblank(gridfile,comm)
        character(*),   intent(in)  :: gridfile
        type(mpi_comm), intent(in)  :: comm

        call tioga_init_f90(comm)

    end subroutine tioga_compute_iblank
    !***************************************************************






end module mod_tioga
