module mod_hole_cutting
#include <messenger.h>
    use mod_tioga,  only: tioga_compute_iblank
    use mpi_f08,    only: mpi_comm
    use mod_io,     only: hole_cutter
    implicit none


contains


    !>  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/27/2019
    !!
    !!
    !------------------------------------------------------------
    subroutine compute_iblank(gridfile,comm)
        character(*),   intent(in)  :: gridfile
        type(mpi_comm), intent(in)  :: comm


        select case (trim(hole_cutter))
            case('tioga','TIOGA')
                call tioga_compute_iblank(gridfile,comm)
            case default
                call chidg_signal(FATAL,"mod_hole_cutting%compute_iblank: invalid value for 'hole_cutter'")
        end select

    end subroutine compute_iblank
    !*************************************************************





end module mod_hole_cutting
