module mod_communication
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_mesh,  only: mesh_t
    use mod_chidg_mpi,  only: IRANK, NRANK, GLOBAL_MASTER
    use mod_chimera,    only: detect_chimera_faces,  &
                              detect_chimera_donors, &
                              compute_chimera_interpolators
    use mpi_f08
    implicit none

contains




    !>  Establish element/face neighbor communication.
    !!
    !!  For a given mesh instance for a domain, attempt to find neighboring elements with which
    !!  to communicate. This searches both in the same domain on the current processor(local), 
    !!  and also across processors, where part of a domain may have been split onto another 
    !!  processor(global)
    !!
    !!  This does NOT establish Chimera communication.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!  @param[inout]   mesh        An array of mesh types on the local processor.
    !!  @param[in]      ChiDG_COMM  An mpi communicator of the relevant processors. 
    !!                              Particularly usefull for testing.
    !!
    !------------------------------------------------------------------------------------------
    subroutine establish_neighbor_communication(mesh,ChiDG_COMM)
        type(mesh_t),   intent(inout)   :: mesh
        type(mpi_comm), intent(in)      :: ChiDG_COMM

        ! All ranks initialize local communication
        call mesh%init_comm_local()

        ! All ranks initialize global communication
        call mesh%init_comm_global(ChiDG_COMM)

    end subroutine establish_neighbor_communication
    !*****************************************************************************************





    !>  Establish chimera face communication.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/22/2016
    !!
    !!
    !!  @param[inout]   mesh        An array of mesh types on the local processor.
    !!  @param[in]      ChiDG_COMM  An mpi communicator of the relevant processors. 
    !!                              Particularly useful for testing.
    !!
    !------------------------------------------------------------------------------------------
    subroutine establish_chimera_communication(mesh,ChiDG_COMM)
        type(mesh_t),   intent(inout)   :: mesh
        type(mpi_comm)                  :: ChiDG_COMM

        integer :: idom, ierr

        !
        ! In case establish_chimera_communication is being called as a reinitialization 
        ! procedure. Calling chimera%clear to wipe out previous data 
        ! and redetect all chimera faces, donors and reinitialize donor data.
        !
        call write_line("Initialize: chimera communication...", io_proc=GLOBAL_MASTER)
        do idom = 1,mesh%ndomains()
            call mesh%domain(idom)%chimera%clear()
        end do !idom
        call MPI_Barrier(ChiDG_COMM,ierr)   ! not sure if this is needed.

        ! Process-local routine, just to flag faces
        call detect_chimera_faces(mesh)

        ! Detect local or global chimera donors
        call detect_chimera_donors(mesh)

        ! Compute chimera interpolation matrices to evaluate donor 
        ! solution at receiver quadrature nodes.
        call compute_chimera_interpolators(mesh)

        ! Barrier
        call MPI_Barrier(ChiDG_COMM,ierr)


    end subroutine establish_chimera_communication
    !*****************************************************************************************










end module mod_communication
