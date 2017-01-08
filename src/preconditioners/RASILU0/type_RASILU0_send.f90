module type_RASILU0_send
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_chidg_mpi,          only: IRANK

    use type_RASILU0_send_comm, only: RASILU0_send_comm_t
    use type_mesh,              only: mesh_t
    use type_chidg_matrix,       only: chidg_matrix_t
    use type_ivector,           only: ivector_t
    implicit none




    !>  A container for containing information about what to send to each neighboring processor.
    !!  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: RASILU0_send_t

        type(RASILU0_send_comm_t),  allocatable :: comm(:)  ! One description for each 
                                                            ! neighboring processor

    contains

        procedure   :: init

        procedure   :: init_wait

    end type RASILU0_send_t
    !******************************************************************************************



    



contains

    
    !>  Find neighboring processors, call initialization for each.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/22/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,A)
        class(RASILU0_send_t),  intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidg_matrix_t),    intent(in)      :: A

        type(ivector_t)             :: send_procs
        integer(ik),    allocatable :: send_procs_dom(:)
        integer(ik)                 :: idom, icomm, nsend_procs, ierr, iproc



        !
        ! Loop through each mesh and accumulate the total number of 
        ! processors we are sending to.
        !
        do idom = 1,size(mesh)
            send_procs_dom = mesh(idom)%get_send_procs_local()

            if ( allocated(send_procs_dom) ) then
                do iproc = 1,size(send_procs_dom)
                    call send_procs%push_back_unique(send_procs_dom(iproc))
                end do !iproc
            end if

        end do !idom



        nsend_procs = send_procs%size()
        allocate(self%comm(nsend_procs), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize communication for each processor
        !
        do icomm = 1,size(self%comm)
            iproc = send_procs%at(icomm)
            call self%comm(icomm)%init(mesh,A,iproc)
        end do



    end subroutine init
    !******************************************************************************************







    !>  Wall MPI_WAITALL on all outstanding nonblocking sends that were initiated during the
    !!  initialization procedure.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/08/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_wait(self)
        class(RASILU0_send_t),  intent(inout)   :: self
        
        integer(ik) :: icomm

        !
        ! Call wait on requests for each send_comm instance.
        !
        do icomm = 1,size(self%comm)

            call self%comm(icomm)%init_wait()

        end do !icomm


    end subroutine init_wait
    !******************************************************************************************











end module type_RASILU0_send
