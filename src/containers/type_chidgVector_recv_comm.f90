module type_chidgVector_recv_comm
#include <messenger.h>
    use mod_kinds,          only: ik
    use type_mesh,          only: mesh_t
    use type_ivector,       only: ivector_t
    use type_blockvector,   only: blockvector_t

    !test
    use mod_chidg_mpi,      only: IRANK
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: chidgVector_recv_comm_t

        integer(ik)                         :: proc
        type(blockvector_t),    allocatable :: dom(:)

    contains

        procedure,  public  :: init

    end type chidgVector_recv_comm_t
    !****************************************************************************








contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self,mesh,iproc,icomm)
        class(chidgVector_recv_comm_t), intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh(:)
        integer(ik),                    intent(in)      :: iproc
        integer(ik),                    intent(in)      :: icomm

        type(ivector_t)             :: dom_recv
        integer(ik)                 :: idom, idom_recv, ndom_recv, ierr
        integer(ik),    allocatable :: comm_procs_dom(:)
        logical                     :: proc_has_domain

        !
        ! Set processor being received from
        !
        self%proc = iproc


        !
        ! Compute number of domains being received from proc
        !
        do idom = 1,size(mesh)

            ! Get comm procs for domain
            comm_procs_dom = mesh(idom)%get_comm_procs()

            ! Is proc communicating any of the current domain
            proc_has_domain = any( iproc == comm_procs_dom )

            ! If so, add domain to list of domains being received from current_proc
            if (proc_has_domain) then
                call dom_recv%push_back(idom)
            end if

        end do ! idom



        !
        ! Allocate number of domains to recv from proc
        !
        ndom_recv = dom_recv%size()
        if (allocated(self%dom)) deallocate(self%dom)
        allocate(self%dom(ndom_recv), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each domain to be received
        !
        do idom_recv = 1,dom_recv%size()
            call self%dom(idom_recv)%init(mesh(dom_recv%at(idom_recv)),iproc,icomm,idom_recv)
        end do ! idom


    end subroutine init
    !*********************************************************************************




end module type_chidgVector_recv_comm
