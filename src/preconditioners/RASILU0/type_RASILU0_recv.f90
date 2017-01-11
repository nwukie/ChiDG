module type_RASILU0_recv
#include <messenger.h>
    use mod_kinds,              only: ik, rk

    use type_mesh,              only: mesh_t
    use type_chidg_matrix,       only: chidg_matrix_t
    use type_chidg_vector,       only: chidg_vector_t
    use type_RASILU0_recv_dom,  only: RASILU0_recv_dom_t
    implicit none






    !>  A data-type for holding matrix data from neighboring processors.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: RASILU0_recv_t

        type(RASILU0_recv_dom_t),   allocatable :: dom(:)   ! A container for each domain of the
                                                            ! local problem to accept neighboring
                                                            ! elements from another processor.
                                                            ! No chimera coupling.

    contains

        procedure   :: init

    end type RASILU0_recv_t
    !******************************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,A,b)
        class(RASILU0_recv_t),  intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: b

        integer(ik) :: idom, ndom, ierr, icomm_vec, idom_vec, ielem_vec, proc, &
                       iblk, icomm, idomain_g, ielem, &
                       dparent_g, eparent_g, dparent_g_vec, eparent_g_vec, parent_proc, &
                       trans_elem, trans_block
        logical     :: check_ok

        ndom = size(mesh)


        !
        ! Allocate overlapping storage for each domain to receive data into
        !
        allocate(self%dom(ndom), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Call initialization on each
        !
        do idom = 1,ndom
            call self%dom(idom)%init(mesh(idom),A%dom(idom))
        end do



        !
        ! Initialize indices of where to find off-processor recv data in the chidgVector
        !
        do idom = 1,ndom
            idomain_g = mesh(idom)%idomain_g

            do icomm = 1,size(self%dom(idom)%comm)
                proc = self%dom(idom)%comm(icomm)%proc

                do ielem = 1,size(self%dom(idom)%comm(icomm)%elem)
                    do iblk = 1,size(self%dom(idom)%comm(icomm)%elem(ielem)%blks)

                        dparent_g = self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%dparent_g()
                        eparent_g = self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%eparent_g()
                        parent_proc = self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%parent_proc()




                        if (parent_proc /= IRANK) then
                            !
                            ! For the given block, search the chidgVector recv component for the location of the off-processor data
                            !
                            do icomm_vec = 1,size(b%recv%comm)
                                do idom_vec = 1,size(b%recv%comm(icomm_vec)%dom)
                                    do ielem_vec = 1,size(b%recv%comm(icomm_vec)%dom(idom_vec)%vecs)
                                        dparent_g_vec = b%recv%comm(icomm_vec)%dom(idom_vec)%vecs(ielem_vec)%dparent_g()
                                        eparent_g_vec = b%recv%comm(icomm_vec)%dom(idom_vec)%vecs(ielem_vec)%eparent_g()

                                        if ( (dparent_g == dparent_g_vec) .and. (eparent_g == eparent_g_vec) )then
                                            self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%recv_comm    = icomm_vec
                                            self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%recv_domain  = idom_vec
                                            self%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk)%recv_element = ielem_vec

                                        end if
                                    end do !ielem_vec
                                end do !idom_vec
                            end do !icomm_vec
                        end if




                    end do !iblk
                end do !ielem
            end do !icomm
        end do !idom







    end subroutine init
    !******************************************************************************************



end module type_RASILU0_recv
