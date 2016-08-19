module mod_spatial
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, DIAG, CHIMERA, INTERIOR
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mpi_f08,                only: MPI_Barrier


    use type_chidg_data,        only: chidg_data_t
    use type_element_info,      only: element_info_t
    use type_face_info,         only: face_info_t
    use type_timer,             only: timer_t
    implicit none


contains





    !>  Spatial loop through domains, elements, and faces. Functions get called for each element/face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Improved layout, added computation of diffusion terms.
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine update_space(data,timing,info)
        implicit none
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        type(timer_t)               :: timer, comm_timer
        integer(ik)                 :: idom, ielem, iface, idiff, ifcn, ibc, ierr, nelem
        logical                     :: interior_face
        logical                     :: chimera_face 
        logical                     :: compute_face 

        logical                     :: compute_function
        logical                     :: linearize_function


        type(element_info_t)        :: elem_info
        type(face_info_t)           :: face_info




            !
            ! Start timer on spatial discretization update
            !
            call timer%start()

        

            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !
            ! Clear function_status data. This tracks if a function has already been called. So, in this way
            ! we can compute a function on a face and apply it to both elements. The function is just registered
            ! as computed for both. So we need to reset all of that data here. This is only tracked for the interior scheme.
            ! Boundary condition evaluations and Chimera faces are not tracked.
            !
            call data%sdata%function_status%clear()




            !
            ! Communicate solution vector
            !
            call comm_timer%start()
            call data%sdata%q%comm_send()
            call data%sdata%q%comm_recv()
            call data%sdata%q%comm_wait()
            call comm_timer%stop()




            !
            ! Loop through given element compute the residual functions and also the linearization of those functions
            !
            call write_line('Updating spatial scheme', io_proc=GLOBAL_MASTER)


            ! Loop through domains
            do idom = 1,data%ndomains()
                associate ( mesh => data%mesh, sdata => data%sdata, eqnset => data%eqnset(idom)%item, prop => data%eqnset(idom)%item%prop)
                nelem = mesh(idom)%nelem

                ! Loop through elements in the current domain
                do ielem = 1,nelem

                    elem_info%idomain_g  = mesh(idom)%elems(ielem)%idomain_g
                    elem_info%idomain_l  = mesh(idom)%elems(ielem)%idomain_l
                    elem_info%ielement_g = mesh(idom)%elems(ielem)%ielement_g
                    elem_info%ielement_l = mesh(idom)%elems(ielem)%ielement_l


                    ! Compute framework diffusion terms and linearization for current element
                    call sdata%update_diffusion(mesh,elem_info)


                    ! 1-6 = linearization of neighbor blocks, 7 = linearization of Q- block(self)
                    do idiff = 1,7


                        ! Faces loop. For the current element, compute the contributions from boundary integrals
                        do iface = 1,NFACES

                            face_info%idomain_g  = mesh(idom)%elems(ielem)%idomain_g
                            face_info%idomain_l  = mesh(idom)%elems(ielem)%idomain_l
                            face_info%ielement_g = mesh(idom)%elems(ielem)%ielement_g
                            face_info%ielement_l = mesh(idom)%elems(ielem)%ielement_l
                            face_info%iface      = iface


                            associate ( face => mesh(idom)%faces(ielem,iface) )

                            !
                            ! Only call the following routines for interior faces -- ftype == 0
                            ! Furthermore, only call the routines if we are computing derivatives for the neighbor of
                            ! iface or for the current element(DIAG). This saves a lot of unnecessary compute_boundary calls.
                            !
                            interior_face = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )
                            chimera_face  = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
                            compute_face  = (interior_face .or. chimera_face) .and. ( (idiff == iface) .or. (idiff == DIAG) )



                            ! Compute boundary fluxes
                            if ( compute_face ) then

                                call eqnset%compute_boundary_advective_flux(mesh, sdata, face_info, idiff)
                                call eqnset%compute_boundary_diffusive_flux(mesh, sdata, face_info, idiff)

                            end if ! compute_face


                            end associate

                        end do  ! faces loop
                        



                        !
                        ! Compute volume fluxes
                        !
                        call eqnset%compute_volume_advective_flux(mesh, sdata, elem_info, idiff)
                        call eqnset%compute_volume_diffusive_flux(mesh, sdata, elem_info, idiff)



                    end do  ! idiff
                end do  ! ielem
                end associate
            end do  ! idom






            !------------------------------------------------------------------------------------------
            !                                      Boundary Scheme
            !------------------------------------------------------------------------------------------
            !
            ! Apply boundary conditions for each domain.
            !
            do idom = 1,data%ndomains()
                call data%bcset(idom)%apply(data%mesh,data%sdata,data%eqnset(idom)%item%prop,idom)
            end do ! idom



            !
            ! Synchronize
            !
            call MPI_Barrier(ChiDG_COMM,ierr)





            call timer%stop()
            call timer%report('Spatial Discretization Time')
            call comm_timer%report('    - Spatial comm time:')
            if (present(timing)) then
                timing = timer%elapsed()
            end if


    end subroutine update_space
    !******************************************************************************************************************















end module mod_spatial
