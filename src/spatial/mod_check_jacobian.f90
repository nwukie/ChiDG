module mod_check_jacobian
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES, DIAG, ZERO, BOUNDARY_ADVECTIVE_FLUX, VOLUME_ADVECTIVE_FLUX

    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_chidg_data,    only: chidg_data_t
    use type_chidg_worker,  only: chidg_worker_t

    use type_blockvector
    use type_densematrix,   only: densematrix_t

    use mod_DNAD_tools,     only: element_compute_seed, face_compute_seed
    implicit none


contains





    !> Compute both the Finite Difference representation and DNAD representation of
    !! the volume advective flux jacobian.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[inout]   domain      Domain structure containing grid, eqnset, etc
    !!  @param[in]      ielem       Element index for computing the jacobian
    !!  @param[in]      iblk        Block index of the linearization
    !!  @param[inout]   blk_dnad    Densematrix for storing the DNAD jacobian
    !!  @param[inout]   blk_fd      Densematrix for storing the Finite Difference jacobian
    !!
    !----------------------------------------------------------------------------------
    subroutine check_jacobian_volume_advective_flux(data,ielem,iblk,blk_dnad,blk_fd)
        type(chidg_data_t),  intent(inout)       :: data
        integer(ik),         intent(in)          :: ielem, iblk
        type(densematrix_t), intent(inout)       :: blk_dnad, blk_fd

        type(blockvector_t), allocatable    :: rhs_r, vec_fd
        real(rk)                            :: qhold, eps
        integer(ik)                         :: nelem, i, iterm, icol, nterms, ivar, iflux, nflux, idom, iel, idepend
        integer(ik)                         :: ielem_p              ! element in which the solution is being perturbed for the FD calculation.
        type(chidg_worker_t)                :: worker


        ! ASSUME ONLY ONE DOMAIN IS PRESENT
        idom = 1

        idepend = 1



        associate ( mesh => data%mesh, eqnset => data%eqnset, sdata => data%sdata, q => data%sdata%q%dom(1), rhs => data%sdata%rhs%dom(1), lhs => data%sdata%lhs%dom(1), prop => data%eqnset(1)%prop )


            nelem = mesh(1)%nelem
            nterms = mesh(1)%nterms_s
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            call worker%init(data%mesh,data%sdata)

            !
            ! Clear all working data
            !
            call sdata%lhs%clear()
            call sdata%rhs%clear()
            call sdata%function_status%clear()

            rhs_r  = sdata%rhs%dom(1)      ! Alloate supporting storage
            vec_fd = sdata%rhs%dom(1)



            !
            ! For the current element, compute the contributions from volume integrals
            !
            idom = 1

            worker%element_info%idomain_g  = 1
            worker%element_info%idomain_l  = 1
            worker%element_info%ielement_g = ielem
            worker%element_info%ielement_l = ielem



            if (allocated(eqnset(1)%volume_advective_operator)) then
                nflux = size(eqnset(1)%volume_advective_operator)
                do iflux = 1,nflux

                    worker%function_info%type    = VOLUME_ADVECTIVE_FLUX
                    worker%function_info%ifcn    = iflux
                    worker%function_info%idiff   = iblk
                    worker%function_info%idepend = 1
                    worker%function_info%seed    = element_compute_seed(mesh,idom,ielem,idepend,iblk)

                    !call data%eqnset(1)%volume_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,elem_info,function_info)
                    call data%eqnset(1)%volume_advective_operator(iflux)%op%compute(worker,prop)
                end do

            else
                call chidg_signal(WARN,"No volume advective flux was found")
            end if
            

           !
           ! Store linearization from DNAD. This is what we want to check against the FD calculation.
           !
           blk_dnad = lhs%lblks(ielem,iblk)
           blk_fd = blk_dnad                       ! Sourced allocation
           blk_fd%mat = ZERO                       ! Zero storage


           !
           ! Store temporary rhs
           !
           rhs_r%vecs(ielem) = rhs%vecs(ielem)


           !
           ! Reset sdata storage
           !
           call lhs%clear()
           call rhs%clear()
           call sdata%function_status%clear()


           !
           ! Select in which element, the solution is being perturbed
           !
           if (iblk == DIAG) then
               ielem_p = ielem
           else
               ielem_p = mesh(1)%faces(ielem,iblk)%get_neighbor_element_l()
           end if




           eps   = 1.e-8_rk ! finite-difference perturbation

           !
           ! Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
           !
           do ivar = 1,eqnset(1)%prop%nequations()
               do iterm = 1,mesh(1)%nterms_s

                   !
                   ! Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                   !
                   qhold = q%vecs(ielem_p)%getterm(ivar,iterm)
                   call q%vecs(ielem_p)%setterm(ivar,iterm,qhold + eps)


                   !
                   ! For the current element, compute the contributions from volume integrals
                   !
                   if (allocated(eqnset(1)%volume_advective_operator)) then
                        nflux = size(eqnset(1)%volume_advective_operator)
                        do iflux = 1,nflux
                            worker%function_info%type    = VOLUME_ADVECTIVE_FLUX
                            worker%function_info%ifcn    = iflux
                            worker%function_info%idiff   = iblk
                            worker%function_info%idepend = 1
                            worker%function_info%seed    = element_compute_seed(mesh,idom,ielem,idepend,iblk)

                            !call eqnset(1)%volume_advective_flux(iflux)%flux%compute(mesh,sdata,prop,elem_info,function_info)
                            call eqnset(1)%volume_advective_operator(iflux)%op%compute(worker,prop)
                        end do

                   else
                       call chidg_signal(WARN,"No volume advective flux was found")
                   end if


                   !
                   ! Return perturbed value to normal state
                   !
                   call q%vecs(ielem_p)%setterm(ivar,iterm,qhold)


                   !
                   ! Compute finite difference jacobian
                   !
                   vec_fd = (rhs - rhs_r)/eps


                   !
                   ! Store to column of blk_fd
                   !
                   icol = (ivar-1)*nterms + iterm                  ! Compute appropriate column for storing linearization
                   blk_fd%mat(:,icol) = vec_fd%vecs(ielem)%vec    ! Store finite difference linearization of the residual


                   !
                   ! Reset sdata storage
                   !
                   call sdata%lhs%clear()
                   call sdata%rhs%clear()
                   call sdata%function_status%clear()


               end do   ! iterm
           end do   ! ivar



        end associate

    end subroutine check_jacobian_volume_advective_flux
    !********************************************************************************************








    !> compute both the finite difference representation and dnad representation of
    !! the boundary advective flux jacobian.
    !!
    !!   @author Nathan A. Wukie
    !!   @date   1/28/2016
    !!
    !!   @param[inout]   domain      domain structure containing grid, eqnset, etc
    !!   @param[in]      ielem       element index for computing the jacobian
    !!   @param[in]      iblk        block index of the linearization
    !!   @param[inout]   blk_dnad    densematrix for storing the dnad jacobian
    !!   @param[inout]   blk_fd      densematrix for storing the finite difference jacobian
    !!
    !--------------------------------------------------------------------------------------------
    subroutine check_jacobian_boundary_advective_flux(data,ielem,iblk,blk_dnad,blk_fd)
        type(chidg_data_t),     intent(inout)   :: data
        integer(ik),            intent(in)      :: ielem, iblk
        type(densematrix_t),    intent(inout)   :: blk_dnad, blk_fd

        type(blockvector_t)  :: rhs_r, vec_fd
        real(rk)    :: qhold, eps
        integer(ik) :: nelem, i, iterm, iface, ivar, icol, nterms, iflux, nflux, idom, idonor
        integer(ik) :: ielem_p              ! ielem_p is the element in which the solution is being perturbed for the finite difference calculation.

        type(chidg_worker_t)    :: worker


        call worker%init(data%mesh,data%sdata)


        !
        ! ASSUMING THERE IS ONLY ONE DOMAIN
        !
        idom   = 1
        idonor = 1

        !
        ! Select in which element, the solution is being perturbed
        ! Select which face the boundary integral is being computed on based on the direction of iblk.
        !
        if (iblk == DIAG) then
            ielem_p = ielem
            iface   = 1
        else
            ielem_p = data%mesh(1)%faces(ielem,iblk)%get_neighbor_element_l()
            iface   = iblk
        end if



        associate ( mesh => data%mesh, q => data%sdata%q%dom(1), rhs => data%sdata%rhs%dom(1), lhs => data%sdata%lhs%dom(1), prop => data%eqnset(1)%prop )

            nelem = mesh(1)%nelem
            nterms = mesh(1)%nterms_s
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !
            ! Zero data storage
            !
            call lhs%clear()      
            call lhs%clear()
            call data%sdata%function_status%clear()


            rhs_r  = rhs          ! Allocate supporing storage containers
            vec_fd = rhs


            worker%face_info%idomain_g  = data%mesh(idom)%elems(ielem)%idomain_g
            worker%face_info%idomain_l  = data%mesh(idom)%elems(ielem)%idomain_l
            worker%face_info%ielement_g = data%mesh(idom)%elems(ielem)%ielement_g
            worker%face_info%ielement_l = data%mesh(idom)%elems(ielem)%ielement_l
            worker%face_info%iface      = iface


            !
            ! For the current element, compute the contributions from boundary integrals
            !
            if (allocated(data%eqnset(1)%boundary_advective_operator)) then
                nflux = size(data%eqnset(1)%boundary_advective_operator)
                do iflux = 1,nflux

                    worker%function_info%type    = BOUNDARY_ADVECTIVE_FLUX
                    worker%function_info%ifcn    = iflux
                    worker%function_info%idiff   = iblk
                    worker%function_info%idepend = idonor
                    worker%function_info%seed    = face_compute_seed(data%mesh,idom,ielem,iface,idonor,iblk)


                    !call data%eqnset(1)%boundary_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,face_info, function_info)
                    call data%eqnset(1)%boundary_advective_operator(iflux)%op%compute(worker,prop)
                end do

            else
                call chidg_signal(WARN,"No boundary advective flux was found")
            end if


            !
            ! Store linearization computed from DNAD
            !
            blk_dnad   = lhs%lblks(ielem,iblk)
            blk_fd     = blk_dnad   ! Sourced allocation
            blk_fd%mat = ZERO       ! Zero storage



            ! Reset sdata storage
            call lhs%clear()
            call rhs%clear()
            call data%sdata%function_status%clear()



            !
            ! For the current element, compute the contributions from boundary integrals
            !
            ! Need to use DIAG to get rhs for finite difference calculation. 
            ! This is because RHS is only stored for DIAG in the integrate procedure.
            if (allocated(data%eqnset(1)%boundary_advective_operator)) then
                nflux = size(data%eqnset(1)%boundary_advective_operator)
                do iflux = 1,nflux

                    worker%function_info%type    = BOUNDARY_ADVECTIVE_FLUX
                    worker%function_info%ifcn    = iflux
                    worker%function_info%idiff   = DIAG
                    worker%function_info%idepend = idonor

                    !call data%eqnset(1)%boundary_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,face_info,function_info)
                    call data%eqnset(1)%boundary_advective_operator(iflux)%op%compute(worker,prop)

                end do

            else
                call chidg_signal(WARN,"No boundary advective flux was found")
            end if


            !
            ! Store temporary rhs
            !
            rhs_r%vecs(ielem) = rhs%vecs(ielem)


            !
            ! Reset sdata storage
            !
            call lhs%clear()
            call rhs%clear()
            call data%sdata%function_status%clear()




            eps   = 1.e-8_rk    ! finite-difference perturbation

            !
            ! Loop through terms, perturb term, compute rhs, compute finite difference jacobian, return term.
            !
            do ivar = 1,data%eqnset(1)%prop%nequations()
                do iterm = 1,data%mesh(1)%nterms_s

                    !
                    ! Perturb the iterm-th term in the solution expansion for variable ivar in element ielem.
                    !
                    qhold = q%vecs(ielem_p)%getterm(ivar,iterm)
                    call q%vecs(ielem_p)%setterm(ivar,iterm,qhold + eps)


                    !
                    ! For the current element, compute the contributions from volume integrals
                    !
                    ! Need to use DIAG to get rhs for finite difference calculation. 
                    ! This is because RHS is only stored for DIAG in the integrate procedure.
                    if (allocated(data%eqnset(1)%boundary_advective_operator)) then
                        nflux = size(data%eqnset(1)%boundary_advective_operator)
                        do iflux = 1,nflux

                            worker%function_info%type    = BOUNDARY_ADVECTIVE_FLUX
                            worker%function_info%ifcn    = iflux
                            worker%function_info%idiff   = DIAG
                            worker%function_info%idepend = idonor

                            !call data%eqnset(1)%boundary_advective_flux(iflux)%flux%compute(data%mesh,data%sdata,prop,face_info,function_info)
                            call data%eqnset(1)%boundary_advective_operator(iflux)%op%compute(worker,prop)
                        end do

                    else
                        print*, data%eqnset(1)%name
                        call chidg_signal(WARN,"No boundary advective flux was found")
                    end if






                    !
                    ! Return perturbed value to normal state
                    !
                    call q%vecs(ielem_p)%setterm(ivar,iterm,qhold)


                    !
                    ! Compute finite difference jacobian
                    !
                    vec_fd = (rhs - rhs_r)/eps


                    !
                    ! Store to column of blk_fd
                    !
                    icol = (ivar-1)*nterms + iterm                  ! Compute appropriate column for storing linearization
                    blk_fd%mat(:,icol) = vec_fd%vecs(ielem)%vec


                    !
                    ! Reset sdata storage
                    !
                    call lhs%clear()
                    call rhs%clear()
                    call data%sdata%function_status%clear()

                end do  ! iterm
            end do  ! ivar


        end associate

    end subroutine check_jacobian_boundary_advective_flux
    !***************************************************************************************








end module mod_check_jacobian
