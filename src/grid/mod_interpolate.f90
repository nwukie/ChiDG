module mod_interpolate
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: CHIMERA, INTERIOR, BOUNDARY, &
                                  LOCAL, NEIGHBOR
    use mod_DNAD_tools,     only: compute_neighbor_domain_l, compute_neighbor_element_l, compute_neighbor_face, &
                                  compute_neighbor_domain_g, compute_neighbor_element_g
    use mod_polynomial,     only: polynomialVal
    use mod_grid_tools_two, only: compute_element_donor
    use mod_chidg_mpi,      only: IRANK
    use DNAD_D

    use type_ivector,       only: ivector_t
    use type_rvector,       only: rvector_t
    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use type_seed,          only: seed_t
    use type_face_info,     only: face_info_t
    use type_chidgVector,   only: chidgVector_t
    implicit none




    interface interpolate_element
        module procedure    interpolate_element_autodiff, interpolate_element_standard
    end interface

    interface interpolate_face
        module procedure    interpolate_face_autodiff,    interpolate_face_standard
    end interface


!    interface interpolate_boundary
!        module procedure    interpolate_boundary_autodiff
!    end interface



contains


    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_element_autodiff(mesh,q,idomain_l,ielement_l,ieqn,var_gq,seed)
        type(mesh_t),        intent(in)      :: mesh(:)
        type(chidgVector_t), intent(in)      :: q
        integer(ik),         intent(in)      :: idomain_l
        integer(ik),         intent(in)      :: ielement_l
        integer(ik),         intent(in)      :: ieqn
        type(AD_D),          intent(inout)   :: var_gq(:)
        type(seed_t),        intent(in)      :: seed

        type(AD_D)  :: qdiff(mesh(idomain_l)%elems(ielement_l)%nterms_s)
        integer(ik) :: nderiv, set_deriv, iterm, igq, i, neqns_seed, nterms_s_seed
        logical     :: linearize_me


        !
        ! Get the number of degrees of freedom for the seed element
        ! and set this as the number of partial derivatives to track
        !
        if (seed%ielement_l == 0) then
            !
            ! If ielem_seed == 0 then we aren't interested in tracking derivatives
            !
            nderiv = 1

        !
        ! Get number of unknowns from element being linearized
        !
        else

            !
            ! Get number of equations and terms in solution expansions
            !
            neqns_seed    = mesh(seed%idomain_l)%elems(seed%ielement_l)%neqns
            nterms_s_seed = mesh(seed%idomain_l)%elems(seed%ielement_l)%nterms_s

            !
            ! Compute number of unknowns in the seed element, which is the number of partial derivatives we are tracking
            !
            nderiv = neqns_seed * nterms_s_seed
        end if



        !
        ! Allocate the derivative array for each autodiff variable
        ! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
        !
        do igq = 1,size(var_gq)
            var_gq(igq) = AD_D(nderiv)
        end do



        !
        ! If the current element is being differentiated (ielem == ielem_seed)
        ! then copy the solution modes to local AD variable and seed derivatives
        !
        linearize_me = ( (idomain_l == seed%idomain_l) .and. (ielement_l == seed%ielement_l) )
        !linearize_me = ( (idom_g_interp == seed%idomain_g) .and. (ielem_g_interp == seed%ielement_g) )

        if (linearize_me) then

            !
            ! Allocate derivative arrays for temporary solution variable
            !
            do iterm = 1,mesh(idomain_l)%elems(ielement_l)%nterms_s
                qdiff(iterm) = AD_D(nderiv)
            end do


            !
            ! Copy the solution variables from 'q' to 'qdiff'
            !
            qdiff = q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn)


            !
            ! Loop through the terms in qdiff
            !
            do iterm = 1,size(qdiff)
                !
                ! For the given term, seed its appropriate derivative
                !
                set_deriv = (ieqn - 1)*mesh(idomain_l)%elems(ielement_l)%nterms_s + iterm
                qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
            end do


            !
            ! Interpolate solution to GQ nodes via matrix-vector multiplication
            !
            var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%gq%vol%val,qdiff)

        else
            !
            ! If the solution variable derivatives dont need initialized
            ! then just use the q(ielem) values and derivatives get
            ! initialized to zero
            !
            var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%gq%vol%val,q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn))
        end if


    end subroutine interpolate_element_autodiff
    !***********************************************************************************************************










    !> Interpolate variable from polynomial expansion to explicit values at quadrature nodes. The automatic
    !! differentiation process really starts here, when the polynomial expansion is evaluated.
    !!
    !! The interpolation process occurs through a matrix-vector multiplication. That is an interpolation matrix
    !! multiplied by a vector of modes from the polynomial expansion. To start the automatic differentiation, 
    !! the derivative arrays of the values for the polynomial modes must be initialized before any computation. 
    !! So, before the matrix-vector multiplication.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      mesh        Array of mesh instances.
    !!  @param[in]      face        Face info, such as indices for locating the face in the mesh.
    !!  @param[in]      q           Solution vector
    !!  @param[in]      ieqn        Index of the equation variable being interpolated
    !!  @param[inout]   var_gq      Array of auto-diff values of the equation evaluated at gq points that is passed back
    !!  @param[in]      source      LOCAL/NEIGHBOR indicating which element to interpolate from
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine interpolate_face_autodiff(mesh,face,q,ieqn,var_gq,source)
        type(mesh_t),           intent(in)              :: mesh(:)
        type(face_info_t),      intent(in)              :: face
        type(chidgVector_t),    intent(in)              :: q
        integer(ik),            intent(in)              :: ieqn
        type(AD_D),             intent(inout)           :: var_gq(:)
        integer(ik),            intent(in)              :: source

        integer(ik)     :: idomain_l, ielement_l, iface
        type(seed_t)    :: seed

        type(AD_D), allocatable  :: qdiff(:)
        real(rk),   allocatable  :: interpolator(:,:)

        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s, ierr, neqns_seed, nterms_s_seed
        integer(ik) :: ndonors, idonor, donor_proc
        integer(ik) :: idom_l_interp, ielem_l_interp, iface_interp, idom_g_interp, ielem_g_interp
        integer(ik) :: recv_comm, recv_domain, recv_element
        integer(ik) :: ChiID
        logical     :: linearize_me             = .false.
        logical     :: chimera_interpolation    = .false.
        logical     :: conforming_interpolation = .false.
        logical     :: parallel_interpolation
        logical     :: parallel_seed


        ! Chimera data
        logical                     :: mask(size(var_gq))    ! node mask for Chimera quadrature points
        type(AD_D),  allocatable    :: var_gq_chimera(:)
        integer(ik), allocatable    :: gq_node_indices(:)
        integer(ik)                 :: inode


        mask = .false.



        idomain_l  = face%idomain_l
        ielement_l = face%ielement_l
        iface      = face%iface
        seed       = face%seed

        !
        ! Test if interpolating from local element
        !
        if ( source == LOCAL ) then
            conforming_interpolation = ( (mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR) .or. &
                                         (mesh(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY) .or. &
                                         (mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA) )         ! including chimera here because in the LOCAL case, it doesn't matter
            parallel_interpolation   = .false.

            ndonors = 1

        !
        ! Test if interpolating from neighbor element(s)
        !
        elseif (source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR )



            !
            ! Test for standard conforming interpolation from neighbor
            !
            if ( conforming_interpolation ) then

                ndonors = 1

            !
            ! Test for chimera interpolation from neighbor
            !
            elseif ( chimera_interpolation ) then

                ChiID   = mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                ndonors = mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()

            else
                call chidg_signal(FATAL,"interpolate_face: invalid value for 'face%ftype'")
            end if


        else
            call chidg_signal(FATAL,"interpolate_face: invalid value for incoming parameter 'source'")
        end if






        ! PROBABLY NEED TO MOVE THIS DOWN BELOW NEXT SECTION SO WE CAN DETERMINE parallel_interpolation
        ! DEPENDING ON CHIMERA DONOR.
        !
        ! Get the number of degrees of freedom for the seed element
        ! and set this as the number of partial derivatives to track
        !
        if (seed%ielement_l == 0) then
            !
            ! If ielem_seed == 0 then we aren't interested in tracking derivatives
            !
            nderiv = 1
        else

            !
            ! Get number of equations and terms in solution expansions
            !
            parallel_seed = (seed%iproc /= IRANK)
            if ( parallel_seed ) then
                neqns_seed    = q%recv%comm(seed%recv_comm)%dom(seed%recv_domain)%vecs(seed%recv_element)%nvars()
                nterms_s_seed = q%recv%comm(seed%recv_comm)%dom(seed%recv_domain)%vecs(seed%recv_element)%nterms()
            else
                neqns_seed    = mesh(seed%idomain_l)%elems(seed%ielement_l)%neqns
                nterms_s_seed = mesh(seed%idomain_l)%elems(seed%ielement_l)%nterms_s
            end if

            !
            ! Compute number of unknowns in the seed element, which is the number of partial derivatives we are tracking
            !
            nderiv = neqns_seed  *  nterms_s_seed
        end if



        !
        ! Allocate the derivative array for each autodiff variable
        ! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
        !
        do igq = 1,size(var_gq)
            allocate(var_gq(igq)%xp_ad_(nderiv))
        end do


        !
        ! For each donor element to the interpolation. (ndonors could be > 1 for Chimera interpolations)
        !
        do idonor = 1,ndonors


            !
            ! Compute neighbor access indices
            !
            if ( source == LOCAL ) then
                !
                ! Interpolate from LOCAL element
                !
                idom_l_interp  = idomain_l
                ielem_l_interp = ielement_l
                idom_g_interp  = mesh(idomain_l)%elems(ielement_l)%idomain_g
                ielem_g_interp = mesh(idomain_l)%elems(ielement_l)%ielement_g
                iface_interp   = iface

                !
                ! Get interpolation matrix from quadrature instance
                !
                interpolator = mesh(idom_l_interp)%faces(ielem_l_interp,iface_interp)%gq%face%val(:,:,iface_interp)


            elseif ( source == NEIGHBOR ) then

                !
                ! Interpolate from conforming NEIGHBOR element
                !
                if ( conforming_interpolation ) then

                    parallel_interpolation   = ( IRANK /= mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_proc )
                    if (parallel_interpolation) then
                        idom_g_interp  = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g  
                        idom_l_interp  = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
                        ielem_g_interp = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g
                        ielem_l_interp = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
                        ! THIS PROBABLY NEEDS IMPROVED
                        iface_interp   = compute_neighbor_face(     mesh,idomain_l,ielement_l,iface,idonor)

                        recv_comm    = mesh(idomain_l)%faces(ielement_l,iface)%recv_comm
                        recv_domain  = mesh(idomain_l)%faces(ielement_l,iface)%recv_domain
                        recv_element = mesh(idomain_l)%faces(ielement_l,iface)%recv_element

                        ! THIS PROBABLY NEEDS IMPROVED
                        ! Using iface_interp to get an interpolation for an opposite face. Assumes same order and eqns.
                        interpolator   = mesh(idomain_l)%faces(ielement_l,iface)%gq%face%val(:,:,iface_interp)

                    else
                        idom_g_interp  = compute_neighbor_domain_g( mesh,idomain_l,ielement_l,iface,idonor)
                        idom_l_interp  = compute_neighbor_domain_l( mesh,idomain_l,ielement_l,iface,idonor)
                        ielem_g_interp = compute_neighbor_element_g(mesh,idomain_l,ielement_l,iface,idonor)
                        ielem_l_interp = compute_neighbor_element_l(mesh,idomain_l,ielement_l,iface,idonor)
                        iface_interp   = compute_neighbor_face(     mesh,idomain_l,ielement_l,iface,idonor)

                        interpolator   = mesh(idom_l_interp)%faces(ielem_l_interp,iface_interp)%gq%face%val(:,:,iface_interp)
                    end if

                !
                ! Interpolate from CHIMERA donor element
                !
                elseif ( chimera_interpolation ) then

                    idom_g_interp   = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)
                    idom_l_interp   = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)
                    ielem_g_interp  = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_element_g%at(idonor)
                    ielem_l_interp  = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)
                    donor_proc      = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_proc%at(idonor)

                    ! Detect parallel interpolation
                    parallel_interpolation = (IRANK /= donor_proc)
                    if (parallel_interpolation) then
                         recv_comm    = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_comm%at(idonor)
                         recv_domain  = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_domain%at(idonor)
                         recv_element = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_element%at(idonor)
                    end if


                    interpolator    = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_interpolator%at(idonor)
                    gq_node_indices = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_gq_indices(idonor)%data()

                    ! Create mask over full GQ vector of only those nodes that are filled by the current element
                    do inode = 1,size(gq_node_indices)
                        mask(gq_node_indices(inode)) = .true.
                    end do

                else
                    call chidg_signal(FATAL,"interpolate_face: neighbor conforming_interpolation nor chimera_interpolation were detected")
                end if

            else
                call chidg_signal(FATAL,"interpolate_face: invalid source. LOCAL or NEIGHBOR.")
            end if




            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            !linearize_me = ( (idom_l_interp == seed%idomain_l) .and. (ielem_l_interp == seed%ielement_l) )
            linearize_me = ( (idom_g_interp == seed%idomain_g) .and. (ielem_g_interp == seed%ielement_g) )

            if ( linearize_me ) then


                !
                ! Allocate AD array to store a copy of the solution which starts the differentiation
                !
                if (parallel_interpolation) then
                    nterms_s = q%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%nterms()
                else
                    nterms_s = mesh(idom_l_interp)%elems(ielem_l_interp)%nterms_s
                end if

                if ( allocated(qdiff) ) deallocate(qdiff)
                allocate(qdiff(nterms_s), stat=ierr)
                if (ierr /= 0) call AllocationError



                !
                ! Allocate derivative arrays for temporary solution variable
                !
                do iterm = 1,nterms_s
                    allocate(qdiff(iterm)%xp_ad_(nderiv))
                end do


                !
                ! Copy the solution variables from 'q' to 'qdiff'
                !
                if (parallel_interpolation) then
                    qdiff = q%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%getvar(ieqn)
                else
                    qdiff = q%dom(idom_l_interp)%vecs(ielem_l_interp)%getvar(ieqn)
                end if



                !
                ! Loop through the terms in qdiff
                !
                do iterm = 1,size(qdiff)
                    ! For the given term, seed its appropriate derivative
                    set_deriv = (ieqn - 1)*nterms_s + iterm
                    qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
                end do



                !
                ! Interpolate solution to GQ nodes via matrix-vector multiplication
                !
                if ( conforming_interpolation ) then
                    var_gq = matmul(interpolator,  qdiff)





                elseif ( chimera_interpolation ) then

                    !
                    ! Allocate var_gq_chimera size to conform with interpolator
                    !
                    if (allocated(var_gq_chimera)) deallocate(var_gq_chimera)
                    allocate(var_gq_chimera(size(interpolator,1)))


                    !
                    ! Allocate the derivative array for each var_gq_chimera element
                    !
                    do igq = 1,size(var_gq_chimera)
                        allocate(var_gq_chimera(igq)%xp_ad_(nderiv))
                    end do


                    !
                    ! Perform interpolation
                    !
                    var_gq_chimera = matmul(interpolator,  qdiff)

                    !
                    ! Scatter chimera nodes to appropriate location in var_gq
                    !
                    var_gq = unpack(var_gq_chimera,mask,var_gq)

                else
                    call chidg_signal(FATAL,"interpolate_face: face interpolation type error")
                end if

            else

                !
                ! If the solution variable derivatives dont need initialized
                ! then just use the q(ielem) values and derivatives get
                ! initialized to zero
                !
                if ( conforming_interpolation ) then
                    if (parallel_interpolation) then
                        var_gq = matmul(interpolator,  q%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%getvar(ieqn))
                    else
                        var_gq = matmul(interpolator,  q%dom(idom_l_interp)%vecs(ielem_l_interp)%getvar(ieqn))
                    end if

                elseif ( chimera_interpolation ) then

                    !
                    ! Allocate var_gq_chimera size to conform with interpolator
                    !
                    if (allocated(var_gq_chimera)) deallocate(var_gq_chimera)
                    allocate(var_gq_chimera(size(interpolator,1)))


                    !
                    ! Allocate the derivative array for each var_gq_chimera element
                    !
                    do igq = 1,size(var_gq_chimera)
                        allocate(var_gq_chimera(igq)%xp_ad_(nderiv))
                    end do


                    !
                    ! Perform interpolation
                    !
                    if (parallel_interpolation) then
                        var_gq_chimera = matmul(interpolator,  q%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%getvar(ieqn))
                    else
                        var_gq_chimera = matmul(interpolator,  q%dom(idom_l_interp)%vecs(ielem_l_interp)%getvar(ieqn))
                    end if

                    !
                    ! Scatter chimera nodes to appropriate location in var_gq
                    !
                    var_gq = unpack(var_gq_chimera,mask,var_gq)

                else
                    call chidg_signal(FATAL,"interpolate_face: face interpolation type error")
                end if

            end if

            ! Reset Chimera mask
            mask = .false.

        end do ! idonor


    end subroutine interpolate_face_autodiff
    !*************************************************************************************************************













    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine interpolate_element_standard(mesh,q,idomain_l,ielement_l,ieqn,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l
        integer(ik),            intent(in)      :: ielement_l
        integer(ik),            intent(in)      :: ieqn
        real(rk),               intent(inout)   :: var_gq(:)

        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the quadrature matrix
        ! with the array of modes for the given variable.
        !
        var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%gq%vol%val, q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn))

    end subroutine interpolate_element_standard
    !*************************************************************************************************************













    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine interpolate_face_standard(mesh,q,idomain_l,ielement_l,iface,ieqn,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l, ielement_l, iface, ieqn
        real(rk),               intent(inout)   :: var_gq(:)


        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the face quadrature matrix
        ! with the array of modes for the given variable
        !
        var_gq = matmul(mesh(idomain_l)%faces(ielement_l,iface)%gq%face%val(:,:,iface), q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn))



    end subroutine interpolate_face_standard
    !*************************************************************************************************************












!
!    !> Interpolate variable from polynomial expansion to explicit values at quadrature nodes. The automatic
!    !! differentiation process really starts here, when the polynomial expansion is evaluated.
!    !!
!    !! The interpolation process occurs through a matrix-vector multiplication. That is an interpolation matrix
!    !! multiplied by a vector of modes from the polynomial expansion. To start the automatic differentiation, 
!    !! the derivative arrays of the values for the polynomial modes must be initialized before any computation. 
!    !! So, before the matrix-vector multiplication.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/1/2016
!    !!
!    !!  @param[in]      mesh        Array of mesh instances.
!    !!  @param[in]      face        Face info, such as indices for locating the face in the mesh.
!    !!  @param[in]      q           Solution vector
!    !!  @param[in]      ieqn        Index of the equation variable being interpolated
!    !!  @param[inout]   var_gq      Array of auto-diff values of the equation evaluated at gq points that is passed back
!    !!  @param[in]      source      LOCAL/NEIGHBOR indicating which element to interpolate from
!    !!
!    !-----------------------------------------------------------------------------------------------------------
!    subroutine interpolate_boundary_autodiff(mesh,face,q,ieqn,points,var)
!        type(mesh_t),           intent(in)              :: mesh(:)
!        type(face_info_t),      intent(in)              :: face
!        type(chidgVector_t),    intent(in)              :: q
!        integer(ik),            intent(in)              :: ieqn
!        type(point_t),          intent(in)              :: points(:)
!        type(AD_D),             intent(inout)           :: var(:)
!
!        integer(ik)     :: idomain_l, ielement_l, iface
!        type(seed_t)    :: seed
!
!        type(AD_D), allocatable  :: qdiff(:)
!        type(AD_D), allocatable  :: tmp(:)
!        real(rk),   allocatable  :: interpolator(:,:)
!
!        real(rk)        :: xi, eta, zeta
!        integer(ik)     :: nderiv, set_deriv, iterm, nterms_s, ierr, neqns_seed, nterms_s_seed, ipnt, nterms
!        integer(ik)     :: idomain_l_seed, ielement_l_seed, idonor, idomain_l_interp, ielement_l_interp, igq
!        integer(ik)     :: idom_d, ielem_d
!        type(point_t)   :: point_comp, node
!        type(ivector_t) :: donor_domains, donor_elements
!        type(rvector_t) :: donor_xi, donor_eta, donor_zeta
!        logical         :: linearize_me             = .false.
!
!
!
!
!
!        idomain_l  = face%idomain_l
!        ielement_l = face%ielement_l
!        iface      = face%iface
!        seed       = face%seed
!
!
!        !
!        ! Get domain/element index that is being differentiated
!        !
!        idomain_l_seed  = seed%idomain_l
!        ielement_l_seed = seed%ielement_l
!
!
!
!        !
!        ! Get the number of degrees of freedom for the seed element
!        ! and set this as the number of partial derivatives to track
!        !
!        if (ielement_l_seed == 0) then
!            !
!            ! If ielem_seed == 0 then we aren't interested in tracking derivatives
!            !
!            nderiv = 1
!
!        else
!            !
!            ! Get number of equations and terms in solution expansions
!            !
!            neqns_seed    = mesh(idomain_l_seed)%elems(ielement_l_seed)%neqns
!            nterms_s_seed = mesh(idomain_l_seed)%elems(ielement_l_seed)%nterms_s
!
!            !
!            ! Compute number of unknowns in the seed element, which is the number of partial derivatives we are tracking
!            !
!            nderiv = neqns_seed  *  nterms_s_seed
!        end if
!
!
!
!        !
!        ! Allocate the derivative array for each autodiff variable
!        ! MIGHT NOT NEED THIS IF IT GETS AUTOMATICALLY ALLOCATED ON ASSIGNMENT -- TEST
!        !
!        do igq = 1,size(var)
!            allocate(var(igq)%xp_ad_(nderiv))
!        end do
!
!
!
!
!        !
!        ! Find elements associated with boundary nodes
!        !
!        do ipnt = 1,size(points)
!
!
!            call compute_element_donor(mesh, points(ipnt), idom_d, ielem_d, point_comp)
!
!
!            call donor_domains%push_back(  idom_d  )
!            call donor_elements%push_back( ielem_d )
!            call donor_xi%push_back(   point_comp%c1_ )
!            call donor_eta%push_back(  point_comp%c2_ )
!            call donor_zeta%push_back( point_comp%c3_ )
!
!        end do ! ipnt
!
!
!
!
!        !
!        ! For each donor element to the interpolation.
!        !
!        do idonor = 1,donor_elements%size()
!
!
!            !
!            ! Get donor element and coordinate
!            !
!            idom_interp  = donor_domains%at(idonor)
!            ielem_interp = donor_elements%at(idonor)
!
!            xi   = donor_xi%at(  idonor)
!            eta  = donor_eta%at( idonor)
!            zeta = donor_zeta%at(idonor)
!
!            call node%set(xi,eta,zeta)
!
!
!            !
!            ! Get interpolation matrix from quadrature instance
!            !
!            nterms = mesh(idom_interp)%elems(ielem_interp)%nterms_s
!            if (allocated(interpolator)) deallocate(interpolator)
!            allocate(interpolator(1,nterms), stat=ierr)
!            if (ierr /= 0) call AllocationError
!
!            do iterm = 1,nterms
!                interpolator(1,iterm) = polynomialVal(3,nterms,iterm,node)
!            end do ! imode
!
!
!
!
!
!            !
!            ! If the current element is being differentiated (ielem == ielem_seed)
!            ! then copy the solution modes to local AD variable and seed derivatives
!            !
!
!            linearize_me = ( (idom_interp == idom_seed) .and. (ielem_interp == ielem_seed) )
!
!            if ( linearize_me ) then
!
!
!                !
!                ! Allocate AD array to store a copy of the solution which starts the differentiation
!                !
!                nterms_s = mesh(idom_interp)%elems(ielem_interp)%nterms_s
!                if ( allocated(qdiff) ) deallocate(qdiff)
!                allocate(qdiff(nterms_s), stat=ierr)
!                if (ierr /= 0) call AllocationError
!
!
!                !
!                ! Allocate derivative arrays for temporary solution variable
!                !
!                do iterm = 1,nterms_s
!                    allocate(qdiff(iterm)%xp_ad_(nderiv))
!                end do
!
!
!                !
!                ! Copy the solution variables from 'q' to 'qdiff'
!                !
!                qdiff = q%dom(idom_interp)%vecs(ielem_interp)%getvar(ieqn)
!
!
!                !
!                ! Loop through the terms in qdiff and seed derivatives
!                !
!                do iterm = 1,size(qdiff)
!                    ! For the given term, seed its appropriate derivative
!                    set_deriv = (ieqn - 1)*nterms_s + iterm
!                    qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
!                end do
!
!
!
!
!                !
!                ! Interpolate solution to GQ nodes via matrix-vector multiplication
!                !
!                !var(idonor) = matmul(interpolator(1,:),  qdiff)
!                tmp = matmul(interpolator,  qdiff)
!                var(idonor) = tmp(1)
!
!            else
!
!                !
!                ! If the solution variable derivatives dont need initialized
!                ! then just use the q(ielem) values and derivatives get
!                ! initialized to zero
!                !
!
!                if ( allocated(tmp) ) deallocate(tmp)
!                allocate(tmp(1))
!                allocate(tmp(1)%xp_ad_(nderiv))
!
!                
!                !var(idonor) = matmul(interpolator(1,:),  q%dom(idom_interp)%vecs(ielem_interp)%getvar(ieqn))
!                tmp = matmul(interpolator,  q%dom(idom_interp)%vecs(ielem_interp)%getvar(ieqn))
!                var(idonor) = tmp(1)
!
!            end if
!
!
!
!
!        end do ! idonor
!
!
!    end subroutine interpolate_boundary_autodiff
!    !*************************************************************************************************************















end module mod_interpolate







