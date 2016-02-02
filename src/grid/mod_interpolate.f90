module mod_interpolate
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: CHIMERA, INTERIOR, BOUNDARY, &
                                  LOCAL, NEIGHBOR
    use DNAD_D
    use type_mesh,          only: mesh_t
    use type_seed,          only: seed_t
    use type_face_info,     only: face_info_t
    use type_chidgVector,   only: chidgVector_t
    use mod_DNAD_tools,     only: compute_neighbor_domain, compute_neighbor_element, compute_neighbor_face
    implicit none




    interface interpolate_element
        module procedure    interpolate_element_autodiff, interpolate_element_standard
    end interface

    interface interpolate_face
        module procedure    interpolate_face_autodiff,    interpolate_face_standard
    end interface





contains


    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------
    subroutine interpolate_element_autodiff(mesh,q,idom,ielem,ieqn,var_gq,seed)
        type(mesh_t),        intent(in)      :: mesh(:)
        type(chidgVector_t), intent(in)      :: q
        integer(ik),         intent(in)      :: idom
        integer(ik),         intent(in)      :: ielem
        integer(ik),         intent(in)      :: ieqn
        type(AD_D),          intent(inout)   :: var_gq(:)
        type(seed_t),        intent(in)      :: seed

        type(AD_D)  :: qdiff(mesh(idom)%elems(ielem)%nterms_s)
        integer(ik) :: nderiv, set_deriv, iterm, igq, i, neqns_seed, nterms_s_seed
        integer(ik) :: idom_seed, ielem_seed
        logical     :: linearize_me


        !
        ! Get domain/element index that is being differentiated
        !
        idom_seed  = seed%idom
        ielem_seed = seed%ielem



        !
        ! Get the number of degrees of freedom for the seed element
        ! and set this as the number of partial derivatives to track
        !
        if (ielem_seed == 0) then
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
            neqns_seed    = mesh(idom_seed)%elems(ielem_seed)%neqns
            nterms_s_seed = mesh(idom_seed)%elems(ielem_seed)%nterms_s

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
        linearize_me = ( (idom == idom_seed) .and. (ielem == ielem_seed) )

        if (linearize_me) then

            !
            ! Allocate derivative arrays for temporary solution variable
            !
            do iterm = 1,mesh(idom)%elems(ielem)%nterms_s
                qdiff(iterm) = AD_D(nderiv)
            end do


            !
            ! Copy the solution variables from 'q' to 'qdiff'
            !
            qdiff = q%dom(idom)%lvecs(ielem)%getvar(ieqn)


            !
            ! Loop through the terms in qdiff
            !
            do iterm = 1,size(qdiff)
                !
                ! For the given term, seed its appropriate derivative
                !
                set_deriv = (ieqn - 1)*mesh(idom)%elems(ielem)%nterms_s + iterm
                qdiff(iterm)%xp_ad_(set_deriv) = 1.0_rk
            end do


            !
            ! Interpolate solution to GQ nodes via matrix-vector multiplication
            !
            var_gq = matmul(mesh(idom)%elems(ielem)%gq%vol%val,qdiff)

        else
            !
            ! If the solution variable derivatives dont need initialized
            ! then just use the q(ielem) values and derivatives get
            ! initialized to zero
            !
            var_gq = matmul(mesh(idom)%elems(ielem)%gq%vol%val,q%dom(idom)%lvecs(ielem)%getvar(ieqn))
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

        integer(ik)     :: idom, ielem, iface
        type(seed_t)    :: seed

        type(AD_D), allocatable  :: qdiff(:)
        real(rk),   allocatable  :: interpolator(:,:)

        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s, ierr, neqns_seed, nterms_s_seed
        integer(ik) :: idom_seed, ielem_seed, ndonors, idonor, idom_interp, ielem_interp, iface_interp
        integer(ik) :: ChiID
        logical     :: linearize_me             = .false.
        logical     :: chimera_interpolation    = .false.
        logical     :: conforming_interpolation = .false.


        ! Chimera data
        logical                     :: mask(size(var_gq))    ! node mask for Chimera quadrature points
        type(AD_D),  allocatable    :: var_gq_chimera(:)
        integer(ik), allocatable    :: gq_node_indices(:)
        integer(ik)                 :: inode


        mask = .false.

        idom  = face%idomain
        ielem = face%ielement
        iface = face%iface
        seed  = face%seed

        !
        ! Test if interpolating from local element
        !
        if ( source == LOCAL ) then
            conforming_interpolation = ( (mesh(idom)%faces(ielem,iface)%ftype == INTERIOR) .or. &
                                         (mesh(idom)%faces(ielem,iface)%ftype == BOUNDARY) .or. &
                                         (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA) )         ! including chimera here because in the LOCAL case, it doesn't matter

            ndonors = 1

        !
        ! Test if interpolating from neighbor element(s)
        !
        elseif (source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )


            !
            ! Test for standard conforming interpolation from neighbor
            !
            if ( conforming_interpolation ) then

                ndonors = 1

            !
            ! Test for chimera interpolation from neighbor
            !
            elseif ( chimera_interpolation ) then

                ChiID   = mesh(idom)%faces(ielem,iface)%ChiID
                ndonors = mesh(idom)%chimera%recv%data(ChiID)%ndonors


            else

                call chidg_signal(FATAL,"interpolate_face: invalid value for 'face%ftype'")

            end if


        else

            call chidg_signal(FATAL,"interpolate_face: invalid value for incoming parameter 'source'")

        end if




        !
        ! Get domain/element index that is being differentiated
        !
        idom_seed  = seed%idom
        ielem_seed = seed%ielem


        !
        ! Get the number of degrees of freedom for the seed element
        ! and set this as the number of partial derivatives to track
        !
        if (ielem_seed == 0) then
            !
            ! If ielem_seed == 0 then we aren't interested in tracking derivatives
            !
            nderiv = 1
        else
            !
            ! Get number of equations and terms in solution expansions
            !
            neqns_seed    = mesh(idom_seed)%elems(ielem_seed)%neqns
            nterms_s_seed = mesh(idom_seed)%elems(ielem_seed)%nterms_s

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
                idom_interp  = idom
                ielem_interp = ielem
                iface_interp = iface

                !
                ! Get interpolation matrix from quadrature instance
                !
                interpolator = mesh(idom_interp)%faces(ielem_interp,iface_interp)%gq%face%val(:,:,iface_interp)


            elseif ( source == NEIGHBOR ) then

                !
                ! Interpolate from conforming NEIGHBOR element
                !
                if ( conforming_interpolation ) then

                    idom_interp  = compute_neighbor_domain( mesh,idom,ielem,iface,idonor)
                    ielem_interp = compute_neighbor_element(mesh,idom,ielem,iface,idonor)
                    iface_interp = compute_neighbor_face(   mesh,idom,ielem,iface,idonor)

                    interpolator = mesh(idom_interp)%faces(ielem_interp,iface_interp)%gq%face%val(:,:,iface_interp)

                !
                ! Interpolate from CHIMERA NEIGHBOR element
                !
                elseif ( chimera_interpolation ) then

                    idom_interp  = mesh(idom)%chimera%recv%data(ChiID)%donor_domain%at(idonor)
                    ielem_interp = mesh(idom)%chimera%recv%data(ChiID)%donor_element%at(idonor)

                    interpolator = mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator%at(idonor)
                    gq_node_indices = mesh(idom)%chimera%recv%data(ChiID)%donor_gq_indices(idonor)%data()

                    ! Create mask over full GQ vector of only those nodes that are filled by the current element
                    do inode = 1,size(gq_node_indices)
                        mask(gq_node_indices(inode)) = .true.
                    end do



                !
                ! Error case
                !
                else
                    call chidg_signal(FATAL,"interpolate_face: neighbor conforming_interpolation nor chimera_interpolation were detected")
                end if
            end if



            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            linearize_me = ( (idom_interp == idom_seed) .and. (ielem_interp == ielem_seed) )

            if ( linearize_me ) then


                !
                ! Allocate AD array to store a copy of the solution which starts the differentiation
                !
                nterms_s = mesh(idom_interp)%elems(ielem_interp)%nterms_s
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
                qdiff = q%dom(idom_interp)%lvecs(ielem_interp)%getvar(ieqn)


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
                    var_gq = matmul(interpolator,  q%dom(idom_interp)%lvecs(ielem_interp)%getvar(ieqn))
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
                    var_gq_chimera = matmul(interpolator,  q%dom(idom_interp)%lvecs(ielem_interp)%getvar(ieqn))

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
    subroutine interpolate_element_standard(mesh,q,idom,ielem,ieqn,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idom
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: ieqn
        real(rk),               intent(inout)   :: var_gq(:)

        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the quadrature matrix
        ! with the array of modes for the given variable.
        !
        var_gq = matmul(mesh(idom)%elems(ielem)%gq%vol%val, q%dom(idom)%lvecs(ielem)%getvar(ieqn))

    end subroutine interpolate_element_standard
    !*************************************************************************************************************













    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine interpolate_face_standard(mesh,q,idom,ielem,iface,ieqn,var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idom, ielem, iface, ieqn
        real(rk),               intent(inout)   :: var_gq(:)


        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the face quadrature matrix
        ! with the array of modes for the given variable
        !
        var_gq = matmul(mesh(idom)%faces(ielem,iface)%gq%face%val(:,:,iface), q%dom(idom)%lvecs(ielem)%getvar(ieqn))



    end subroutine interpolate_face_standard
    !*************************************************************************************************************
















end module mod_interpolate







