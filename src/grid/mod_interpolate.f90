!>  ChiDG Interpolation Utilities
!!
!!  Purpose: 
!!  ---------------------------------------
!!  Facilitate interpolation from modal representations 
!!  to node discrete representations of variables.
!!
!!
!!  Public routines:
!!  ---------------------------------------
!!      interpolate_element_autodiff
!!      interpolate_face_autodiff 
!!
!!      interpolate_element_standard
!!      interpolate_face_standard
!!
!!
!!  Module utility routines:
!!  ---------------------------------------
!!      get_face_interpolation_info
!!      get_face_interpolation_interpolator
!!      get_face_interpolation_mask
!!      get_face_interpolation_comm
!!      get_face_interpolation_ndonors
!!      get_face_interpolation_style
!!      get_interpolation_nderiv
!!
!!
!---------------------------------------------------------------------------------------------
module mod_interpolate
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: CHIMERA, INTERIOR, BOUNDARY, &
                                  ME, NEIGHBOR, ONE, ZERO
                                  
    use mod_polynomial,     only: polynomialVal
    use mod_chidg_mpi,      only: IRANK
    use mod_DNAD_tools,     only: compute_neighbor_face
    use DNAD_D

    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use type_recv,          only: recv_t
    use type_chidg_vector,   only: chidg_vector_t
    implicit none



contains




    !>  Interpolate polynomial expansion to volume quadrature node set, initializing the
    !!  automatic differentiation process if needed.
    !!
    !!  This returns an array of automatic differentiation(AD) values at the quadrature nodes,
    !!  that also have their derivatives initialized.
    !!
    !!  Some interpolation parameters to note that are used inside here:
    !!      - interpolation_type:   'value', 'ddx', 'ddy', 'ddz'
    !!
    !!  Partial mode expansions can be used to perform interpolations. To do this, pass the 
    !!  Pmin and Pmax optional arguments. Pmin is the minimum 1D polynomial order to be used
    !!  in the interpolation. Pmax is the maximum 1D polynomial order to be used in the 
    !!  interpolation. One might pass Pmin=0, Pmax=1. This would return an interpolation that
    !!  used only the first 8 modes in the polynomial expansion. The constant mode and all
    !!  the linear modes, but no modes of higher order.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!  @param[in]  itime - Index for the time step in solution
    !!
    !-----------------------------------------------------------------------------------------
    function interpolate_element_autodiff(mesh,q,elem_info,fcn_info,ieqn,itime,interpolation_type,Pmin,Pmax) result(var_gq)
        type(mesh_t),           intent(in)              :: mesh(:)
        type(chidg_vector_t),   intent(in)              :: q
        type(element_info_t),   intent(in)              :: elem_info
        type(function_info_t),  intent(in)              :: fcn_info
        integer(ik),            intent(in)              :: ieqn
        integer(ik),            intent(in)              :: itime
        character(*),           intent(in)              :: interpolation_type
        integer(ik),            intent(in), optional    :: Pmin
        integer(ik),            intent(in), optional    :: Pmax


        character(:),   allocatable :: user_msg

        type(AD_D)              :: var_gq(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)
        type(AD_D)              :: qdiff(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s)
        real(rk),   allocatable :: qtmp(:)

        integer(ik) :: nderiv, set_deriv, iterm, mode_min, mode_max
        logical     :: differentiate_me

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l )

        !
        ! Get number of derivatives to initialize for automatic differentiation
        !
        nderiv = get_interpolation_nderiv(mesh,fcn_info)


        !
        ! Allocate derivative arrays for temporary solution variable
        !
        do iterm = 1,mesh(idom)%elems(ielem)%nterms_s
            qdiff(iterm) = AD_D(nderiv)
        end do


        !
        ! Copy the solution variables from 'q' to 'qtmp'
        !
        qtmp = q%dom(idom)%vecs(ielem)%getvar(ieqn,itime)

        !
        ! If Pmin,Pmax are present, use a subset of the expansion modes.
        !
        ! If not present, use the entire expansion.
        !
        qdiff = ZERO
        if (present(Pmin) .and. present(Pmax)) then
            mode_min = (Pmin+1)*(Pmin+1)*(Pmin+1)
            mode_max = (Pmax+1)*(Pmax+1)*(Pmax+1)
            qdiff(mode_min:mode_max) = qtmp(mode_min:mode_max)
        else
            qdiff = qtmp
        end if


        !
        ! If the current element is being differentiated (ielem == ielem_seed)
        ! then copy the solution modes to local AD variable and seed derivatives
        !
        differentiate_me = ( (elem_info%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (elem_info%ielement_g == fcn_info%seed%ielement_g) )

        if (differentiate_me) then
            ! Differentiating element, initialize appropriate derivatives to ONE
            do iterm = 1,size(qdiff)
                set_deriv = (ieqn - 1)*mesh(idom)%elems(ielem)%nterms_s + iterm
                qdiff(iterm)%xp_ad_(set_deriv) = ONE
            end do
        else
            ! Not differentiating element, set all derivatives to ZERO
            do iterm = 1,size(qdiff)
                qdiff(iterm)%xp_ad_ = ZERO
            end do
        end if



        !
        ! Select interpolation type
        !
        select case (interpolation_type)
            case('value')
                var_gq = matmul(mesh(idom)%elems(ielem)%gq%vol%val,qdiff)

            case('ddx')
                var_gq = matmul(mesh(idom)%elems(ielem)%ddx,qdiff)

            case('ddy')
                var_gq = matmul(mesh(idom)%elems(ielem)%ddy,qdiff)

            case('ddz')
                var_gq = matmul(mesh(idom)%elems(ielem)%ddz,qdiff)

            case default
                user_msg = "interpolate_element_autodiff: The 'interpolation_type' incoming parameter was not&
                            a valid string. Valid strings for 'interpolation_type' include &
                            'value', 'ddx', 'ddy', 'ddz'."
                call chidg_signal_one(FATAL,user_msg,interpolation_type)
        end select


        end associate

    end function interpolate_element_autodiff
    !*****************************************************************************************










    !>  Interpolate variable from polynomial expansion to explicit values at quadrature 
    !!  nodes. The automatic differentiation process really starts here, when the polynomial 
    !!  expansion is evaluated.
    !!
    !!  The interpolation process occurs through a matrix-vector multiplication. That is an 
    !!  interpolation matrix multiplied by a vector of modes from the polynomial expansion. 
    !!  To start the automatic differentiation, the derivative arrays of the values for the 
    !!  polynomial modes must be initialized before any computation. So, before the 
    !!  matrix-vector multiplication.
    !!
    !!  Some interpolation parameters to note that a user might select:
    !!      - interpolation_type:   'value', 'ddx', 'ddy', 'ddz'
    !!      - interpolation_source: ME, NEIGHBOR
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      mesh                    Array of mesh instances.
    !!  @param[in]      face                    Face indices for locating the face in mesh.
    !!  @param[in]      q                       Solution vector
    !!  @param[in]      ieqn                    Index of field being interpolated
    !!  @param[inout]   var_gq                  Autodiff values of field evaluated at gq points
    !!  @param[in]      interpolation_type      Interpolate 'value', 'ddx', 'ddy', 'ddz'
    !!  @param[in]      interpolation_source    ME/NEIGHBOR indicating element to interpolate from
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!  @param[in]      itime                   Index for time step in solution
    !!
    !------------------------------------------------------------------------------------------
    function interpolate_face_autodiff(mesh,q,face_info,fcn_info, ieqn, itime, interpolation_type, interpolation_source) result(var_gq)
        type(mesh_t),           intent(in)              :: mesh(:)
        type(chidg_vector_t),   intent(in)              :: q
        type(face_info_t),      intent(in)              :: face_info
        type(function_info_t),  intent(in)              :: fcn_info
        integer(ik),            intent(in)              :: ieqn
        integer(ik),            intent(in)              :: itime
        character(*),           intent(in)              :: interpolation_type
        integer(ik),            intent(in)              :: interpolation_source

        type(face_info_t)   :: iface_info
        type(recv_t)        :: recv_info

        type(AD_D)                      :: var_gq(mesh(face_info%idomain_l)%elems(face_info%ielement_l)%gq%face%nnodes)
        type(AD_D),         allocatable :: qdiff(:)
        real(rk),           allocatable :: qtmp(:)
        real(rk),           allocatable :: interpolator(:,:)
        character(:),       allocatable :: interpolation_style

        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s, ierr
        logical     :: differentiate_me, conforming_interpolation, chimera_interpolation, parallel_interpolation


        ! Chimera data
        integer(ik)                 :: ndonors, idonor
        logical,    allocatable     :: mask(:)          ! node mask for distributing Chimera quadrature points
        type(AD_D), allocatable     :: var_gq_chimera(:)


        !
        ! Get number of donors for the interpolation
        !
        ndonors = get_face_interpolation_ndonors(mesh,face_info,interpolation_source)


        !
        ! Get interpolation style. Conforming or Chimera
        !
        interpolation_style = get_face_interpolation_style(mesh,face_info,interpolation_source)
        conforming_interpolation = (interpolation_style == 'conforming')
        chimera_interpolation    = (interpolation_style == 'chimera')


        !
        ! Get number of derivatives to initialize for automatic differentiation
        !
        nderiv = get_interpolation_nderiv(mesh,fcn_info)


        !
        ! For each donor element to the interpolation. 
        ! (ndonors could be > 1 for Chimera interpolations)
        !
        do idonor = 1,ndonors

            !
            ! Get face info for face being interpolated to(ME, NEIGHBOR), 
            ! interpolation matrix, and recv data for parallel access
            !
            iface_info   = get_face_interpolation_info(        mesh,face_info,interpolation_source,idonor)
            mask         = get_face_interpolation_mask(        mesh,face_info,interpolation_source,idonor)
            recv_info    = get_face_interpolation_comm(        mesh,face_info,interpolation_source,idonor)
            interpolator = get_face_interpolation_interpolator(mesh,face_info,interpolation_source,idonor,interpolation_type,iface_info)

            parallel_interpolation = (recv_info%comm /= 0)

        

            !
            ! Allocate AD array to store a copy of the solution which starts the differentiation
            !
            if (parallel_interpolation) then
                nterms_s = q%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%nterms()
            else
                nterms_s = mesh(iface_info%idomain_l)%elems(iface_info%ielement_l)%nterms_s
            end if


            !
            ! Allocate solution and derivative arrays for temporary solution variable
            !
            if ( allocated(qdiff) ) deallocate(qdiff)
            allocate(qdiff(nterms_s), stat=ierr)
            if (ierr /= 0) call AllocationError

            do iterm = 1,nterms_s
                qdiff(iterm) = AD_D(nderiv)
            end do


            !
            ! Copy the solution variables from 'q' to 'qdiff'
            !
            if (parallel_interpolation) then
                !qdiff = q%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%getvar(ieqn,itime)
                qtmp = q%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%getvar(ieqn,itime)
            else
                !qdiff = q%dom(iface_info%idomain_l)%vecs(iface_info%ielement_l)%getvar(ieqn,itime)
                qtmp = q%dom(iface_info%idomain_l)%vecs(iface_info%ielement_l)%getvar(ieqn,itime)
            end if


            !
            ! Copy correct number of modes. We use this because there is the possibility
            ! that 'q' could be an auxiliary vector with nterms > nterms_s. For example,
            ! if wall distance was computed for P1, and we are running Navier Stokes P0.
            ! So we take only up to the modes that we need.
            !
            qdiff = qtmp(1:nterms_s)



            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            differentiate_me = ( (iface_info%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                                 (iface_info%ielement_g == fcn_info%seed%ielement_g) )

            if ( differentiate_me ) then
                ! Loop through the terms in qdiff, seed appropriate derivatives to ONE
                do iterm = 1,size(qdiff)
                    ! For the given term, seed its appropriate derivative
                    set_deriv = (ieqn - 1)*nterms_s + iterm
                    qdiff(iterm)%xp_ad_(set_deriv) = ONE
                end do

            else
                ! Loop through the terms in qdiff. Set all derivatives to ZERO
                do iterm = 1,size(qdiff)
                    qdiff(iterm)%xp_ad_ = ZERO
                end do

            end if



            !
            ! Interpolate solution to GQ nodes via matrix-vector multiplication
            !
            if ( conforming_interpolation ) then
                var_gq = matmul(interpolator,  qdiff)


            elseif ( chimera_interpolation ) then
                ! Perform interpolation
                var_gq_chimera = matmul(interpolator,  qdiff)

                ! Scatter chimera nodes to appropriate location in var_gq
                var_gq = unpack(var_gq_chimera,mask,var_gq)

            else
                call chidg_signal(FATAL,"interpolate_face: face interpolation type error")
            end if



        end do ! idonor



    end function interpolate_face_autodiff
    !*****************************************************************************************













    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !----------------------------------------------------------------------------------------
    function interpolate_element_standard(mesh,q,idomain_l,ielement_l,ieqn,itime,interpolation_type) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidg_vector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l
        integer(ik),            intent(in)      :: ielement_l
        integer(ik),            intent(in)      :: ieqn
        integer(ik),            intent(in)      :: itime
        character(len=*),       intent(in)      :: interpolation_type

        real(rk),   dimension(mesh(idomain_l)%elems(ielement_l)%gq%vol%nnodes) :: var_gq


        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the quadrature matrix
        ! with the array of modes for the given variable.
        !
        select case (interpolation_type)
            case('value')
                var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%gq%vol%val, q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn,itime))
            case('ddx')
                var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%ddx, q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn,itime))
            case('ddy')
                var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%ddy, q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn,itime))
            case('ddz')
                var_gq = matmul(mesh(idomain_l)%elems(ielement_l)%ddz, q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn,itime))
            case default
                call chidg_signal(FATAL,"interpolate_element_standard: invalid interpolation_type. Options are 'value', 'ddx', 'ddy', 'ddz'.")
        end select

    end function interpolate_element_standard
    !*****************************************************************************************













    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    function interpolate_face_standard(mesh,q,idomain_l,ielement_l,iface,ieqn,itime) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidg_vector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l, ielement_l, iface, ieqn
        integer(ik),            intent(in)      :: itime

        real(rk),   dimension(mesh(idomain_l)%elems(ielement_l)%gq%face%nnodes) :: var_gq

        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the face quadrature matrix
        ! with the array of modes for the given variable
        !
        var_gq = matmul(mesh(idomain_l)%faces(ielement_l,iface)%gq%face%val(:,:,iface), q%dom(idomain_l)%vecs(ielement_l)%getvar(ieqn,itime))


    end function interpolate_face_standard
    !*****************************************************************************************












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
!    !!  @param[in]      source      ME/NEIGHBOR indicating which element to interpolate from
!    !!
!    !-----------------------------------------------------------------------------------------------------------
!    subroutine interpolate_boundary_autodiff(mesh,face,q,ieqn,points,var)
!        type(mesh_t),           intent(in)              :: mesh(:)
!        type(face_info_t),      intent(in)              :: face
!        type(chidg_vector_t),    intent(in)              :: q
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






    !>  This routine returns face_info data for the face being interpolated to.
    !!
    !!  Given a face, interpolation type, and interpolation source, this routine returns
    !!      - a face_info_t of the face to be interpolated to
    !!      - an interpolation matrix
    !!      - receiver info to be used if the interpolation is remote
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_face_interpolation_info(mesh,face_info,interpolation_source,idonor) result(iface_info)
        type(mesh_t),       intent(in)                  :: mesh(:)
        type(face_info_t),  intent(in)                  :: face_info
        integer(ik),        intent(in)                  :: interpolation_source
        integer(ik),        intent(in)                  :: idonor

        type(face_info_t)   :: iface_info
        integer(ik)         :: ChiID
        logical             :: conforming_interpolation, chimera_interpolation, parallel_interpolation


        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )

        !
        ! Compute neighbor access indices
        !
        if ( interpolation_source == ME ) then

            ! Interpolate from ME element
            iface_info%idomain_l  = face_info%idomain_l
            iface_info%ielement_l = face_info%ielement_l
            iface_info%idomain_g  = face_info%idomain_g
            iface_info%ielement_g = face_info%ielement_g
            iface_info%iface      = face_info%iface

        elseif ( interpolation_source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )

            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then
                iface_info%idomain_g  = mesh(idom)%faces(ielem,iface)%ineighbor_domain_g  
                iface_info%idomain_l  = mesh(idom)%faces(ielem,iface)%ineighbor_domain_l
                iface_info%ielement_g = mesh(idom)%faces(ielem,iface)%ineighbor_element_g
                iface_info%ielement_l = mesh(idom)%faces(ielem,iface)%ineighbor_element_l
                iface_info%iface      = compute_neighbor_face(mesh,idom,ielem,iface,idonor)         ! THIS PROBABLY NEEDS IMPROVED

            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh(idom)%faces(ielem,iface)%ChiID
                iface_info%idomain_g  = mesh(idom)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)
                iface_info%idomain_l  = mesh(idom)%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)
                iface_info%ielement_g = mesh(idom)%chimera%recv%data(ChiID)%donor_element_g%at(idonor)
                iface_info%ielement_l = mesh(idom)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)
            else
                call chidg_signal(FATAL,"get_face_interpolation_info: neighbor conforming_interpolation nor chimera_interpolation were detected")
            end if

        else
            call chidg_signal(FATAL,"get_face_interpolation_info: invalid source. ME or NEIGHBOR.")
        end if



        end associate


    end function get_face_interpolation_info
    !*****************************************************************************************











    !>  This returns an interpolation matrix that is used to actually perform the interpolation 
    !!  from a modal expansion to a set of quadrature nodes.
    !!
    !!  The interpolation from a modal expansion to a set of quadrature nodes takes the form 
    !!  of a matrix-vector multiplication, where the vector entries are modal coefficients of 
    !!  the polynomial expansion, and the matrix contains entries that evaluate the modes of 
    !!  the polynomial expansion at the quadrature nodes. This routine returns the 
    !!  interpolation matrix. Additionally, an interpolation could be computing the actual 
    !!  value of the expansion at the nodes('value'), or it could be computing derivatives 
    !!  ('ddx', 'ddy', 'ddz'). The interpolation_type specifies what kind of interpolation 
    !!  to perform.
    !!
    !!  Given a face, interpolation type, and interpolation source, this routine returns
    !!      - a face_info_t of the face to be interpolated to
    !!      - an interpolation matrix
    !!      - receiver info to be used if the interpolation is remote
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_face_interpolation_interpolator(mesh,source_face,interpolation_source,idonor,interpolation_type,donor_face) result(interpolator)
        type(mesh_t),       intent(in)  :: mesh(:)
        type(face_info_t),  intent(in)  :: source_face
        integer(ik),        intent(in)  :: interpolation_source
        integer(ik),        intent(in)  :: idonor
        type(face_info_t),  intent(in)  :: donor_face
        character(len=*),   intent(in)  :: interpolation_type

        real(rk),       allocatable :: interpolator(:,:)
        integer(ik),    allocatable :: gq_node_indices(:)
        integer(ik)                 :: inode, ChiID, donor_proc
        logical                     :: conforming_interpolation, chimera_interpolation, parallel_interpolation


        associate( idom => source_face%idomain_l, ielem => source_face%ielement_l, iface => source_face%iface )

        !
        ! Compute neighbor access indices
        !
        if ( interpolation_source == ME ) then
            select case(interpolation_type)
                case('value')
                    interpolator = mesh(idom)%faces(ielem,iface)%gq%face%val(:,:,iface)
                case('ddx')
                    interpolator = mesh(idom)%faces(ielem,iface)%ddx
                case('ddy')
                    interpolator = mesh(idom)%faces(ielem,iface)%ddy
                case('ddz')
                    interpolator = mesh(idom)%faces(ielem,iface)%ddz
                case default
                    call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'ddx', 'ddy', 'ddz'.")
            end select



        elseif ( interpolation_source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then

                parallel_interpolation   = ( IRANK /= mesh(idom)%faces(ielem,iface)%ineighbor_proc )

                ! If parallel use iface_interp to get an interpolation for an opposite face. Assumes same order and eqns.
                if (parallel_interpolation) then
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh(idom)%faces(ielem,iface)%gq%face%val(:,:,donor_face%iface)    ! THIS PROBABLY NEEDS IMPROVED
                        case('ddx')
                            interpolator = mesh(idom)%faces(ielem,iface)%neighbor_ddx
                        case('ddy')
                            interpolator = mesh(idom)%faces(ielem,iface)%neighbor_ddy
                        case('ddz')
                            interpolator = mesh(idom)%faces(ielem,iface)%neighbor_ddz
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'ddx', 'ddy', 'ddz'.")
                    end select
                else
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%gq%face%val(:,:,donor_face%iface)
                        case('ddx')
                            interpolator = mesh(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%ddx
                        case('ddy')
                            interpolator = mesh(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%ddy
                        case('ddz')
                            interpolator = mesh(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%ddz
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'ddx', 'ddy', 'ddz'.")
                    end select
                end if


            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh(idom)%faces(ielem,iface)%ChiID
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator%at(idonor)
                        case('ddx')
                            interpolator = mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator_ddx%at(idonor)
                        case('ddy')
                            interpolator = mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator_ddy%at(idonor)
                        case('ddz')
                            interpolator = mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator_ddz%at(idonor)
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'ddx', 'ddy', 'ddz'.")
                    end select

            else
                call chidg_signal(FATAL,"get_face_interpolation_interpolator: neighbor conforming_interpolation nor chimera_interpolation were detected")
            end if



        else
            call chidg_signal(FATAL,"get_face_interpolation_interpolator: invalid source. ME or NEIGHBOR.")
        end if



        end associate


    end function get_face_interpolation_interpolator
    !*****************************************************************************************











    !>  This routine returns a logical mask that indicates which nodes in the full 
    !!  quadrature node set are to be filled by a chimera donor. This defines the scatter 
    !!  pattern for the nodes returned by a chimera donor to the nodes in the full node set.
    !!
    !!  If the interpolation is not a CHIMERA interpolation, then the mask is simply allocated
    !!  to length one and has no meaning.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_face_interpolation_mask(mesh,face_info,interpolation_source,idonor) result(mask)
        type(mesh_t),       intent(in)                  :: mesh(:)
        type(face_info_t),  intent(in)                  :: face_info
        integer(ik),        intent(in)                  :: interpolation_source
        integer(ik),        intent(in)                  :: idonor

        logical, allocatable :: mask(:) !< This gets returned if CHIMERA interpolation


        integer(ik), allocatable    :: gq_node_indices(:)
        integer(ik)                 :: inode, ChiID, nnodes, ierr
        logical                     :: chimera_interpolation

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )


        if ( interpolation_source == NEIGHBOR ) then


            nnodes = mesh(idom)%faces(ielem,iface)%gq%face%nnodes
            allocate(mask(nnodes), stat=ierr) 
            mask = .false.
            if (ierr /= 0) call AllocationError



            chimera_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )

            !
            ! Interpolate from CHIMERA donor element
            !
            if ( chimera_interpolation ) then
                ChiID = mesh(idom)%faces(ielem,iface)%ChiID
                gq_node_indices = mesh(idom)%chimera%recv%data(ChiID)%donor_gq_indices(idonor)%data()

                ! Create mask over full GQ vector of only those nodes that are filled by the current element
                do inode = 1,size(gq_node_indices)
                    mask(gq_node_indices(inode)) = .true.
                end do

            end if


        else
            ! Otherwise, just allocate something so it doens't segfault when it tries to return.
            allocate(mask(1))

        end if

        end associate

    end function get_face_interpolation_mask
    !*****************************************************************************************












    !>  This routine returns a recv_t communication structure indicating whether the 
    !!  interpolation donor is LOCAL or REMOTE. 
    !!
    !!  The recv_t contains recv_comm, recv_domain, and recv_element components. If these are 
    !!  set, then the interpolation is remote and these indices specify where in the 'recv' 
    !!  container to find the solution modes for the interpolation. If these indices are not 
    !!  set, then the interpolation is LOCAL and the main interpolation routine can use the 
    !!  local element indices to locate the solution modes.
    !!  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_face_interpolation_comm(mesh,face_info,interpolation_source,idonor) result(recv_info)
        type(mesh_t),       intent(in)                  :: mesh(:)
        type(face_info_t),  intent(in)                  :: face_info
        integer(ik),        intent(in)                  :: interpolation_source
        integer(ik),        intent(in)                  :: idonor


        type(recv_t)    :: recv_info                !< This gets set if REMOTE interpolation
        integer(ik)     :: ChiID, donor_proc
        logical         :: conforming_interpolation, chimera_interpolation, parallel_interpolation


        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )

        !
        ! Initialize recv_info container to null, indicating LOCAL interpolation. 
        ! Always the case, if interpolation_source=ME.
        !
        recv_info = recv_t(0,0,0)



        if ( interpolation_source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then

                parallel_interpolation   = ( IRANK /= mesh(idom)%faces(ielem,iface)%ineighbor_proc )

                if (parallel_interpolation) then
                    recv_info%comm    = mesh(idom)%faces(ielem,iface)%recv_comm
                    recv_info%domain  = mesh(idom)%faces(ielem,iface)%recv_domain
                    recv_info%element = mesh(idom)%faces(ielem,iface)%recv_element
                end if



            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh(idom)%faces(ielem,iface)%ChiID
                donor_proc = mesh(idom)%chimera%recv%data(ChiID)%donor_proc%at(idonor)

                parallel_interpolation = (IRANK /= donor_proc)
                if (parallel_interpolation) then
                     recv_info%comm    = mesh(idom)%chimera%recv%data(ChiID)%donor_recv_comm%at(idonor)
                     recv_info%domain  = mesh(idom)%chimera%recv%data(ChiID)%donor_recv_domain%at(idonor)
                     recv_info%element = mesh(idom)%chimera%recv%data(ChiID)%donor_recv_element%at(idonor)
                end if

            else
                call chidg_signal(FATAL,"get_face_interpolation_comm: neighbor conforming_interpolation nor chimera_interpolation were detected")
            end if


        end if

        end associate

    end function get_face_interpolation_comm
    !*****************************************************************************************














    !>  Return the number of elements contributing to the interpolation. 
    !!
    !!  For standard conforming interpolations this will be ==1. 
    !!  For Chimera interpolations, this could be >=1.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/17/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_face_interpolation_ndonors(mesh,face_info,interpolation_source) result(ndonors)
        type(mesh_t),       intent(in)  :: mesh(:)
        type(face_info_t),  intent(in)  :: face_info
        integer(ik),        intent(in)  :: interpolation_source

        integer(ik) :: ndonors, ChiID
        logical     :: chimera_interpolation, conforming_interpolation

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )


        !
        ! Test if interpolating from local element
        !
        if ( interpolation_source == ME ) then

            ndonors = 1


        !
        ! Test if interpolating from neighbor element(s)
        !
        elseif (interpolation_source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Test for standard conforming interpolation from neighbor
            if ( conforming_interpolation ) then
                ndonors = 1


            ! Test for chimera interpolation from neighbor
            elseif ( chimera_interpolation ) then
                ChiID   = mesh(idom)%faces(ielem,iface)%ChiID
                ndonors = mesh(idom)%chimera%recv%data(ChiID)%ndonors()



            else
                call chidg_signal(FATAL,"get_face_interpolation_ndonors: invalid value for 'face%ftype'")
            end if

        else
            call chidg_signal(FATAL,"get_face_interpolation_ndonors: invalid value for incoming parameter 'source'")
        end if

        end associate

    end function get_face_interpolation_ndonors
    !*****************************************************************************************






    !>  Return the style of the interpolation.
    !!
    !!  Interpolation styles include:
    !!      - 'conforming'    conforming interpolation
    !!      - 'chimera'       chimera interpolation
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_face_interpolation_style(mesh,face_info,interpolation_source) result(interpolation_style)
        type(mesh_t),       intent(in)  :: mesh(:)
        type(face_info_t),  intent(in)  :: face_info
        integer(ik),        intent(in)  :: interpolation_source

        character(len=:),   allocatable :: interpolation_style
        logical                         :: conforming_interpolation, chimera_interpolation


        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )

        if ( interpolation_source == ME ) then
            conforming_interpolation = ( (mesh(idom)%faces(ielem,iface)%ftype == INTERIOR) .or. &
                                         (mesh(idom)%faces(ielem,iface)%ftype == BOUNDARY) .or. &
                                         (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA) )         ! including chimera here because in the ME case, it doesn't matter
        elseif (interpolation_source == NEIGHBOR ) then
            chimera_interpolation    = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA  )
            conforming_interpolation = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )
        else
            call chidg_signal(FATAL,"get_face_interpolation_style: Invalid interpolation_source. ME or NEIGHBOR")
        end if


        
        if (conforming_interpolation) then
            interpolation_style = "conforming"
        else if (chimera_interpolation) then
            interpolation_style = "chimera"
        else
            call chidg_signal(FATAL,"get_face_interpolation_style: Error in selecting interpolation style")
        end if

        end associate

    end function get_face_interpolation_style
    !******************************************************************************************




    !>  Determine the number of derivatives being computed for the automatic differentiation 
    !!  process.
    !!
    !!  This is the number of degrees of freedom in the element being differentiated with 
    !!  respect to.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_interpolation_nderiv(mesh,function_info) result(nderiv)
        type(mesh_t),           intent(in)  :: mesh(:)
        type(function_info_t),  intent(in)  :: function_info

        integer(ik) :: nderiv, neqns_seed, nterms_s_seed
        logical     :: parallel_seed

        !
        ! If ielem_seed == 0 then we aren't interested in tracking derivatives. So set it 
        ! to lowest number possible while still having something allocated in the AD type 
        ! so the operations are valid.
        !
        if (function_info%seed%ielement_l == 0) then
            nderiv = 1

        else

            !
            ! Compute number of unknowns in the seed element, which is the number of 
            ! partial derivatives we are tracking.
            !
            neqns_seed    = function_info%seed%neqns
            nterms_s_seed = function_info%seed%nterms_s

            nderiv = neqns_seed  *  nterms_s_seed

        end if


    end function get_interpolation_nderiv
    !*****************************************************************************************



end module mod_interpolate







