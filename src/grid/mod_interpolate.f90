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
!!      !interpolate_edge_autodiff 
!!      interpolate_general_autodiff
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
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: CHIMERA, INTERIOR, BOUNDARY, &
                                      ME, NEIGHBOR, ONE, ZERO
                                  
    use mod_chidg_mpi,          only: IRANK
    use mod_DNAD_tools,         only: compute_neighbor_face
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel
    use mod_polynomial,         only: polynomial_val
    use DNAD_D

    use type_mesh,              only: mesh_t
    use type_element_info,      only: element_info_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_edge_info,         only: edge_info_t
    use type_seed,              only: seed_t
    use type_function_info,     only: function_info_t
    use type_recv,              only: recv_t
    use type_chidg_vector,      only: chidg_vector_t
    implicit none



contains




    !>  Interpolate polynomial expansion to volume quadrature node set, initializing the
    !!  automatic differentiation process if needed.
    !!
    !!  This returns an array of automatic differentiation(AD) values at the quadrature nodes,
    !!  that also have their derivatives initialized.
    !!
    !!  Some interpolation parameters to note that are used inside here:
    !!      - interpolation_type:   'value', 'grad1', 'grad2', 'grad3'
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
    !!  @param[in]  interpolator: optionally, pass in another interpolator to interpolate to 
    !!              a set of nodes other than the default set of nodes.
    !!
    !-----------------------------------------------------------------------------------------
    function interpolate_element_autodiff(mesh,vector,elem_info,fcn_info,ifield,itime,interpolation_type,interpolator,mode_min,mode_max) result(var_gq)
        type(mesh_t),           intent(in)              :: mesh
        type(chidg_vector_t),   intent(in)              :: vector
        type(element_info_t),   intent(in)              :: elem_info
        type(function_info_t),  intent(in)              :: fcn_info
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime
        character(*),           intent(in)              :: interpolation_type
        real(rk),               intent(in), optional    :: interpolator(:,:)
        integer(ik),            intent(in), optional    :: mode_min
        integer(ik),            intent(in), optional    :: mode_max


        character(:),   allocatable :: user_msg

        type(AD_D), allocatable :: qdiff(:), var_gq(:)
        real(rk),   allocatable :: qtmp(:)

        integer(ik) :: nderiv, set_deriv, iterm, nterms, ierr
        logical     :: differentiate_me

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l )


        user_msg = "interpolate_element_autodiff: currently only supports overriding the &
                    default interpolator for 'value' interpolations."
        if (present(interpolator) .and. (interpolation_type /= 'value')) call chidg_signal(FATAL,user_msg)



        !
        ! Get number of derivatives to initialize for automatic differentiation
        !
        nderiv = get_interpolation_nderiv(mesh,fcn_info%seed)


        !
        ! Determine nterms from interpolator matrix
        !
        nterms = mesh%domain(idom)%elems(ielem)%basis_s%nterms_i()

        ! Allocate qdiff buffer and initialize derivative allocations + zero
        allocate(qdiff(nterms), stat=ierr)
        if (ierr /= 0) call AllocationError
        do iterm = 1,nterms
            qdiff(iterm) = AD_D(nderiv)
        end do
        qdiff = ZERO


        !
        ! Retrieve modal coefficients representing ifield in vector to 'qtmp'
        !
        qtmp = vector%dom(idom)%vecs(ielem)%getvar(ifield,itime)


        !
        ! 1: If mode_min,mode_max are present, use a subset of the expansion modes.
        ! 2: If not present, use the entire expansion.
        !
        if (present(mode_min) .and. present(mode_max)) then
            qdiff(mode_min:mode_max) = qtmp(mode_min:mode_max)
        else
            qdiff(1:nterms) = qtmp(1:nterms)
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
                set_deriv = (ifield - 1)*mesh%domain(idom)%elems(ielem)%nterms_s + iterm
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
                if (present(interpolator)) then
                    var_gq = matmul(interpolator,qdiff)
                else
                    var_gq = matmul(mesh%domain(idom)%elems(ielem)%basis_s%interpolator_element('Value'),qdiff)
                end if

            case('grad1')
                var_gq = matmul(mesh%domain(idom)%elems(ielem)%grad1,qdiff)

            case('grad2')
                var_gq = matmul(mesh%domain(idom)%elems(ielem)%grad2,qdiff)

            case('grad3')
                var_gq = matmul(mesh%domain(idom)%elems(ielem)%grad3,qdiff)

            case default
                user_msg = "interpolate_element_autodiff: The 'interpolation_type' incoming&
                            parameter was not a valid string. Valid strings for &
                            'interpolation_type' include 'value', 'grad1', 'grad2', 'grad3'."
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
    !!      - interpolation_type:   'value', 'grad1', 'grad2', 'grad3'
    !!      - interpolation_source: ME, NEIGHBOR
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      mesh                    Array of mesh instances.
    !!  @param[in]      face                    Face indices for locating the face in mesh.
    !!  @param[in]      vector                  A chidg_vector containing a modal representation of fields on mesh
    !!  @param[in]      ifield                  Index of field being interpolated
    !!  @param[inout]   var_gq                  Autodiff values of field evaluated at gq points
    !!  @param[in]      interpolation_type      Interpolate 'value', 'grad1', 'grad2', 'grad3'
    !!  @param[in]      interpolation_source    ME/NEIGHBOR indicating element to interpolate from
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !!  @param[in]      itime                   Index for time step in solution
    !!
    !------------------------------------------------------------------------------------------
    function interpolate_face_autodiff(mesh,vector,face_info,fcn_info, ifield, itime, interpolation_type, interpolation_source) result(var_gq)
        type(mesh_t),           intent(in)              :: mesh
        type(chidg_vector_t),   intent(in)              :: vector
        type(face_info_t),      intent(in)              :: face_info
        type(function_info_t),  intent(in)              :: fcn_info
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime
        character(*),           intent(in)              :: interpolation_type
        integer(ik),            intent(in)              :: interpolation_source

        type(face_info_t)   :: iface_info
        type(recv_t)        :: recv_info

        type(AD_D),         allocatable :: qdiff(:), var_gq(:)
        real(rk),           allocatable :: qtmp(:)
        real(rk),           allocatable :: interpolator(:,:)
        character(:),       allocatable :: interpolation_style

        integer(ik) :: nderiv, set_deriv, iterm, igq, nterms_s, ierr, nnodes
        logical     :: differentiate_me, conforming_interpolation, chimera_interpolation, parallel_interpolation


        ! Chimera data
        integer(ik)                 :: ndonors, idonor
        logical,    allocatable     :: mask(:)          ! node mask for distributing Chimera quadrature points
        type(AD_D), allocatable     :: var_gq_chimera(:)

        
        !
        ! Allocate output array
        !
        nnodes   = mesh%domain(face_info%idomain_l)%elems(face_info%ielement_l)%basis_s%nnodes_face()
        nterms_s = mesh%domain(face_info%idomain_l)%elems(face_info%ielement_l)%basis_s%nterms_i()
        allocate(var_gq(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

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
        nderiv = get_interpolation_nderiv(mesh,fcn_info%seed)




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
            ! Deduce nterms from interpolator matrix
            !
            nterms_s = size(interpolator,2)


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
            ! Retrieve modal coefficients for ifield from vector
            !
            if (parallel_interpolation) then
                qtmp = vector%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%getvar(ifield,itime)
            else
                qtmp = vector%dom(iface_info%idomain_l)%vecs(iface_info%ielement_l)%getvar(ifield,itime)
            end if


            !
            ! Copy correct number of modes. We use this because there is the possibility
            ! that 'q' could be an auxiliary vector with nterms > nterms_s. For example,
            ! if wall distance was computed for P1, and we are running Navier Stokes P0.
            ! So we take only up to the modes that we need.
            !
            qdiff(1:nterms_s) = qtmp(1:nterms_s)



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
                    set_deriv = (ifield - 1)*nterms_s + iterm
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









!    !>  Interpolate variable from polynomial expansion to explicit values at quadrature 
!    !!  nodes. The automatic differentiation process really starts here, when the polynomial 
!    !!  expansion is evaluated.
!    !!
!    !!  The interpolation process occurs through a matrix-vector multiplication. That is an 
!    !!  interpolation matrix multiplied by a vector of modes from the polynomial expansion. 
!    !!  To start the automatic differentiation, the derivative arrays of the values for the 
!    !!  polynomial modes must be initialized before any computation. So, before the 
!    !!  matrix-vector multiplication.
!    !!
!    !!  Some interpolation parameters to note that a user might select:
!    !!      - interpolation_type:   'value', 'grad1', 'grad2', 'grad3'
!    !!      - interpolation_source: ME, NEIGHBOR
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/1/2016
!    !!
!    !!  @param[in]      mesh                    Array of mesh instances.
!    !!  @param[in]      face                    Face indices for locating the face in mesh.
!    !!  @param[in]      vector                  A chidg_vector containing a modal representation of fields on mesh
!    !!  @param[in]      ifield                  Index of field being interpolated
!    !!  @param[inout]   var_gq                  Autodiff values of field evaluated at gq points
!    !!  @param[in]      interpolation_type      Interpolate 'value', 'grad1', 'grad2', 'grad3'
!    !!  @param[in]      interpolation_source    ME/NEIGHBOR indicating element to interpolate from
!    !!
!    !!  @author Mayank Sharma + Matteo Ugolotti
!    !!  @date   11/5/2016
!    !!
!    !!  @param[in]      itime                   Index for time step in solution
!    !!
!    !------------------------------------------------------------------------------------------
!    function interpolate_edge_autodiff(mesh,vector,edge_info,seed, ifield, itime, interpolation_type) result(var_gq)
!        type(mesh_t),           intent(in)              :: mesh
!        type(chidg_vector_t),   intent(in)              :: vector
!        type(edge_info_t),      intent(in)              :: edge_info
!        type(seed_t),           intent(in)              :: seed
!        integer(ik),            intent(in)              :: ifield
!        integer(ik),            intent(in)              :: itime
!        character(*),           intent(in)              :: interpolation_type
!
!        type(AD_D),         allocatable :: qdiff(:), var_gq(:)
!        real(rk),           allocatable :: qtmp(:)
!        real(rk),           allocatable :: interpolator(:,:)
!
!        integer(ik) :: nderiv, set_deriv, iterm, nterms_s, ierr, nnodes
!        logical     :: differentiate_me
!
!        
!        !
!        ! Allocate output array
!        !
!        nnodes   = mesh%domain(edge_info%idomain_l)%elems(edge_info%ielement_l)%basis_s%nnodes_face()
!        nterms_s = mesh%domain(edge_info%idomain_l)%elems(edge_info%ielement_l)%basis_s%nterms_i()
!        allocate(var_gq(nnodes), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        !
!        ! Get number of derivatives to initialize for automatic differentiation
!        !
!        nderiv = get_interpolation_nderiv(mesh,seed)
!
!
!        !
!        ! Get edge interpolation matrix
!        !
!        interpolator = get_edge_interpolation_interpolator(mesh,edge_info,interpolation_type)
!
!
!        !
!        ! Deduce nterms from interpolator matrix
!        !
!        nterms_s = size(interpolator,2)
!
!
!        !
!        ! Allocate solution and derivative arrays for temporary solution variable
!        !
!        if ( allocated(qdiff) ) deallocate(qdiff)
!        allocate(qdiff(nterms_s), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        do iterm = 1,nterms_s
!            qdiff(iterm) = AD_D(nderiv)
!        end do
!
!
!        !
!        ! Retrieve modal coefficients for ifield from vector
!        !
!        qtmp = vector%dom(edge_info%idomain_l)%vecs(edge_info%ielement_l)%getvar(ifield,itime)
!
!
!        !
!        ! Copy correct number of modes. We use this because there is the possibility
!        ! that 'q' could be an auxiliary vector with nterms > nterms_s. For example,
!        ! if wall distance was computed for P1, and we are running Navier Stokes P0.
!        ! So we take only up to the modes that we need.
!        !
!        qdiff(1:nterms_s) = qtmp(1:nterms_s)
!
!
!
!        !
!        ! If the current element is being differentiated (ielem == ielem_seed)
!        ! then copy the solution modes to local AD variable and seed derivatives
!        !
!        differentiate_me = ( (edge_info%idomain_g  == seed%idomain_g ) .and. &
!                             (edge_info%ielement_g == seed%ielement_g) )
!
!        if ( differentiate_me ) then
!            ! Loop through the terms in qdiff, seed appropriate derivatives to ONE
!            do iterm = 1,size(qdiff)
!                ! For the given term, seed its appropriate derivative
!                set_deriv = (ifield - 1)*nterms_s + iterm
!                qdiff(iterm)%xp_ad_(set_deriv) = ONE
!            end do
!
!        else
!            ! Loop through the terms in qdiff. Set all derivatives to ZERO
!            do iterm = 1,size(qdiff)
!                qdiff(iterm)%xp_ad_ = ZERO
!            end do
!
!        end if
!
!
!
!        !
!        ! Interpolate solution to GQ nodes via matrix-vector multiplication
!        !
!        var_gq = matmul(interpolator,  qdiff)
!
!
!    end function interpolate_edge_autodiff
!    !*****************************************************************************************










    !> Compute variable at quadrature nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !----------------------------------------------------------------------------------------
    function interpolate_element_standard(mesh,q,idomain_l,ielement_l,ifield,itime,interpolation_type) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh
        type(chidg_vector_t),   intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l
        integer(ik),            intent(in)      :: ielement_l
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime
        character(*),           intent(in)      :: interpolation_type

        real(rk),   allocatable :: var_gq(:)


        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the quadrature matrix
        ! with the array of modes for the given variable.
        !
        select case (interpolation_type)
            case('value')
                var_gq = matmul(mesh%domain(idomain_l)%elems(ielement_l)%basis_s%interpolator_element('Value'), q%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime))
            case('grad1')
                var_gq = matmul(mesh%domain(idomain_l)%elems(ielement_l)%grad1,      q%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime))
            case('grad2')
                var_gq = matmul(mesh%domain(idomain_l)%elems(ielement_l)%grad2,      q%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime))
            case('grad3')
                var_gq = matmul(mesh%domain(idomain_l)%elems(ielement_l)%grad3,      q%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime))
            case default
                call chidg_signal(FATAL,"interpolate_element_standard: invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
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
    function interpolate_face_standard(mesh,q,idomain_l,ielement_l,iface,ifield,itime) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh
        type(chidg_vector_t),   intent(in)      :: q
        integer(ik),            intent(in)      :: idomain_l, ielement_l, iface, ifield
        integer(ik),            intent(in)      :: itime

        real(rk),   allocatable :: var_gq(:)

        !
        ! Use quadrature instance to compute variable at quadrature nodes.
        ! This takes the form of a matrix multiplication of the face quadrature matrix
        ! with the array of modes for the given variable
        !
        var_gq = matmul(mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%interpolator_face('Value',iface), q%dom(idomain_l)%vecs(ielement_l)%getvar(ifield,itime))


    end function interpolate_face_standard
    !*****************************************************************************************











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
    !!  @param[in]      source      ME/NEIGHBOR indicating which element to interpolate from
    !!
    !-----------------------------------------------------------------------------------------------------------
    function interpolate_general_autodiff(mesh,vector,fcn_info,ifield,itime,interpolation_type,physical_nodes,try_offset,donors,donor_coords) result(var)
        type(mesh_t),           intent(in)              :: mesh
        type(chidg_vector_t),   intent(in)              :: vector
        type(function_info_t),  intent(in)              :: fcn_info
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime
        character(*),           intent(in)              :: interpolation_type
        real(rk),               intent(in)              :: physical_nodes(:,:)
        real(rk),               intent(in), optional    :: try_offset(3)
        type(element_info_t),   intent(in), optional    :: donors(:)
        real(rk),               intent(in), optional    :: donor_coords(:,:)

        type(element_info_t)    :: donor
        type(recv_t)            :: recv_info

        type(AD_D), allocatable  :: qdiff(:), var(:), tmp(:)
        real(rk),   allocatable  :: interpolator(:,:)

        real(rk)        :: donor_coord(3), donor_volume
        integer(ik)     :: nderiv, set_deriv, iterm, ierr, nterms_s, inode, nnodes
        logical         :: differentiate_me = .false.
        logical         :: donor_found      = .false.
        logical         :: parallel_donor   = .false.

        
        !
        ! 1: Check incoming node array makes sense
        ! 2: Check interpolation_type
        !
        if (size(physical_nodes,2) /= 3) call chidg_signal(FATAL,'interpolate_general_autodiff: size(physical_nodes,2) /= 3.')
        if (interpolation_type /= 'value') call chidg_signal(FATAL,"interpolation_general_autodiff: currently only supports interpolation_type == 'value'.")


        !
        ! Get number of derivatives to initialize for automatic differentiation
        !
        nderiv = get_interpolation_nderiv(mesh,fcn_info%seed)


        !
        ! Allocate result and derivatives
        !
        nnodes = size(physical_nodes,1)
        allocate(var(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        do inode = 1,size(var)
            allocate(var(inode)%xp_ad_(nderiv))
        end do




        !
        ! Find donor elements for incoming nodes
        !
        do inode = 1,size(physical_nodes,1)

            if (present(donors) .and. present(donor_coords)) then
                donor       = donors(inode)
                donor_coord = donor_coords(inode,:)

            else
                !
                ! Try processor-LOCAL elements
                !
                call find_gq_donor(mesh,                                &
                                   physical_nodes(inode,1:3),           &
                                   [ZERO,ZERO,ZERO],                    &
                                   face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                   donor,                               &
                                   donor_coord,                         &
                                   donor_found,                         &
                                   donor_volume=donor_volume)

                !
                ! Try LOCAL elements with try_offset if still not found and try_offset 
                ! is present
                !
                if ( (.not. donor_found) .and. (present(try_offset)) ) then
                    call find_gq_donor(mesh,                                &
                                       physical_nodes(inode,1:3),           &
                                       try_offset,                          &
                                       face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                       donor,                               &
                                       donor_coord,                         &
                                       donor_found,                         &
                                       donor_volume=donor_volume)

                end if



                !
                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                !
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                                &
                                                physical_nodes(inode,1:3),           &
                                                [ZERO,ZERO,ZERO],                    &
                                                face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                                donor,                               &
                                                donor_coord,                         &
                                                donor_found,                         &
                                                donor_volume=donor_volume)
                end if

                
                !
                ! Try PARALLEL_ELEMENTS with try_offset if still not found and
                ! try_offset is present
                !
                if ( (.not. donor_found) .and. (present(try_offset)) ) then
                    call find_gq_donor_parallel(mesh,                                &
                                                physical_nodes(inode,1:3),           &
                                                try_offset,                          &
                                                face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                                donor,                               &
                                                donor_coord,                         &
                                                donor_found,                         &
                                                donor_volume=donor_volume)
                end if 




                ! Abort if we didn't find a donor
                if (.not. donor_found) call chidg_signal(FATAL,"interpolate_general_autodiff: no donor element found for interpolation node.")

            end if ! Find Donor


            ! Check parallel or local donor
            parallel_donor = (donor%iproc /= IRANK)



            !
            ! Get nterms
            !
            if (parallel_donor) then
                nterms_s = mesh%parallel_element(donor%pelem_ID)%nterms_s
            else
                nterms_s = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%nterms_s
            end if

            


            !
            ! Construct interpolator
            !
            if (allocated(interpolator)) deallocate(interpolator)
            allocate(interpolator(1,nterms_s), stat=ierr)
            if (ierr /= 0) call AllocationError


            do iterm = 1,nterms_s
                interpolator(1,iterm) = polynomial_val(3,nterms_s,iterm,donor_coord)
            end do ! iterm



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
            ! Retrieve modal coefficients for ifield from vector
            !
            if (parallel_donor) then
                recv_info%comm    = mesh%parallel_element(donor%pelem_ID)%recv_comm
                recv_info%domain  = mesh%parallel_element(donor%pelem_ID)%recv_domain
                recv_info%element = mesh%parallel_element(donor%pelem_ID)%recv_element
                qdiff = vector%recv%comm(recv_info%comm)%dom(recv_info%domain)%vecs(recv_info%element)%getvar(ifield,itime)
            else
                qdiff = vector%dom(donor%idomain_l)%vecs(donor%ielement_l)%getvar(ifield,itime)
            end if


            !
            ! If the current element is being differentiated (ielem == ielem_seed)
            ! then copy the solution modes to local AD variable and seed derivatives
            !
            differentiate_me = ( (donor%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                                 (donor%ielement_g == fcn_info%seed%ielement_g) )

            if ( differentiate_me ) then
                ! Loop through the terms in qdiff, seed appropriate derivatives to ONE
                do iterm = 1,size(qdiff)
                    ! For the given term, seed its appropriate derivative
                    set_deriv = (ifield - 1)*nterms_s + iterm
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
            tmp = matmul(interpolator,  qdiff)
            var(inode) = tmp(1)



        end do ! inode


    end function interpolate_general_autodiff
    !*************************************************************************************************************










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
!                interpolator(1,iterm) = polynomial_val(3,nterms,iterm,node)
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
        type(mesh_t),       intent(in)                  :: mesh
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

            chimera_interpolation    = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )

            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then
                iface_info%idomain_g  = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_g  
                iface_info%idomain_l  = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_l
                iface_info%ielement_g = mesh%domain(idom)%faces(ielem,iface)%ineighbor_element_g
                iface_info%ielement_l = mesh%domain(idom)%faces(ielem,iface)%ineighbor_element_l
                iface_info%iface      = compute_neighbor_face(mesh,idom,ielem,iface,idonor)         ! THIS PROBABLY NEEDS IMPROVED

            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
                iface_info%idomain_g  = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_g
                iface_info%idomain_l  = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_l
                iface_info%ielement_g = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_g
                iface_info%ielement_l = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_l
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
    !!  ('grad1', 'grad2', 'grad3'). The interpolation_type specifies what kind of interpolation 
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
        type(mesh_t),       intent(in)  :: mesh
        type(face_info_t),  intent(in)  :: source_face
        integer(ik),        intent(in)  :: interpolation_source
        integer(ik),        intent(in)  :: idonor
        type(face_info_t),  intent(in)  :: donor_face
        character(*),       intent(in)  :: interpolation_type

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
                    interpolator = mesh%domain(idom)%faces(ielem,iface)%basis_s%interpolator_face('Value',iface)
                case('grad1')
                    interpolator = mesh%domain(idom)%faces(ielem,iface)%grad1
                case('grad2')
                    interpolator = mesh%domain(idom)%faces(ielem,iface)%grad2
                case('grad3')
                    interpolator = mesh%domain(idom)%faces(ielem,iface)%grad3
                case default
                    call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
            end select




        elseif ( interpolation_source == NEIGHBOR ) then

            chimera_interpolation    = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then

                parallel_interpolation   = ( IRANK /= mesh%domain(idom)%faces(ielem,iface)%ineighbor_proc )

                ! If parallel use iface_interp to get an interpolation for an opposite face. Assumes same order and eqns.
                if (parallel_interpolation) then
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh%domain(idom)%faces(ielem,iface)%basis_s%interpolator_face('Value',donor_face%iface)    ! THIS PROBABLY NEEDS IMPROVED
                        case('grad1')
                            interpolator = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad1
                        case('grad2')
                            interpolator = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad2
                        case('grad3')
                            interpolator = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad3
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
                    end select
                else
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh%domain(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%basis_s%interpolator_face('Value',donor_face%iface)
                        case('grad1')
                            interpolator = mesh%domain(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%grad1
                        case('grad2')
                            interpolator = mesh%domain(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%grad2
                        case('grad3')
                            interpolator = mesh%domain(donor_face%idomain_l)%faces(donor_face%ielement_l,donor_face%iface)%grad3
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
                    end select
                end if


            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
                    select case(interpolation_type)
                        case('value')
                            interpolator = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%value
                        case('grad1')
                            interpolator = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad1
                        case('grad2')
                            interpolator = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad2
                        case('grad3')
                            interpolator = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad3
                        case default
                            call chidg_signal(FATAL,"get_face_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
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
        type(mesh_t),       intent(in)                  :: mesh
        type(face_info_t),  intent(in)                  :: face_info
        integer(ik),        intent(in)                  :: interpolation_source
        integer(ik),        intent(in)                  :: idonor

        logical, allocatable :: mask(:) !< This gets returned if CHIMERA interpolation


        integer(ik), allocatable    :: gq_node_indices(:)
        integer(ik)                 :: inode, ChiID, nnodes, ierr
        logical                     :: chimera_interpolation

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )


        if ( interpolation_source == NEIGHBOR ) then


            nnodes = mesh%domain(idom)%faces(ielem,iface)%basis_s%nnodes_face()
            allocate(mask(nnodes), stat=ierr) 
            mask = .false.
            if (ierr /= 0) call AllocationError



            chimera_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )

            !
            ! Interpolate from CHIMERA donor element
            !
            if ( chimera_interpolation ) then
                ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
                gq_node_indices = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%node_index(:)

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
        type(mesh_t),       intent(in)                  :: mesh
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

            chimera_interpolation    = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Interpolate from conforming NEIGHBOR element
            if ( conforming_interpolation ) then

                parallel_interpolation   = ( IRANK /= mesh%domain(idom)%faces(ielem,iface)%ineighbor_proc )

                if (parallel_interpolation) then
                    recv_info%comm    = mesh%domain(idom)%faces(ielem,iface)%recv_comm
                    recv_info%domain  = mesh%domain(idom)%faces(ielem,iface)%recv_domain
                    recv_info%element = mesh%domain(idom)%faces(ielem,iface)%recv_element
                end if



            ! Interpolate from CHIMERA donor element
            elseif ( chimera_interpolation ) then
                ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
                donor_proc = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%iproc

                parallel_interpolation = (IRANK /= donor_proc)
                if (parallel_interpolation) then
                     recv_info%comm    = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%recv_comm
                     recv_info%domain  = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%recv_domain
                     recv_info%element = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%recv_element
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
        type(mesh_t),       intent(in)  :: mesh
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

            chimera_interpolation    = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
            conforming_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )


            ! Test for standard conforming interpolation from neighbor
            if ( conforming_interpolation ) then
                ndonors = 1


            ! Test for chimera interpolation from neighbor
            elseif ( chimera_interpolation ) then
                ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
                ndonors = mesh%domain(idom)%chimera%recv(ChiID)%ndonors()


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
        type(mesh_t),       intent(in)  :: mesh
        type(face_info_t),  intent(in)  :: face_info
        integer(ik),        intent(in)  :: interpolation_source

        character(len=:),   allocatable :: interpolation_style
        logical                         :: conforming_interpolation, chimera_interpolation


        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )

        if ( interpolation_source == ME ) then
            conforming_interpolation = ( (mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR) .or. &
                                         (mesh%domain(idom)%faces(ielem,iface)%ftype == BOUNDARY) .or. &
                                         (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA) )         ! including chimera here because in the ME case, it doesn't matter
        elseif (interpolation_source == NEIGHBOR ) then
            chimera_interpolation    = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA  )
            conforming_interpolation = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )
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









!    !>  This returns an interpolation matrix that is used to actually perform the interpolation 
!    !!  from a modal expansion to a set of interpolation nodes.
!    !!
!    !!  The interpolation from a modal expansion to a set of quadrature nodes takes the form 
!    !!  of a matrix-vector multiplication, where the vector entries are modal coefficients of 
!    !!  the polynomial expansion, and the matrix contains entries that evaluate the modes of 
!    !!  the polynomial expansion at the interpolation nodes. This routine returns the 
!    !!  interpolation matrix. Additionally, an interpolation could be computing the actual 
!    !!  value of the expansion at the nodes('value'), or it could be computing derivatives 
!    !!  ('grad1', 'grad2', 'grad3'). The interpolation_type specifies what kind of interpolation 
!    !!  to perform.
!    !!
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   12/6/2017
!    !!
!    !!
!    !----------------------------------------------------------------------------------------
!    function get_edge_interpolation_interpolator(mesh,edge_info,interpolation_type) result(interpolator)
!        type(mesh_t),       intent(in)  :: mesh
!        type(edge_info_t),  intent(in)  :: edge_info
!        character(*),       intent(in)  :: interpolation_type
!
!        real(rk), allocatable :: interpolator(:,:)
!
!        associate( idom  => edge_info%idomain_l,  &
!                   ielem => edge_info%ielement_l, &
!                   iedge => edge_info%iedge )
!
!        !
!        ! Compute neighbor access indices
!        !
!        select case(interpolation_type)
!            case('value')
!                interpolator = mesh%domain(idom)%elems(ielem)%basis_s%interpolator_edge('Value',iedge)
!            case('grad1')
!                interpolator = mesh%domain(idom)%elems(ielem)%edge_grad1(:,:,iedge)
!            case('grad2')
!                interpolator = mesh%domain(idom)%elems(ielem)%edge_grad2(:,:,iedge)
!            case('grad3')
!                interpolator = mesh%domain(idom)%elems(ielem)%edge_grad3(:,:,iedge)
!            case default
!                call chidg_signal(FATAL,"get_edge_interpolation_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
!        end select
!
!
!        end associate
!
!
!    end function get_edge_interpolation_interpolator
!    !*****************************************************************************************







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
    function get_interpolation_nderiv(mesh,seed) result(nderiv)
        type(mesh_t),           intent(in)  :: mesh
        type(seed_t),           intent(in)  :: seed

        integer(ik) :: nderiv

        !
        ! If ielem_seed == 0 then we aren't interested in tracking derivatives. 
        !   !So set it to lowest number possible while still having something 
        !   !allocated in the AD type so the operations are valid.
        !
        !   Actually, allocating with size 0 is okay in fortran. Just need to be
        !   careful not to try and access anything as var%xp_ad_(1)
        !
        if (seed%ielement_l == 0) then
            !nderiv = 1
            nderiv = 0

        else

            !
            ! Compute number of unknowns in the seed element, which is the number of 
            ! partial derivatives we are tracking.
            !
            !neqns_seed    = seed%neqns
            !nterms_s_seed = seed%nterms_s
            !nderiv        = neqns_seed  *  nterms_s_seed
            nderiv = seed%neqns * seed%nterms_s

        end if


    end function get_interpolation_nderiv
    !*****************************************************************************************



end module mod_interpolate







