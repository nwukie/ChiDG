!>  ChiDG Differentiation Utilities
!!
!!  Purpose: 
!!  ---------------------------------------
!!  Allow differentiation of geometrical 
!!  entitites such as jinv, normals, interpolators
!!  for grid-nodes sensitivity calculations
!!  If grid-sensitivity is not necessary it
!!  returns AD_D type with number of derivatives
!!  equal to the number of degree of freedom 
!!  of the face/element (typically nterms_s*nfields)
!!  with null derivatives
!!
!!
!!  Public routines for automatic differentiation 
!!  ---------------------------------------
!!      get_nderiv 
!!      differentiate_modal_coefficients 
!!      differentiate_jinv
!!      differentiate_normal
!!      differentiate_unit_normal
!!      differentiate_coordinate
!!      differentiate_element_size
!!      differentiate_element_interpolator
!!      differentiate_face_interior_interpolator
!!      differentiate_face_parallel_interpolator
!!      differentiate_face_local_interpolator
!!      differentiate_face_chimera_interpolator  
!!
!---------------------------------------------------------------------------------------------
module mod_differentiate
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: CHIMERA, INTERIOR, BOUNDARY, &
                                      ME, NEIGHBOR, ONE, ZERO, TWO,&
                                      NO_DIFF, dQ_DIFF, dX_DIFF,   &
                                      dBC_DIFF, dD_DIFF
                                  
    use mod_chidg_mpi,          only: IRANK
    use mod_DNAD_tools,         only: compute_neighbor_face
    use DNAD_D

    use type_mesh,              only: mesh_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t
    use type_recv,              only: recv_t
    use type_chidg_vector,      only: chidg_vector_t
    use ieee_arithmetic
    implicit none


contains


    
    !>  Return the number of derivatives based on the type of differentiation
    !!
    !!  NOTE: this is an updated copy of get_interpolation_nderiv 
    !!        that was originally in mod_interpolate.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/6/2018
    !!
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/10/2019
    !!
    !!  Introduced differentiation wrt distance field
    !!
    !--------------------------------------------------------------------------------
    function get_nderiv(function_info) result(nderiv)
        type(function_info_t),  intent(in)  :: function_info

        integer(ik) :: nderiv, nfields_seed, nterms_s_seed, nnodes_r, directions
        logical     :: parallel_seed
        
        if (function_info%dtype == NO_DIFF) then
            nderiv = 0

        else if (function_info%dtype == dQ_DIFF) then
            ! Compute number of unknowns in the seed element, which is the number of 
            ! partial derivatives we are tracking.
            nfields_seed  = function_info%seed%nfields
            nterms_s_seed = function_info%seed%nterms_s

            nderiv = nfields_seed * nterms_s_seed

        else if (function_info%dtype == dD_DIFF) then
            ! Compute number of unknowns in the seed element, which is the number of 
            ! partial derivatives we are tracking.
            !
            ! For wall_distance field we have only one variables, whereas the 
            ! nterms are equal to primal problem
            nfields_seed    = 1
            nterms_s_seed = function_info%seed%nterms_s

            nderiv = nfields_seed  *  nterms_s_seed

        else if (function_info%dtype == dX_DIFF) then
            ! The number of derivatives is equal to the number of the reference nodes
            ! of the element seed times the number of directions (3)
            nnodes_r   = function_info%seed%nnodes_r 
            directions = 3
            
            nderiv = nnodes_r * directions
        
        else if (function_info%dtype == dBC_DIFF) then
            ! The number of derivatives is equal to 1, only one BC property
            nderiv = 1 

        else
             call chidg_signal_one(FATAL,"mod_differentiate%get_nderiv: unexpected differentiate type",function_info%dtype)
        end if



    end function get_nderiv
    !********************************************************************************

   
   
   
   
   
   
   
    
    !>  Initialize the derivatives of the modal coefficients. 
    !!
    !!  NOTE: This becomes necessary in case since we can differentiate
    !!        not only wrt the primary modal coefficient but also respect
    !!        to the auxiliary modal coefficients. 
    !!
    !!  Here the possible combinations
    !!      primary variables   wrt primary variables   -> initialize non null derivatives
    !!      primary variables   wrt auxiliary variables -> initialize null derivatives
    !!      auxiliary variables wrt primary variables   -> initialize null derivatives
    !!      auxiliary variables wrt auxiliary variables -> initialize non null derivatives
    !!
    !!      everything that is differentiated wrt boundary conditions or it is NOT
    !!      differentiated will have NULL derivatives
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/15/2019
    !!
    !!  TODO: create tests
    !!
    !--------------------------------------------------------------------------------
    function differentiate_modal_coefficients(q_in,fcn_info,elem_info,ifield,itime,nterms_s,min_mode,max_mode) result(q_out)
        real(rk),                   intent(in)  :: q_in(:)
        type(function_info_t),      intent(in)  :: fcn_info
        type(element_info_t),       intent(in)  :: elem_info
        integer(ik),                intent(in)  :: ifield
        integer(ik),                intent(in)  :: itime
        integer(ik),                intent(in)  :: nterms_s
        integer(ik),    optional,   intent(in)  :: min_mode
        integer(ik),    optional,   intent(in)  :: max_mode

        type(AD_D),     allocatable     :: q_out(:)
        logical                         :: differentiate
        integer(ik)                     :: iterm, ierr, set_deriv, nderiv

        ! Get number of derivatives
        nderiv = get_nderiv(fcn_info)

        ! Allocate q_out buffer and initialize derivative allocations + zero 
        allocate(q_out(nterms_s), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Allocate AD_D derivatives
        do iterm = 1,nterms_s
            q_out(iterm) = AD_D(nderiv)
        end do

        ! Assign real values of the modes to AD_D type
        if (present(min_mode) .and.present(max_mode)) then
            q_out(min_mode:max_mode) = q_in
        else
            q_out = q_in
        end if


        ! Select cases
        ! NOTE: Nathan does not like this 'overloading' of 'type' and 'dtype'. Removed for now
        !differentiate = ( (fcn_info%type        == fcn_info%dtype)           .and. &
        !                  (elem_info%idomain_g  == fcn_info%seed%idomain_g)  .and. &
        !                  (elem_info%ielement_g == fcn_info%seed%ielement_g) .and. &
        !                  (itime                == fcn_info%seed%itime) ) 


        differentiate = ( (fcn_info%dtype       == dQ_DIFF)                  .and. &
                          (elem_info%idomain_g  == fcn_info%seed%idomain_g)  .and. &
                          (elem_info%ielement_g == fcn_info%seed%ielement_g) .and. &
                          (itime                == fcn_info%seed%itime) ) 


        if ( differentiate ) then 
            ! Loop through the terms in q_out, seed appropriate derivatives to ONE
            do iterm = 1,nterms_s
                ! For the given term, seed its appropriate derivative
                set_deriv = (ifield - 1)*nterms_s + iterm
                q_out(iterm)%xp_ad_(set_deriv) = ONE
            end do
        else 
            ! Loop through the terms in q_out. Set all derivatives to ZERO
            do iterm = 1,nterms_s
                q_out(iterm)%xp_ad_ = ZERO 
            end do
        end if

    end function differentiate_modal_coefficients
    !********************************************************************************






    !>  Return the br2_face or br2_vol matrix with derivatives
    !!
    !!  @author Matteo Ugolotti
    !!  @date   1/31/2019
    !!
    !--------------------------------------------------------------------------------
    function differentiate_br2(mesh,elem_info,iface,fcn_info,source,exception) result(br2_matrix)
        type(mesh_t),                       intent(in)  :: mesh
        type(element_info_t),               intent(in)  :: elem_info
        integer(ik),                        intent(in)  :: iface
        type(function_info_t),              intent(in)  :: fcn_info
        character(*),                       intent(in)  :: source
        logical,                optional,   intent(in)  :: exception
    
        type(AD_D), allocatable, dimension(:,:)     :: br2_matrix
        real(rk),   allocatable, dimension(:,:)     :: br2_rk
        real(rk),   allocatable, dimension(:,:,:,:) :: br2_dx

        type(element_info_t)    :: source_info
        integer(ik)             :: idomain_l_n, ielement_l_n, idomain_g_n, ielement_g_n,    &
                                   iface_n, iproc_n, nderiv, br2_size1, br2_size2, nnodes_r,&
                                   ierr, irow, icol, idir, inode, deriv_index 
        logical                 :: differentiate_me, local_neighbor, remote_neighbor


        associate( idom_l  => elem_info%idomain_l,   &
                   ielem_l => elem_info%ielement_l  )

        !
        ! If 'face exterior' check if the neighbor is local or remote and assign source element
        !
        if (source == 'face exterior') then

            idomain_l_n  = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_domain_l
            ielement_l_n = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_element_l
            idomain_g_n  = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_domain_g
            ielement_g_n = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_element_g
            iface_n      = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_face
            iproc_n      = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_proc


            local_neighbor  = (iproc_n == IRANK)
            remote_neighbor = (iproc_n /= IRANK)

            source_info%idomain_l  = idomain_l_n
            source_info%idomain_g  = idomain_g_n
            source_info%ielement_l = ielement_l_n
            source_info%ielement_g = ielement_g_n
        
        else

            source_info = elem_info

        end if


        !
        ! If the face is a CHIMERA face, set source_info to zero.
        ! This is meant to handle the approximation on chimera faces
        ! for the BR2 scheme. We assumed that there exists a reflected
        ! element on the chimera surface and we use the br2_face of the
        ! interior element. However, we want to avoid that the BR2_face
        ! is differentiated wrt to interior in this case.
        !
        if (present(exception)) then
            if (exception)  then
                source_info%idomain_l  = 0 
                source_info%idomain_g  = 0 
                source_info%ielement_l = 0 
                source_info%ielement_g = 0 
            end if
        end if


        ! Get real BR2 values
        if (source == 'face interior') then
            br2_rk = mesh%domain(idom_l)%faces(ielem_l,iface)%br2_face
        else if (source == 'face exterior' .and. local_neighbor) then
            br2_rk = mesh%domain(idomain_l_n)%faces(ielement_l_n,iface_n)%br2_face
        else if (source == 'face exterior' .and. remote_neighbor) then
            br2_rk = mesh%domain(idom_l)%faces(ielem_l,iface)%neighbor_br2_face
        else if (source == 'element') then
            br2_rk = mesh%domain(idom_l)%faces(ielem_l,iface)%br2_vol
        else
            call chidg_signal(FATAL,"mod_differentiate differentiate_br2: Invalid value for 'source'. Options are 'face', 'element'")
        end if


        ! Get number of derivatives and interpolation nodes to initialize the AD_D vector
        nderiv        = get_nderiv(fcn_info)
        br2_size1     = size(br2_rk,1)
        br2_size2     = size(br2_rk,2)
        nnodes_r      = fcn_info%seed%nnodes_r


        ! Allocate resultant vector
        allocate(br2_matrix(br2_size1,br2_size2), stat=ierr)
        if (ierr/=0) call AllocationError
        do irow = 1,br2_size1
            do icol = 1,br2_size2
                br2_matrix(irow,icol) = AD_D(nderiv)
            end do
        end do
        

        ! Assign real values of the br2_matrix
        br2_matrix = br2_rk
        differentiate_me = ( (source_info%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (source_info%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (fcn_info%dtype         == dX_DIFF) )

        if (differentiate_me) then
            
            !
            ! If differentiate, select derivatives
            !
            if (source == 'face interior') then
                br2_dx = mesh%domain(idom_l)%faces(ielem_l,iface)%dbr2_f_dx
            else if (source == 'face exterior' .and. local_neighbor) then
                br2_dx = mesh%domain(idomain_l_n)%faces(ielement_l_n,iface_n)%dbr2_f_dx
            else if (source == 'face exterior' .and. remote_neighbor) then
                br2_dx = mesh%domain(idom_l)%faces(ielem_l,iface)%neighbor_dbr2_f_dx
            else if (source == 'element') then
                br2_dx = mesh%domain(idom_l)%faces(ielem_l,iface)%dbr2_v_dx
            else
                call chidg_signal(FATAL,"mod_differentiate differentiate_br2: Invalid value for 'source'. Options are 'face', 'element'")
            end if
        
        
            do irow = 1,br2_size1
                do icol = 1,br2_size2
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            
                            deriv_index = (idir-1)*nnodes_r + inode
                            br2_matrix(irow,icol)%xp_ad_(deriv_index) = br2_dx(irow,icol,inode,idir)

                        end do !inode
                    end do !idir
                end do !icol
            end do !irow
       
        end if


        end associate

    end function differentiate_br2
    !********************************************************************************






    !>  Return the inverse jacobian with derivatives
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/6/2018
    !!
    !--------------------------------------------------------------------------------
    function differentiate_jinv(mesh,elem_info,iface,fcn_info,source) result(jinv)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        integer(ik),            intent(in)  :: iface
        type(function_info_t),  intent(in)  :: fcn_info
        character(*),           intent(in)  :: source
    
        type(AD_D), allocatable, dimension(:)   :: jinv
        
        real(rk),   allocatable, dimension(:)       :: jinv_rk, djinv_d1, djinv_d2, &
                                                       djinv_d3
        real(rk),   allocatable, dimension(:,:,:)   :: jinv_dx
        integer(ik)             :: nderiv, nnodes_interp, nnodes_r, &
                                   inode, d1_start, d1_end,         &
                                   d2_start, d2_end, d3_start,      &
                                   d3_end, ierr
        logical                 :: mismatch, differentiate_me

        associate( idom_l  => elem_info%idomain_l,   &
                   ielem_l => elem_info%ielement_l  )

        ! Get real jinv values
        if (source == 'face') then
            jinv_rk = mesh%domain(idom_l)%faces(ielem_l,iface)%jinv
        else if (source == 'element') then
            jinv_rk = mesh%domain(idom_l)%elems(ielem_l)%jinv
        else
            call chidg_signal(FATAL,"mod_differentiate differentiate_jinv: Invalid value for 'source'. Options are 'face', 'element'")
        end if


        ! Get number of derivatives and interpolation nodes to initialize the AD_D vector
        nderiv        = get_nderiv(fcn_info)
        nnodes_interp = size(jinv_rk)
        nnodes_r      = fcn_info%seed%nnodes_r

        ! Allocate resultant vector
        allocate(jinv(nnodes_interp), stat=ierr)
        if (ierr/=0) call AllocationError
        do inode = 1,nnodes_interp
            jinv(inode) = AD_D(nderiv)
        end do
        
        ! Assign real values of the jinv
        jinv = jinv_rk
        
        differentiate_me = ( (elem_info%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (elem_info%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (dX_DIFF              == fcn_info%dtype))

        if (differentiate_me) then

            do inode = 1,nnodes_interp
        
                ! Get real jinv values
                if (source == 'face') then
                    jinv_dx = mesh%domain(idom_l)%faces(ielem_l,iface)%djinv_dx
                else if (source == 'element') then
                    jinv_dx = mesh%domain(idom_l)%elems(ielem_l)%djinv_dx
                else
                    call chidg_signal(FATAL,"mod_differentiate differentiate_jinv: Invalid value for 'source'. Options are 'face', 'element'")
                end if
        
                ! retrieve derivatives wrt grid nodes in direction 1,2 and 3 for the specific 
                ! interpolation node                
                djinv_d1 = jinv_dx(inode,:,1)
                djinv_d2 = jinv_dx(inode,:,2)
                djinv_d3 = jinv_dx(inode,:,3)
                
                ! Define indeces to start and end the vector allocation
                ! Remember the order: x-derivatives, y-derivatives and z-derivatives
                d1_start = 1
                d1_end   = nnodes_r
                d2_start = nnodes_r +1
                d2_end   = 2*nnodes_r
                d3_start = 2*nnodes_r + 1
                d3_end   = 3*nnodes_r

                ! Check if the number of reference grid nodes is correct
                mismatch = (size(djinv_d1) /= nnodes_r .or. &
                            size(djinv_d2) /= nnodes_r .or. &
                            size(djinv_d3) /= nnodes_r      )

                if (mismatch) then 
                    call chidg_signal(FATAL,"differentiate_jinv: mismatch between djinv_dx reference nodes and element grid nodes .")
                end if

                ! Store each directional derivaties in the correct order.
                jinv(inode)%xp_ad_(d1_start:d1_end) = djinv_d1
                jinv(inode)%xp_ad_(d2_start:d2_end) = djinv_d2
                jinv(inode)%xp_ad_(d3_start:d3_end) = djinv_d3
        
            end do

        end if


        end associate

    end function differentiate_jinv 
    !********************************************************************************









    !>  Return the face normal with derivatives
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/9/2019
    !!
    !!  Reconstructed after bug
    !!
    !!
    !--------------------------------------------------------------------------------
    function differentiate_normal(mesh,elem_info,fcn_info,iface,direction) result(norm_gq)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: iface
        integer(ik),            intent(in)  :: direction
    
        type(AD_D), allocatable, dimension(:)       :: norm_gq
        real(rk),   allocatable, dimension(:)       :: norm_rk
        real(rk),   allocatable, dimension(:,:,:,:) :: norm_dx

        type(element_info_t)    :: source_info
        integer(ik)             :: nderiv, interp_nodes, nnodes_r,            &
                                   d1_start, d1_end, d2_start, d2_end, d3_start, d3_end,& 
                                   ierr, irow, icol, idir, inode, deriv_index 
        logical                 :: diff_interior, differentiate, mismatch


        associate( idom  => elem_info%idomain_l,   &
                   ielem => elem_info%ielement_l )


        ! Get real normal values, of the current face
        norm_rk = mesh%domain(idom)%faces(ielem,iface)%norm(:,direction)


        ! Get number of derivatives and interpolation nodes to initialize the AD_D vector
        nderiv        = get_nderiv(fcn_info)
        interp_nodes  = size(norm_rk)
        nnodes_r      = fcn_info%seed%nnodes_r


        ! Allocate resultant vector
        allocate(norm_gq(interp_nodes), stat=ierr)
        if (ierr/=0) call AllocationError
        do irow = 1,interp_nodes
            norm_gq(irow) = AD_D(nderiv)
        end do
        
        ! Assign real value
        norm_gq = norm_rk 
        
        differentiate = (fcn_info%dtype == dX_DIFF)
        diff_interior = ( (elem_info%idomain_g  == fcn_info%seed%idomain_g) .and. &
                          (elem_info%ielement_g == fcn_info%seed%ielement_g) )

        if (differentiate .and. diff_interior) then
            
            ! Retrieve derivatives wrt to interior
            norm_dx = mesh%domain(idom)%faces(ielem,iface)%dnorm_dx
        
            do inode = 1,interp_nodes

                ! Define indeces to start and end the vector allocation
                ! Remember the order: x-derivatives, y-derivatives and z-derivatives
                d1_start = 1
                d1_end   = nnodes_r
                d2_start = nnodes_r +1
                d2_end   = 2*nnodes_r
                d3_start = 2*nnodes_r + 1
                d3_end   = 3*nnodes_r

                ! Check if the number of reference grid nodes is correct
                mismatch = (size(norm_dx(inode,direction,:,1)) /= nnodes_r .or. &
                            size(norm_dx(inode,direction,:,2)) /= nnodes_r .or. &
                            size(norm_dx(inode,direction,:,3)) /= nnodes_r      )

                if (mismatch) then 
                    call chidg_signal(FATAL,"differentiate_normal: mismatch between dnorm_dx reference nodes and element grid nodes.")
                end if

                ! Store each directional derivaties in the correct order.
                norm_gq(inode)%xp_ad_(d1_start:d1_end) = norm_dx(inode,direction,:,1)
                norm_gq(inode)%xp_ad_(d2_start:d2_end) = norm_dx(inode,direction,:,2)
                norm_gq(inode)%xp_ad_(d3_start:d3_end) = norm_dx(inode,direction,:,3)
            
            end do
       
        end if

        end associate

    end function differentiate_normal
    !********************************************************************************




    !>  Return the face normal with derivatives
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/24/2020
    !!
    !!  Reconstructed after bug
    !!
    !--------------------------------------------------------------------------------
    function differentiate_normal_ale(mesh,elem_info,fcn_info,iface,direction) result(norm_ale_gq)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: iface
        integer(ik),            intent(in)  :: direction
    
        type(AD_D), allocatable, dimension(:)       :: norm_ale_gq
        real(rk),   allocatable, dimension(:)       :: norm_ale_rk
        real(rk),   allocatable, dimension(:,:,:,:) :: norm_ale_dx

        type(element_info_t)    :: source_info
        integer(ik)             :: nderiv, interp_nodes, nnodes_r,            &
                                   d1_start, d1_end, d2_start, d2_end, d3_start, d3_end,& 
                                   ierr, irow, icol, idir, inode, deriv_index 
        logical                 :: diff_interior, differentiate, mismatch


        associate( idom  => elem_info%idomain_l,   &
                   ielem => elem_info%ielement_l )


        ! Get real normal values, of the current face
        norm_ale_rk = mesh%domain(idom)%faces(ielem,iface)%norm_def(:,direction)


        ! Get number of derivatives and interpolation nodes to initialize the AD_D vector
        nderiv        = get_nderiv(fcn_info)
        interp_nodes  = size(norm_ale_rk)
        nnodes_r      = fcn_info%seed%nnodes_r


        ! Allocate resultant vector
        allocate(norm_ale_gq(interp_nodes), stat=ierr)
        if (ierr/=0) call AllocationError
        do irow = 1,interp_nodes
            norm_ale_gq(irow) = AD_D(nderiv)
        end do
        
        ! Assign real value
        norm_ale_gq = norm_ale_rk 
        
        differentiate = (fcn_info%dtype == dX_DIFF)
        diff_interior = ( (elem_info%idomain_g  == fcn_info%seed%idomain_g) .and. &
                          (elem_info%ielement_g == fcn_info%seed%ielement_g) )

        if (differentiate .and. diff_interior) then
            
            ! Retrieve derivatives wrt to interior
            norm_ale_dx = mesh%domain(idom)%faces(ielem,iface)%dnorm_ale_dx
        
            do inode = 1,interp_nodes

                ! Define indeces to start and end the vector allocation
                ! Remember the order: x-derivatives, y-derivatives and z-derivatives
                d1_start = 1
                d1_end   = nnodes_r
                d2_start = nnodes_r +1
                d2_end   = 2*nnodes_r
                d3_start = 2*nnodes_r + 1
                d3_end   = 3*nnodes_r

                ! Check if the number of reference grid nodes is correct
                mismatch = (size(norm_ale_dx(inode,direction,:,1)) /= nnodes_r .or. &
                            size(norm_ale_dx(inode,direction,:,2)) /= nnodes_r .or. &
                            size(norm_ale_dx(inode,direction,:,3)) /= nnodes_r      )

                if (mismatch) then 
                    call chidg_signal(FATAL,"differentiate_normal: mismatch between dnorm_dx reference nodes and element grid nodes.")
                end if

                ! Store each directional derivaties in the correct order.
                norm_ale_gq(inode)%xp_ad_(d1_start:d1_end) = norm_ale_dx(inode,direction,:,1)
                norm_ale_gq(inode)%xp_ad_(d2_start:d2_end) = norm_ale_dx(inode,direction,:,2)
                norm_ale_gq(inode)%xp_ad_(d3_start:d3_end) = norm_ale_dx(inode,direction,:,3)
            
            end do
       
        end if

        end associate

    end function differentiate_normal_ale
    !********************************************************************************





    !>  Return the face unit normals with derivatives
    !!
    !!  Note: this procedure is not explicitly tested because it relies on differentiate_normal
    !!        and DNAD procedures that are tested.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/6/2018
    !!
    !--------------------------------------------------------------------------------
    function differentiate_unit_normal(mesh,elem_info,fcn_info,iface,direction) result(unorm_gq)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: iface
        integer(ik),            intent(in)  :: direction

        type(AD_D), dimension(:), allocatable :: unorm_gq

        type(AD_D), dimension(:), allocatable :: norm1_gq, norm2_gq, norm3_gq, norm_mag

        norm1_gq = differentiate_normal(mesh,elem_info,fcn_info,iface,1)
        norm2_gq = differentiate_normal(mesh,elem_info,fcn_info,iface,2)
        norm3_gq = differentiate_normal(mesh,elem_info,fcn_info,iface,3)
        
        norm_mag = sqrt(norm1_gq**TWO + norm2_gq**TWO + norm3_gq**TWO)

        select case (direction)
            case(1)
                unorm_gq = norm1_gq/norm_mag
            case(2)
                unorm_gq = norm2_gq/norm_mag
            case(3)
                unorm_gq = norm3_gq/norm_mag
            case default
                call chidg_signal_one(FATAL,"differentiate_unit_normal: Invalid direction for selecting coordinate.",direction)
        end select

    end function differentiate_unit_normal
    !********************************************************************************




    !>  Return the face ALE unit normals with derivatives
    !!
    !!  Note: this procedure is not explicitly tested because it relies on differentiate_normal
    !!        and DNAD procedures that are tested.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/27/2020
    !!
    !--------------------------------------------------------------------------------
    function differentiate_unit_normal_ale(mesh,elem_info,fcn_info,iface,direction) result(unorm_ale_gq)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: iface
        integer(ik),            intent(in)  :: direction

        type(AD_D), dimension(:), allocatable :: unorm_ale_gq

        type(AD_D), dimension(:), allocatable :: norm1_ale_gq, norm2_ale_gq, norm3_ale_gq, norm_mag

        norm1_ale_gq = differentiate_normal_ale(mesh,elem_info,fcn_info,iface,1)
        norm2_ale_gq = differentiate_normal_ale(mesh,elem_info,fcn_info,iface,2)
        norm3_ale_gq = differentiate_normal_ale(mesh,elem_info,fcn_info,iface,3)
        
        norm_mag = sqrt(norm1_ale_gq**TWO + norm2_ale_gq**TWO + norm3_ale_gq**TWO)

        select case (direction)
            case(1)
                unorm_ale_gq = norm1_ale_gq/norm_mag
            case(2)
                unorm_ale_gq = norm2_ale_gq/norm_mag
            case(3)
                unorm_ale_gq = norm3_ale_gq/norm_mag
            case default
                call chidg_signal_one(FATAL,"differentiate_unit_normal_ale: Invalid direction for selecting coordinate.",direction)
        end select

    end function differentiate_unit_normal_ale
    !********************************************************************************





    !>  Return the coordinates at GQ nodes  with derivatives
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/6/2018
    !!
    !--------------------------------------------------------------------------------
    function differentiate_coordinate(mesh,elem_info,fcn_info,iface,direction,source) result(coords)
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: iface
        character(*),           intent(in)  :: direction
        character(*),           intent(in)  :: source

        type(AD_D),     allocatable, dimension(:)   :: coords

        character(:),   allocatable                 :: user_msg
        real(rk),       allocatable, dimension(:)   :: gq_1, gq_2, gq_3
        logical                                     :: differentiate_me
        integer(ik)                                 :: nderiv, nnodes_interp, inode, ierr

        associate( idom_l   => elem_info%idomain_l,     &
                   idom_g   => elem_info%idomain_g,     &
                   ielem_l  => elem_info%ielement_l,    &
                   ielem_g  => elem_info%ielement_g )

        ! Get coordinates
        if ( (source == 'boundary') .or. (source == 'face interior') .or. (source == 'face exterior') ) then
            gq_1 = mesh%domain(idom_l)%faces(ielem_l,iface)%interp_coords_def(:,1)
            gq_2 = mesh%domain(idom_l)%faces(ielem_l,iface)%interp_coords_def(:,2)
            gq_3 = mesh%domain(idom_l)%faces(ielem_l,iface)%interp_coords_def(:,3)
        else if ( (source == 'volume') .or. (source == 'element') ) then
            gq_1 = mesh%domain(idom_l)%elems(ielem_l)%interp_coords_def(:,1)
            gq_2 = mesh%domain(idom_l)%elems(ielem_l)%interp_coords_def(:,2)
            gq_3 = mesh%domain(idom_l)%elems(ielem_l)%interp_coords_def(:,3)
        else
            user_msg = "differentiate_coordinate: Invalid source for returning coordinate. Options are 'boundary' and 'volume'."
            call chidg_signal_one(FATAL,user_msg,source)
        end if


        ! Initialize number of derivatives
        nderiv          = get_nderiv(fcn_info)
        nnodes_interp   = size(gq_1)
        allocate(coords(nnodes_interp), stat=ierr)
        if (ierr/=0) call AllocationError
        do inode = 1,nnodes_interp
            coords(inode) = AD_D(nderiv)
        end do


        ! Define coordinate to return with derivatives set to ZERO (ok for standard linearization).
        select case (direction)
            case ('1')
                coords = gq_1
            case ('2')
                coords = gq_2
            case ('3')
                coords = gq_3
            case default
                call chidg_signal_one(FATAL,"differentiate_coordinate: Invalid string for selecting coordinate.",direction)
        end select


        ! If grid-node sensitivities is required overwrite the derivatives
        differentiate_me = ( (idom_g  == fcn_info%seed%idomain_g ) .and. &
                             (ielem_g == fcn_info%seed%ielement_g) .and. &
                              fcn_info%dtype == dX_DIFF )
        
        if (differentiate_me) then
            do inode = 1,nnodes_interp
                coords(inode)%xp_ad_ = mesh%domain(idom_l)%elems(ielem_l)%basis_c%get_coord_derivatives(source,direction,inode,iface)
            end do
        end if

        end associate

    end function differentiate_coordinate
    !********************************************************************************






!    !>  Return the element size with derivatives 
!    !!  DOUBT: not sure of this
!    !!  TO BE TESTED
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   8/6/2018
!    !!
!    !--------------------------------------------------------------------------------
!    function differentiate_element_size(mesh,face_info,face_type,fcn_info,source) result(h)
!        type(mesh_t),           intent(in)  :: mesh
!        type(face_info_t),      intent(in)  :: face_info
!        integer(ik),            intent(in)  :: face_type
!        type(function_info_t),  intent(in)  :: fcn_info
!        character(*),           intent(in)  :: source
!
!
!        integer(ik) :: ineighbor_domain_l, ineighbor_element_l, h_index(3,2),   &
!                       nderiv, nnodes_r, idir, imin_l, imax_l,                  &
!                       iminimum, imaximum, ineighbor_domain_g, ineighbor_element_g
!        type(AD_D)  :: h(3)
!        logical     :: proc_local, chimera_face, diff_me
!        real(rk)    :: h_rk(3)
!
!
!        associate( idom_l   =>  face_info%idomain_l,    &
!                   idom_g   =>  face_info%idomain_g,    &
!                   ielem_l  =>  face_info%ielement_l,   &
!                   ielem_g  =>  face_info%ielement_g,   &
!                   iface    =>  face_info%iface )
!
!        ! Get number of derivatives and interpolation nodes to initialize the AD_D vector
!        nderiv         = get_nderiv(fcn_info)
!        nnodes_r       = fcn_info%seed%nnodes_r
!
!        ! Allocate resultant vectors
!        do idir = 1,3
!            h(idir) = AD_D(nderiv)
!        end do
!
!        ! Find h based on face source
!        if (source == 'interior') then
!
!            h_rk    = mesh%domain(idom_l)%elems(ielem_l)%h
!            h_index = mesh%domain(idom_l)%elems(ielem_l)%h_index
!            diff_me = ((idom_g  == fcn_info%seed%idomain_g ) .and. &
!                       (ielem_g == fcn_info%seed%ielement_g) .and. &
!                        fcn_info%dtype == dX_DIFF )
!        
!        else if (source == 'exterior') then
!
!            ! If Chimera face, use interior element size. APPROXIMATION
!            ! This is because the size of the bounding box for the
!            ! exterior chimera element might be of the same size of the current
!            ! element.
!            ! For dx-sensitivities, set derivatives to zero (ASSUMPTION)
!            chimera_face = (face_type == CHIMERA)
!            if (chimera_face) then
!
!                h_rk    = mesh%domain(idom_l)%elems(ielem_l)%h
!                diff_me = .false.
!
!            ! If conforming face, check for processor status of neighbor.
!            else
!
!                proc_local = (mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_proc  ==  IRANK)
!                if (proc_local) then
!
!                    ineighbor_domain_l  = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_domain_l
!                    ineighbor_domain_g  = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_domain_g
!                    ineighbor_element_l = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_element_l
!                    ineighbor_element_g = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_element_g
!                    h_rk    = mesh%domain(ineighbor_domain_l)%elems(ineighbor_element_l)%h
!                    h_index = mesh%domain(ineighbor_domain_l)%elems(ineighbor_element_l)%h_index
!                    diff_me = ((ineighbor_domain_g  == fcn_info%seed%idomain_g ) .and. &
!                               (ineighbor_element_g == fcn_info%seed%ielement_g) .and. &
!                                fcn_info%dtype == dX_DIFF )
!
!                else
!
!                    ineighbor_domain_g  = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_domain_g
!                    ineighbor_element_g = mesh%domain(idom_l)%faces(ielem_l,iface)%ineighbor_element_g
!                    h_rk    = mesh%domain(idom_l)%faces(ielem_l,iface)%neighbor_h
!                    h_index = mesh%domain(idom_l)%faces(ielem_l,iface)%neighbor_h_index
!                    diff_me = ((ineighbor_domain_g  == fcn_info%seed%idomain_g ) .and. &
!                               (ineighbor_element_g == fcn_info%seed%ielement_g) .and. &
!                                fcn_info%dtype == dX_DIFF )
!
!                end if
!
!            end if
!
!
!        else
!            call chidg_signal(FATAL,"differentiate_element_size(source): Invalid value for 'source'. Options are 'interior', 'exterior'")
!        end if
!
!
!        ! Set the real values of the volume
!        h = h_rk
!
!        if (diff_me) then
!
!            do idir = 1,3
!                ! Get local indeces of min and max
!                imin_l = h_index(idir,1)
!                imax_l = h_index(idir,2)
!                ! Define the location in the derivative vector to locate the 
!                ! derivatives (=1 for max, =-1 for min)
!                iminimum = (idir-1)*nnodes_r + imin_l
!                imaximum = (idir-1)*nnodes_r + imax_l
!                if (imin_l /= 0 ) h(idir)%xp_ad_(iminimum) = -ONE
!                if (imax_l /= 0 ) h(idir)%xp_ad_(imaximum) = ONE
!            end do
!
!        end if
!
!        end associate
!
!    end function differentiate_element_size
!    !********************************************************************************







    !>  Set derivatives for element interpolation
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/29/2018
    !!
    !-----------------------------------------------------------------------------------------
    function differentiate_element_interpolator(interpolation_type,mesh,elem_info,fcn_info,itime) result(interpolator)
        character(*),           intent(in)  :: interpolation_type
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: elem_info
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: itime
        
        type(AD_D), allocatable :: interpolator(:,:)

        integer(ik)             :: nderiv, nterms_s, nnodes_i, nnodes_r, igq, iterm, inode, idir, &
                                   deriv_index, ierr
        real(rk), allocatable   :: interpolator_rk(:,:), interpolator_dx(:,:,:,:)
        logical                 :: differentiate_me
        
        associate( idomain => elem_info%idomain_l, ielement => elem_info%ielement_l )
         
        ! Access grad and dgrad_dx matrices for give indeces
        select case(interpolation_type)
            case('value')
                interpolator_rk = mesh%domain(idomain)%elems(ielement)%basis_s%interpolator_element('Value') 
            case('grad1')
                interpolator_rk = mesh%domain(idomain)%elems(ielement)%grad1
                interpolator_dx = mesh%domain(idomain)%elems(ielement)%dgrad1_dx
            case('grad2')
                interpolator_rk = mesh%domain(idomain)%elems(ielement)%grad2
                interpolator_dx = mesh%domain(idomain)%elems(ielement)%dgrad2_dx
            case('grad3')
                interpolator_rk = mesh%domain(idomain)%elems(ielement)%grad3
                interpolator_dx = mesh%domain(idomain)%elems(ielement)%dgrad3_dx
            case default
                call chidg_signal(FATAL,"differentiate_element_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
        end select

        ! Compute number of derivatives and reference nodes 
        nderiv   = get_nderiv(fcn_info)
        nterms_s = size(interpolator_rk,2)
        nnodes_i = size(interpolator_rk,1)
        nnodes_r = fcn_info%seed%nnodes_r 

        ! Allocate the interpolator
        if (allocated(interpolator)) deallocate(interpolator)
        allocate(interpolator(nnodes_i,nterms_s), stat=ierr)
        if (ierr /=0) call AllocationError

        ! For each entry of the interpolator define the number of derivatives in the AD_D
        do igq = 1,nnodes_i
            do iterm = 1,nterms_s
                interpolator(igq,iterm) = AD_D(nderiv)
            end do
        end do
        
        ! Set the interpolator%x_ad_ equal to the real interpolator
        interpolator = interpolator_rk

        differentiate_me = ( (elem_info%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (elem_info%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (itime                == fcn_info%seed%itime) )

        ! Specify the derivatives if node grid sensitivities are required.
        if (fcn_info%dtype == dX_DIFF .and. (interpolation_type /= 'value') .and. differentiate_me) then

            do igq = 1,nnodes_i
                do iterm = 1,nterms_s
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            
                            deriv_index = (idir-1)*nnodes_r + inode
                            interpolator(igq,iterm)%xp_ad_(deriv_index) = interpolator_dx(igq,iterm,inode,idir)

                        end do !inode
                    end do !idir
                end do !iterm
            end do !igq

        end if

        ! Sanity check: this subroutine should not be called if fcn_info%dtype == 'dQ' or 'NO'
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF) then
            call chidg_signal(FATAL,"differentiate_element_interpolator: this function should not be called &
                                     for dtype = 'dQ' or 'none'. Implementation error.") 
        end if


        end associate


    end function differentiate_element_interpolator
    !*****************************************************************************************






    !>  Set derivatives for interior face interpolation
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/29/2018
    !!
    !-----------------------------------------------------------------------------------------
    function differentiate_face_interior_interpolator(interpolation_type,mesh,source_elem,source_iface,donor_elem,donor_iface,fcn_info,itime) result(interpolator)
        character(*),           intent(in)  :: interpolation_type
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: source_elem
        integer(ik),            intent(in)  :: source_iface
        type(element_info_t),   intent(in)  :: donor_elem
        integer(ik),            intent(in)  :: donor_iface
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: itime
        
        type(AD_D), allocatable :: interpolator(:,:)

        integer(ik)             :: nderiv, nterms_s, nnodes_f, nnodes_r, igq, iterm, inode, idir, &
                                   ierr, deriv_index
        real(rk), allocatable   :: interpolator_rk(:,:), interpolator_dx(:,:,:,:)
        logical                 :: differentiate_me
         
        associate( idom => source_elem%idomain_l, ielem => source_elem%ielement_l, iface => source_iface )
        

        ! Access grad and dgrad_dx matrices for give indeces
        select case(interpolation_type)
            case('value')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%basis_s%interpolator_face('Value',iface) 
            case('grad1')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%grad1
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%dgrad1_dx
            case('grad2')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%grad2
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%dgrad2_dx
            case('grad3')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%grad3
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%dgrad3_dx
            case default
                call chidg_signal(FATAL,"differentiate_face_interior_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
        end select

        ! Compute number of derivatives and reference nodes 
        nderiv   = get_nderiv(fcn_info)
        nterms_s = size(interpolator_rk,2)
        nnodes_f = size(interpolator_rk,1)
        nnodes_r = fcn_info%seed%nnodes_r 

        ! Allocate the interpolator
        if (allocated(interpolator)) deallocate(interpolator)
        allocate(interpolator(nnodes_f,nterms_s), stat=ierr)
        if (ierr /=0) call AllocationError

        ! For each entry of the interpolator define the number of derivatives in the AD_D
        do igq = 1,nnodes_f
            do iterm = 1,nterms_s
                interpolator(igq,iterm) = AD_D(nderiv)
            end do
        end do
        
        ! Set the interpolator%x_ad_ equal to the real interpolator
        interpolator = interpolator_rk

        ! Interpolation_source = ME means that donor_elem is equal to worker%elem_info.
        ! This is defined by mod_interpolate%get_face_interpolation_info.
        differentiate_me = ( (donor_elem%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (donor_elem%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (itime                 == fcn_info%seed%itime) )

        ! Specify the derivatives if node grid sensitivities are require.
        if (fcn_info%dtype == dX_DIFF .and. (interpolation_type /= 'value') .and. differentiate_me) then

            do igq = 1,nnodes_f
                do iterm = 1,nterms_s
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            
                            deriv_index = (idir-1)*nnodes_r + inode
                            interpolator(igq,iterm)%xp_ad_(deriv_index) = interpolator_dx(igq,iterm,inode,idir)

                        end do !inode
                    end do !idir
                end do !iterm
            end do !igq

        end if
    
        ! Sanity check: this subroutine should not be called if fcn_info%dtype == 'dQ' or 'NO'
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF) then
            call chidg_signal(FATAL,"differentiate_face_interior_interpolator: this function should not be called &
                                     for dtype = 'dQ' or 'NO'. Implementation error.") 
        end if

        end associate


    end function differentiate_face_interior_interpolator
    !*****************************************************************************************






    !>  Set derivatives for exterior parallel face interpolation
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/29/2018
    !!
    !-----------------------------------------------------------------------------------------
    function differentiate_face_parallel_interpolator(interpolation_type,mesh,source_elem,source_iface,donor_elem,donor_iface,fcn_info,itime) result(interpolator)
        character(*),           intent(in)  :: interpolation_type
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: source_elem
        integer(ik),            intent(in)  :: source_iface
        type(element_info_t),   intent(in)  :: donor_elem
        integer(ik),            intent(in)  :: donor_iface
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: itime
        
        type(AD_D), allocatable :: interpolator(:,:)

        integer(ik)                 :: nderiv, nterms_s, nnodes_f, nnodes_r, igq, iterm, inode, idir, &
                                       deriv_index, ierr
        real(rk),       allocatable :: interpolator_rk(:,:), interpolator_dx(:,:,:,:)
        logical                     :: differentiate_me

        
        associate( idom => source_elem%idomain_l, ielem => source_elem%ielement_l, iface => source_iface )
        
        
        ! Access grad and dgrad_dx matrices for give indeces
        select case(interpolation_type)
            case('value')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%basis_s%interpolator_face('Value',donor_iface)
            case('grad1')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad1
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%neighbor_dgrad1_dx
            case('grad2')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad2
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%neighbor_dgrad2_dx
            case('grad3')
                interpolator_rk = mesh%domain(idom)%faces(ielem,iface)%neighbor_grad3
                interpolator_dx = mesh%domain(idom)%faces(ielem,iface)%neighbor_dgrad3_dx
            case default
                call chidg_signal(FATAL,"differentiate_parallel_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
        end select


        ! Compute number of derivatives and reference nodes 
        nderiv   = get_nderiv(fcn_info)
        nterms_s = size(interpolator_rk,2)
        nnodes_f = size(interpolator_rk,1)
        nnodes_r = fcn_info%seed%nnodes_r 


        ! Allocate the interpolator
        if (allocated(interpolator)) deallocate(interpolator)
        allocate(interpolator(nnodes_f,nterms_s), stat=ierr)
        if (ierr /=0) call AllocationError

        ! For each entry of the interpolator define the number of derivatives in the AD_D
        do igq = 1,nnodes_f
            do iterm = 1,nterms_s
                interpolator(igq,iterm) = AD_D(nderiv)
            end do
        end do
        
        ! Set the interpolator%x_ad_ equal to the real interpolator
        interpolator = interpolator_rk
        differentiate_me = ( (donor_elem%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (donor_elem%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (itime                 == fcn_info%seed%itime) )

        ! Specify the derivatives if node grid sensitivities are require.
        if (fcn_info%dtype == dX_DIFF .and. (interpolation_type /= 'value') .and. differentiate_me) then

            do igq = 1,nnodes_f
                do iterm = 1,nterms_s
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            deriv_index = (idir-1)*nnodes_r + inode
                            interpolator(igq,iterm)%xp_ad_(deriv_index) = interpolator_dx(igq,iterm,inode,idir)
                        end do !inode
                    end do !idir
                end do !iterm
            end do !igq

        end if


        ! Sanity check: this subroutine should not be called if fcn_info%dtype == 'dQ' or 'NO'
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF) then
            call chidg_signal(FATAL,"differentiate_parallel_interpolator: this function should not be called &
                                     for dtype = 'dQ' or 'NO'. Implementation error.") 
        end if

        end associate

    end function differentiate_face_parallel_interpolator
    !*****************************************************************************************








    !>  Set derivatives for exterior local face interpolation
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/29/2018
    !!
    !-----------------------------------------------------------------------------------------
    function differentiate_face_local_interpolator(interpolation_type,mesh,source_elem,source_iface,donor_elem,donor_iface,fcn_info,itime) result(interpolator)
        character(*),           intent(in)  :: interpolation_type
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: source_elem
        integer(ik),            intent(in)  :: source_iface
        type(element_info_t),   intent(in)  :: donor_elem
        integer(ik),            intent(in)  :: donor_iface
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: itime
        
        type(AD_D), allocatable :: interpolator(:,:)

        integer(ik)                 :: nderiv, nterms_s, nnodes_f, nnodes_r, igq, iterm, inode, idir, &
                                       deriv_index, inode_l, iindex, ierr
        real(rk),       allocatable :: interpolator_rk(:,:), interpolator_dx(:,:,:,:)
        logical                     :: differentiate_me
         
        associate( idom         => source_elem%idomain_l,   & 
                   ielem        => source_elem%ielement_l,  &
                   iface        => source_iface,       &
                   donor_idom   => donor_elem%idomain_l,    &
                   donor_ielem  => donor_elem%ielement_l) 

        ! Access grad and dgrad_dx matrices for give indeces of the local neighbor face
        select case(interpolation_type)
            case('value')
                interpolator_rk = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%basis_s%interpolator_face('Value',donor_iface)
            case('grad1')
                interpolator_rk = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%grad1
                interpolator_dx = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%dgrad1_dx
            case('grad2')
                interpolator_rk = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%grad2
                interpolator_dx = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%dgrad2_dx
            case('grad3')
                interpolator_rk = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%grad3
                interpolator_dx = mesh%domain(donor_idom)%faces(donor_ielem,donor_iface)%dgrad3_dx
            case default
                call chidg_signal(FATAL,"differentiate_local_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
        end select

        ! Compute number of derivatives and reference nodes 
        nderiv   = get_nderiv(fcn_info)
        nterms_s = size(interpolator_rk,2)
        nnodes_f = size(interpolator_rk,1)
        nnodes_r = fcn_info%seed%nnodes_r 

        ! Allocate the interpolator
        if (allocated(interpolator)) deallocate(interpolator)
        allocate(interpolator(nnodes_f,nterms_s), stat=ierr)
        if (ierr /=0) call AllocationError

        ! For each entry of the interpolator define the number of derivatives in the AD_D
        do igq = 1,nnodes_f
            do iterm = 1,nterms_s
                interpolator(igq,iterm) = AD_D(nderiv)
            end do
        end do
        
        ! Set the interpolator%x_ad_ equal to the real interpolator
        interpolator = interpolator_rk
        differentiate_me = ( (donor_elem%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (donor_elem%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (itime                 == fcn_info%seed%itime) )

        ! Specify the derivatives if node grid sensitivities are require.
        if (fcn_info%dtype == dX_DIFF .and. (interpolation_type /= 'value') .and. differentiate_me) then

            do igq = 1,nnodes_f
                do iterm = 1,nterms_s
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            deriv_index = (idir-1)*nnodes_r + inode
                            interpolator(igq,iterm)%xp_ad_(deriv_index) = interpolator_dx(igq,iterm,inode,idir)
                        end do !inode
                    end do !idir
                end do !iterm
            end do !igq

        end if

        ! Sanity check: this subroutine should not be called if fcn_info%dtype == 'dQ' or 'NO'
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF) then
            call chidg_signal(FATAL,"differentiate_face_local_interpolator: this function should not be called &
                                     for dtype = 'dQ' or 'NO'. Implementation error.") 
        end if

        end associate

    end function differentiate_face_local_interpolator
    !*****************************************************************************************






    !>  Set derivatives for exterior chimera face interpolation
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/29/2018
    !!
    !-----------------------------------------------------------------------------------------
    function differentiate_face_chimera_interpolator(interpolation_type,mesh,source_elem,source_iface,donor_elem,ChiID,idonor,fcn_info,itime) result(interpolator)
        character(*),           intent(in)  :: interpolation_type
        type(mesh_t),           intent(in)  :: mesh
        type(element_info_t),   intent(in)  :: source_elem
        integer(ik),            intent(in)  :: source_iface
        type(element_info_t),   intent(in)  :: donor_elem
        integer(ik),            intent(in)  :: ChiID
        integer(ik),            intent(in)  :: idonor
        type(function_info_t),  intent(in)  :: fcn_info
        integer(ik),            intent(in)  :: itime
        
        type(AD_D), allocatable :: interpolator(:,:)

        integer(ik)                 :: nderiv, nterms_s, nnodes_f, igq, iterm, ierr, &
                                       deriv_index, idir, inode, nnodes_r
        real(rk),       allocatable :: interpolator_rk(:,:),interpolator_dx(:,:,:,:)
        logical                     :: differentiate_me
         
        associate( idom => source_elem%idomain_l, ielem => source_elem%ielement_l, iface => source_iface )
        
        
        ! Access grad and dgrad_dx matrices for give indeces of the local neighbor face
        select case(interpolation_type)
            case('value')
                interpolator_rk = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%value
            case('grad1')
                interpolator_rk = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad1
                interpolator_dx = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%dgrad1_dx
            case('grad2')
                interpolator_rk = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad2
                interpolator_dx = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%dgrad2_dx
            case('grad3')
                interpolator_rk = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad3
                interpolator_dx = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%dgrad3_dx
            case default
                call chidg_signal(FATAL,"differentiate_chimera_interpolator: Invalid interpolation_type. Options are 'value', 'grad1', 'grad2', 'grad3'.")
        end select

        ! Compute number of derivatives and reference nodes 
        nderiv   = get_nderiv(fcn_info)
        nterms_s = size(interpolator_rk,2)
        nnodes_f = size(interpolator_rk,1)
        nnodes_r = fcn_info%seed%nnodes_r 
        
        ! Allocate the interpolator
        if (allocated(interpolator)) deallocate(interpolator)
        allocate(interpolator(nnodes_f,nterms_s), stat=ierr)
        if (ierr /=0) call AllocationError

        ! For each entry of the interpolator define the number of derivatives in the AD_D
        do igq = 1,nnodes_f
            do iterm = 1,nterms_s
                interpolator(igq,iterm) = AD_D(nderiv)
            end do
        end do
        
        ! Set the interpolator%x_ad_ equal to the real interpolator
        interpolator = interpolator_rk
        differentiate_me = ( (donor_elem%idomain_g  == fcn_info%seed%idomain_g ) .and. &
                             (donor_elem%ielement_g == fcn_info%seed%ielement_g) .and. &
                             (itime                 == fcn_info%seed%itime) )

        ! Specify the derivatives if node grid sensitivities are require.
        if (fcn_info%dtype == dX_DIFF .and. (interpolation_type /= 'value') .and. differentiate_me) then

            do igq = 1,nnodes_f
                do iterm = 1,nterms_s
                    do idir = 1,3
                        do inode = 1,nnodes_r
                            deriv_index = (idir-1)*nnodes_r + inode
                            interpolator(igq,iterm)%xp_ad_(deriv_index) = interpolator_dx(igq,iterm,inode,idir)
                        end do !inode
                    end do !idir
                end do !iterm
            end do !igq

        end if

        ! Sanity check: this subroutine should not be called if fcn_info%dtype == 'dQ' or 'NO'
        if (fcn_info%dtype == dQ_DIFF .or. fcn_info%dtype == dBC_DIFF .or. fcn_info%dtype == NO_DIFF) then
            call chidg_signal(FATAL,"differentiate_face_chimera_interpolator: this function should not be called &
                                     for dtype = 'dQ' or 'NO'. Implementation error.") 
        end if


        end associate



    end function differentiate_face_chimera_interpolator
    !*****************************************************************************************




end module mod_differentiate
