module type_BR2_lift
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_constants,      only: BR2_INTERIOR_LOCATION, BR2_INTERIOR, BR2_EXTERIOR, CHIMERA, INTERIOR, &
                                  BOUNDARY_DIFFUSIVE_FLUX, DIAG, HALF, ME, NEIGHBOR
    use mod_interpolate,    only: interpolate_face_autodiff
    use mod_DNAD_tools,     only: face_compute_seed, compute_neighbor_face
    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use type_chidgVector,   only: chidgVector_t
    use DNAD_D
    implicit none







    !>  Container holding a lifting operator for a single equation.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !------------------------------------------------------------------------------------
    type, public :: BR2_lift_t

        type(AD_D), allocatable :: lift(:,:)    ! nterms, 3 coordinate directions

    contains

        procedure   :: init
        procedure   :: update

    end type BR2_lift_t
    !************************************************************************************






contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self,mesh,elem_info,iface,idiff)
        class(BR2_lift_t),      intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: idiff

        integer(ik) :: nterms_s, idom_n, ielem_n, idonor, ierr, ChiID
        logical     :: differentiating_interior, conforming_face, chimera_face, boundary_face
        logical     :: allocate_lift, reallocate_lift


        differentiating_interior = (idiff == BR2_INTERIOR_LOCATION)

        !
        ! Get number of equations and number of terms in solution expansion
        !
        if (differentiating_interior) then
            nterms_s = mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%nterms_s

        else


            ! Get number of dependent elements
            conforming_face = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == INTERIOR)
            chimera_face    = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == CHIMERA )

            if (conforming_face) then
                idom_n  = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ineighbor_domain_l
                ielem_n = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ineighbor_element_l

                nterms_s = mesh(idom_n)%elems(ielem_n)%nterms_s

            else if (chimera_face) then
                idonor = idiff - 1  ! idiff = interior, idonor1, idonor2, etc...
                ChiID = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ChiID

                nterms_s = mesh(elem_info%idomain_l)%chimera%recv%data(ChiID)%donor_nterms_s%at(idonor)

            else
                call chidg_signal(FATAL,"BR2_lift%init: Invalid face type for BR2_lift initialization")
            end if


        end if




        !
        ! (Re)Allocate storage for differentiated polynomial expansion of lifting operator vector
        !
        ! Lifting operator is a vector, so (nterms,ncoords)
        !
        allocate_lift = (.not. allocated(self%lift))

        if (allocate_lift) then
            allocate(self%lift(nterms_s,3), stat=ierr)
            if (ierr /= 0) call AllocationError

        else

            reallocate_lift = (nterms_s /= size(self%lift,1))

            if (reallocate_lift) then
                deallocate(self%lift)
                allocate(self%lift(nterms_s,3), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        end if


    end subroutine init
    !*****************************************************************************************

















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine update(self,mesh,elem_info,iface,idiff_BR2,ieqn,q,BR2_TYPE)
        class(BR2_lift_t),      intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: idiff_BR2
        integer(ik),            intent(in)      :: ieqn
        type(chidgVector_t),    intent(in)      :: q
        integer(ik),            intent(in)      :: BR2_TYPE

        type(face_info_t)       :: face_info
        type(function_info_t)   :: function_info
        integer(ik)             :: idiff_dir, idepend, idom_n, ielem_n, iface_n
        logical                 :: conforming_face, chimera_face

        type(AD_D), allocatable, dimension(:)   ::  &
                var_m,      var_p,      diff,       &
                diff_x,     diff_y,     diff_z,     &
                integral_x, integral_y, integral_z

        !
        ! Set up face info for interpolation
        !
        face_info%idomain_g  = elem_info%idomain_g
        face_info%idomain_l  = elem_info%idomain_l
        face_info%ielement_g = elem_info%ielement_g
        face_info%ielement_l = elem_info%ielement_l
        face_info%iface      = iface


        !
        ! Set up function info for differentiation
        !
        if (idiff_BR2 == 1) then
            idiff_dir = DIAG
            idepend   = 1
        else
            idiff_dir = iface
            idepend = idiff_BR2 - 1
        end if



        function_info%type    = BOUNDARY_DIFFUSIVE_FLUX ! I don't think this matters here, just want to give a valid entry.
        function_info%ifcn    = 1                       ! I don't think this matters here, just want to give a valid entry.
        function_info%idiff   = idiff_dir
        function_info%idepend = idepend
        function_info%seed    = face_compute_seed(mesh,face_info%idomain_l,face_info%ielement_l,iface,idepend,idiff_dir)



        !
        ! Interpolate variable to face quadrature nodes
        !
        conforming_face = ( mesh(face_info%idomain_l)%faces(face_info%ielement_l,iface)%ftype == INTERIOR )
        chimera_face    = ( mesh(face_info%idomain_l)%faces(face_info%ielement_l,iface)%ftype == CHIMERA  )

        if (conforming_face .or. chimera_face) then


            if (BR2_TYPE == BR2_INTERIOR) then
                var_m = interpolate_face_autodiff(mesh,q,face_info,function_info, ieqn, 'value', ME)
                var_p = interpolate_face_autodiff(mesh,q,face_info,function_info, ieqn, 'value', NEIGHBOR)
            else if (BR2_TYPE == BR2_EXTERIOR) then
                var_p = interpolate_face_autodiff(mesh,q,face_info,function_info, ieqn, 'value', ME)
                var_m = interpolate_face_autodiff(mesh,q,face_info,function_info, ieqn, 'value', NEIGHBOR)
            else
                call chidg_signal(FATAL,"BR2_lift%update: Invalid BR2_TYPE selection. Options are BR2_INTERIOR and BR2_EXTERIOR.")
            end if


            diff = HALF*(var_p - var_m)



            !
            ! Dot with normal vector
            !
            if (BR2_TYPE == BR2_INTERIOR) then
                associate ( norm     => mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%norm,                         &
                            weights  => mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%gq%face%weights(:,iface),     &
                            valtrans => mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%gq%face%val_trans(:,:,iface), &
                            invmass  => mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%invmass )

                diff_x = diff * norm(:,1) * weights
                diff_y = diff * norm(:,2) * weights
                diff_z = diff * norm(:,3) * weights

                integral_x = matmul(valtrans,diff_x)
                integral_y = matmul(valtrans,diff_y)
                integral_z = matmul(valtrans,diff_z)

                self%lift(:,1) = matmul(invmass,integral_x)
                self%lift(:,2) = matmul(invmass,integral_y)
                self%lift(:,3) = matmul(invmass,integral_z)

                end associate




            else if (BR2_TYPE == BR2_EXTERIOR) then


                if (conforming_face) then
                    idom_n  = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ineighbor_domain_l
                    ielem_n = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ineighbor_element_l
                    iface_n = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%get_neighbor_face()

                    associate ( norm     => mesh(idom_n)%faces(ielem_n,iface_n)%norm,                           &
                                weights  => mesh(idom_n)%faces(ielem_n,iface_n)%gq%face%weights(:,iface_n),     &
                                valtrans => mesh(idom_n)%faces(ielem_n,iface_n)%gq%face%val_trans(:,:,iface_n), &
                                invmass  => mesh(idom_n)%elems(ielem_n)%invmass )

                    diff_x = diff * norm(:,1) * weights
                    diff_y = diff * norm(:,2) * weights
                    diff_z = diff * norm(:,3) * weights

                    integral_x = matmul(valtrans,diff_x)
                    integral_y = matmul(valtrans,diff_y)
                    integral_z = matmul(valtrans,diff_z)

                    self%lift(:,1) = matmul(invmass,integral_x)
                    self%lift(:,2) = matmul(invmass,integral_y)
                    self%lift(:,3) = matmul(invmass,integral_z)

                    end associate

                else if (chimera_face) then
                    call chidg_signal(FATAL,"BR2_lift%update: chimera not yet implemented")
                end if

            end if

        end if


    end subroutine update
    !***************************************************************************************************
























end module type_BR2_lift
