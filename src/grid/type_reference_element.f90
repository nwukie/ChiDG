module type_reference_element
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, ZERO, ONE
    use mod_gauss_legendre, only: quadrature_nodes, quadrature_weights
    use mod_nodes_uniform,  only: uniform_nodes
    implicit none


!    call reference_element%init_element(element = 'Hexahedral8')
!    call reference_element%init_nodes(nodes='uniform', level=1, dim=3)
!    call reference_element%init_interpolator(nterms=8)
!
!    nodes   = reference_element%nodes()
!    weights = reference_element%weights()
!    val     = reference_element%interpolator('Value')
!    ddxi    = reference_element%interpolator('Grad1')
!    ddeta   = reference_element%interpolator('Grad2')
!    ddzeta  = reference_element%interpolator('Grad3')




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/26/2017
    !!
    !---------------------------------------------------------------------
    type, public :: reference_element_t

        character(:),   allocatable :: type
        character(:),   allocatable :: polynomial

        real(rk),       allocatable :: e_nodes_(:,:)        ! nodes_3d(nnodes, 3). xi = nodes(:,1), eta = nodes(:,2), zeta = nodes(:,3)
        real(rk),       allocatable :: e_weights_(:)        ! weights_3d(nnodes)
        real(rk),       allocatable :: f_nodes_(:,:,:)      ! nodes_2d(nnodes, 3, nfaces). xi = nodes(:,1,iface), eta = nodes(:,2,iface), zeta = nodes(:,3,iface)
        real(rk),       allocatable :: f_weights_(:)        ! weights_2d(nnodes)


        real(rk),       allocatable :: e_val(:,:)           ! element interpolator, expansion value  to elem nodes
        real(rk),       allocatable :: e_ddxi(:,:)          ! element interpolator, expansion ddxi   to elem nodes
        real(rk),       allocatable :: e_ddeta(:,:)         ! element interpolator, expansion ddeta  to elem nodes
        real(rk),       allocatable :: e_ddzeta(:,:)        ! element interpolator, expansion ddzeta to elem nodes
        real(rk),       allocatable :: f_val(:,:,:)         ! face    interpolator, expansion value  to face nodes
        real(rk),       allocatable :: f_ddxi(:,:,:)        ! face    interpolator, expansion ddxi   to face nodes
        real(rk),       allocatable :: f_ddeta(:,:,:)       ! face    interpolator, expansion ddeta  to face nodes
        real(rk),       allocatable :: f_ddzeta(:,:,:)      ! face    interpolator, expansion ddzeta to face nodes
        real(rk),       allocatable :: nodes_to_modes(:,:)  ! linear projector

        logical :: nodes_initialized   = .false.
        logical :: weights_initialized = .false.
    contains

        ! Initialize
        procedure   :: init_element
        procedure   :: init_polynomial
        procedure   :: init_nodes

        ! Produce 
        generic     :: nodes        => nodes_element,   nodes_face
        procedure   :: nodes_element
        procedure   :: nodes_face
        generic     :: weights      => weights_element, weights_face
        procedure   :: weights_element
        procedure   :: weights_face
        !procedure   :: interpolator => interpolate_element, interpolate_face

    end type reference_element_t
    !********************************************************************






contains




    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    subroutine init_element(self,element_type)
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: element_type



    end subroutine init_element
    !*******************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    subroutine init_nodes(self, nodes, level)
        class(reference_element_t), intent(inout)   :: self
        character(*),               intent(in)      :: nodes
        integer(ik),                intent(in)      :: level

        integer(ik)             :: iface, ierr
        real(rk),   allocatable :: nodes_3d(:,:), nodes_2d(:,:)
        real(rk),   allocatable :: weights_3d(:), weights_2d(:)

        !
        ! Reset node/weight initialization flags
        !
        self%nodes_initialized   = .false.
        self%weights_initialized = .false.


        !
        ! Generate 3D/2D node sets
        !
        select case(trim(nodes))
            case('Uniform')
                nodes_3d = uniform_nodes(level,dim=3)
                nodes_2d = uniform_nodes(level,dim=2)


            case('Quadrature')
                nodes_3d   = quadrature_nodes(  level=level,dim=3)
                nodes_2d   = quadrature_nodes(  level=level,dim=2)
                weights_3d = quadrature_weights(level=level,dim=3)
                weights_2d = quadrature_weights(level=level,dim=2)

                self%weights_initialized = .true.
            case default
                call chidg_signal_one(FATAL,"reference_element%init_nodes: requested node set not implemented.", trim(nodes))
        end select


        !
        ! Set nodes (if provided)
        !
        if (allocated(nodes_3d) .and. allocated(nodes_2d)) then
            ! Set 3D element nodes
            self%e_nodes_ = nodes_3d

            ! Set 2D face nodes
            if (allocated(self%f_nodes_)) deallocate(self%f_nodes_)
            allocate(self%f_nodes_(size(nodes_2d,1),3,NFACES), stat=ierr)
            if (ierr /= 0) call AllocationError

            iface = 1
            self%f_nodes_(:,1,iface) = -ONE
            self%f_nodes_(:,2,iface) = nodes_2d(:,1)
            self%f_nodes_(:,3,iface) = nodes_2d(:,2)

            iface = 2
            self%f_nodes_(:,1,iface) = ONE
            self%f_nodes_(:,2,iface) = nodes_2d(:,1)
            self%f_nodes_(:,3,iface) = nodes_2d(:,2)

            iface = 3
            self%f_nodes_(:,1,iface) = nodes_2d(:,1)
            self%f_nodes_(:,2,iface) = -ONE
            self%f_nodes_(:,3,iface) = nodes_2d(:,2)

            iface = 4
            self%f_nodes_(:,1,iface) = nodes_2d(:,1)
            self%f_nodes_(:,2,iface) = ONE
            self%f_nodes_(:,3,iface) = nodes_2d(:,2)

            iface = 5
            self%f_nodes_(:,1,iface) = nodes_2d(:,1)
            self%f_nodes_(:,2,iface) = nodes_2d(:,2)
            self%f_nodes_(:,3,iface) = -ONE

            iface = 6
            self%f_nodes_(:,1,iface) = nodes_2d(:,1)
            self%f_nodes_(:,2,iface) = nodes_2d(:,2)
            self%f_nodes_(:,3,iface) = ONE


            self%nodes_initialized = .true.

        end if !nodes



        !
        ! Set weights (if provided)
        !
        if (allocated(self%e_weights_)) deallocate(self%e_weights_)
        if (allocated(self%f_weights_)) deallocate(self%f_weights_)
        if (allocated(weights_3d) .and. allocated(weights_2d)) then
            ! Set 3D element weights
            self%e_weights_ = weights_3d

            ! Set 2D element weights
            self%f_weights_ = weights_2d

            self%weights_initialized = .true.
        end if





    end subroutine init_nodes
    !*******************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    subroutine init_polynomial(self, polynomial, nterms)
        class(reference_element_t), intent(inout)   :: self
        character(*),               intent(in)      :: polynomial
        integer(ik),                intent(in)      :: nterms




    end subroutine init_polynomial
    !*******************************************************************





    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function nodes_element(self) result(nodes_)
        class(reference_element_t), intent(in)  :: self

        real(rk),       allocatable :: nodes_(:,:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return element nodes without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%nodes_initialized) call chidg_signal(FATAL,user_msg)

        nodes_ = self%e_nodes_

    end function nodes_element
    !*******************************************************************


    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function nodes_face(self,iface) result(nodes_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iface

        real(rk),       allocatable :: nodes_(:,:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return face nodes without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%nodes_initialized) call chidg_signal(FATAL,user_msg)

        nodes_ = self%f_nodes_(:,:,iface)

    end function nodes_face
    !*******************************************************************




    !>
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------
    function weights_element(self) result(weights_)
        class(reference_element_t), intent(in)  :: self

        real(rk),       allocatable :: weights_(:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return element weights without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%weights_initialized) call chidg_signal(FATAL,user_msg)

        weights_ = self%e_weights_

    end function weights_element
    !*******************************************************************




    !>
    !!
    !!  NOTE: we don't actually use the 'iface' argument in the call.
    !!  it is just there to be consistent with the other calls and it
    !!  also helps differentiate this routine from the weights_element
    !!  routine for the generic interface.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!
    !!
    !-------------------------------------------------------------------
    function weights_face(self,iface) result(weights_)
        class(reference_element_t), intent(in)  :: self
        integer(ik),                intent(in)  :: iface

        real(rk),       allocatable :: weights_(:)
        character(:),   allocatable :: user_msg

        user_msg = "reference_element_t: tried to return face weights without being &
                    initialized. Make sure init_nodes is getting called."
        if (.not. self%weights_initialized) call chidg_signal(FATAL,user_msg)

        weights_ = self%f_weights_

    end function weights_face
    !*******************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolate_element(self,selector) result(interpolator_)  
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: selector

        real(rk),   allocatable :: interpolator_(:,:)


        select case(trim(selector))
            case('Value')
                interpolator_ = self%e_val
            case('Grad1')
                interpolator_ = self%e_ddxi
            case('Grad2')
                interpolator_ = self%e_ddeta
            case('Grad3')
                interpolator_ = self%e_ddzeta
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolate: Invalid selector for element interpolator.", trim(selector))
        end select
        
    end function interpolate_element
    !*******************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------
    function interpolate_face(self,selector,iface) result(interpolator_)  
        class(reference_element_t), intent(in)  :: self
        character(*),               intent(in)  :: selector
        integer(ik),                intent(in)  :: iface

        real(rk),   allocatable :: interpolator_(:,:)


        select case(trim(selector))
            case('Value')
                interpolator_ = self%f_val(:,:,iface)
            case('Grad1')
                interpolator_ = self%f_ddxi(:,:,iface)
            case('Grad2')
                interpolator_ = self%f_ddeta(:,:,iface)
            case('Grad3')
                interpolator_ = self%f_ddzeta(:,:,iface)
            case default
                call chidg_signal_one(FATAL,"reference_element%interpolate: Invalid selector for face interpolator.", trim(selector))
        end select
        
    end function interpolate_face
    !*******************************************************************












end module type_reference_element
