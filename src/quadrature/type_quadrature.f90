module type_quadrature
    use mod_kinds,              only: rk,ik
    use type_volumeQuadrature,  only: volumeQuadrature_t
    use type_faceQuadrature,    only: faceQuadrature_t
    implicit none

    private



    !> Primary quadrature container. Has volume and face quadrature components. Facilitates initialization
    !! of volume/face quadrature components.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    type, public :: quadrature_t

        integer(ik)                 :: spacedim     !< Number of spatial dimensions the quadrature set is initialialized for.
        integer(ik)                 :: nnodes_v     !< Number of nodes in the volume quadrature set.
        integer(ik)                 :: nnodes_f     !< Number of nodes in the face quadrature set.
        integer(ik)                 :: nterms       !< Number of terms in the expansion being integrated.

        type(volumeQuadrature_t)    :: vol          !< Volume quadrature instance
        type(faceQuadrature_t)      :: face         !< Face quadrature instance

        logical                     :: isInitialized = .false.

    contains

        procedure   :: init

    end type quadrature_t
    !**********************************************************************************************************




contains






    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine init(self,spacedim,nnodes_face,nnodes_vol,nterms)
        class(quadrature_t), intent(inout) :: self
        integer(ik),         intent(in)    :: spacedim
        integer(ik),         intent(in)    :: nnodes_face
        integer(ik),         intent(in)    :: nnodes_vol
        integer(ik),         intent(in)    :: nterms

        !
        ! Set nnodes, nterms for the quadrature instance.
        !
        self%spacedim = spacedim
        self%nnodes_v = nnodes_vol
        self%nnodes_f = nnodes_face
        self%nterms   = nterms

        !
        ! Initialize volume/face quadrature instances
        !
        call self%vol%init(spacedim,nnodes_vol,nterms)
        call self%face%init(spacedim,nnodes_face,nterms)


        self%isInitialized = .true.

    end subroutine init
    !**********************************************************************************************************
    







end module type_quadrature
