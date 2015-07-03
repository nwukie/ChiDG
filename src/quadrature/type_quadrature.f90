module type_quadrature
    use mod_kinds,              only: rk,ik

    use type_volumeQuadrature,  only: volumeQuadrature_t
    use type_faceQuadrature,    only: faceQuadrature_t

    implicit none
    private

    type, public :: quadrature_t
        integer(ik) :: nnodes_vol
        integer(ik) :: nterms

        type(volumeQuadrature_t)    :: vol  !> Volume quadrature instance
        type(faceQuadrature_t)      :: face !> Face quadrature instance

        logical                     :: isInitialized = .false.

    contains
        procedure   :: init
        final       :: destructor
    end type quadrature_t

contains

    subroutine init(self,nnodes_face,nnodes_vol,nterms)
        class(quadrature_t), intent(inout) :: self
        integer(ik),         intent(in)    :: nnodes_face
        integer(ik),         intent(in)    :: nnodes_vol
        integer(ik),         intent(in)    :: nterms

        self%nnodes_vol = nnodes_vol
        self%nterms     = nterms

        call self%vol%init(nnodes_vol,nterms)       !> Initialize volume quadrature
        call self%face%init(nnodes_face,nterms)     !> Initialize face quadrature


        self%isInitialized = .true.
    end subroutine
    
    subroutine destructor(self)
        type(quadrature_t), intent(inout) :: self
    end subroutine

end module type_quadrature
