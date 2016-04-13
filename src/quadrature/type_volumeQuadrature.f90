module type_volumeQuadrature
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                      XI_DIR,ETA_DIR,ZETA_DIR,ZERO, TWO_DIM, THREE_DIM, ONE

    use mod_GaussLegendre,      only: gl_nodes, gl_weights
    use mod_polynomial,         only: polynomialVal, dpolynomialVal
    use mod_quadrature_tools,   only: compute_nnodes1d_volume
    use type_point,             only: point_t

    implicit none
    private

    !> Type defining volume quadrature
    !!      - Contains nodes, weights, and matrices for Gauss-quadrature. Used for
    !!        integration and projection routines. An instance is defined for a specified
    !!        number of quadrature nodes and specified number of terms in a polynomial expansion
    !!
    !!  @author Nathan A. Wukie
    !-----------------------------------------------------------------------------------------------------------
    type, public :: volumeQuadrature_t
        ! Number of volume quadrature nodes
        integer(ik)                 :: nnodes       !< Number of volume quadrature nodes
        integer(ik)                 :: nterms       !< Number of terms in the polynomial expansion

        ! Array of points and weights for each node
        type(point_t),  allocatable :: nodes(:)     !< Array of quadrature node points
        real(rk),       allocatable :: weights(:)   !< Array of quadrature node weights

        ! Polynomial values and derivatives for each node
        real(rk),       allocatable :: val(:,:)     !< Matrix used to interpolate an expansion to quadrature nodes
        real(rk),       allocatable :: ddxi(:,:)    !< Matrix used to interpolate partial derivatives(ddxi) to quadrature nodes
        real(rk),       allocatable :: ddeta(:,:)   !< Matrix used to interpolate partial derivatives(ddeta) to quadrature nodes
        real(rk),       allocatable :: ddzeta(:,:)  !< Matrix used to interpolate partial derivatives(ddzeta) to quadrature nodes

        ! Reference mass matrix
        real(rk),       allocatable :: mass(:,:)    !< Mass matrix for the reference element
        real(rk),       allocatable :: dmass(:)     !< Diagonal of the mass matrix for the reference element

    contains

        procedure :: init

    end type volumeQuadrature_t
    !***********************************************************************************************************





contains



    !> Initialization routine for volumeQuadrature_t instance.
    !!      - Allocates storage for member data. Initializes ndoes, weights, and interpolation matrices.
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  nnodes  Number of nodes used for Gauss-quadrature
    !!  @param[in]  nterms  Number of terms in the associated polynomial expansion
    !!
    !!
    !!  TODO: TEST SPACEDIM
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine init(self,spacedim,nnodes,nterms)
        class(volumeQuadrature_t),  intent(inout)   :: self
        integer(ik),                intent(in)      :: spacedim
        integer(ik),                intent(in)      :: nnodes
        integer(ik),                intent(in)      :: nterms


        integer(ik)                                 :: ixi,ieta,izeta,inode,iterm,nnodes1d,ierr
        integer(ik)                                 :: nnodes_xi, nnodes_eta, nnodes_zeta
        real(rk)                                    :: xi,eta,zeta
        real(rk), dimension(:), allocatable         :: xi_vals,eta_vals,zeta_vals
        real(rk), dimension(:), allocatable         :: xi_weights,eta_weights,zeta_weights
        real(rk), dimension(nterms,nnodes)          :: temp
        type(point_t)                               :: node

        self%nnodes = nnodes
        self%nterms = nterms


        !
        ! allocate quadrature storage
        !
        allocate(self%nodes(nnodes),self%weights(nnodes),       &
                 self%val(   nnodes,nterms),                    &
                 self%ddxi(  nnodes,nterms),                    &
                 self%ddeta( nnodes,nterms),                    &
                 self%ddzeta(nnodes,nterms),                    &
                 self%mass(nterms,nterms),                      &
                 self%dmass(nterms), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize quadrature node coordinates for face sets
        !
        ! find number nodes in 1D polynomial
        nnodes1d = compute_nnodes1d_volume(spacedim,nnodes)

        allocate(xi_vals(nnodes1d),eta_vals(nnodes1d),zeta_vals(nnodes1d))
        allocate(xi_weights(nnodes1d),eta_weights(nnodes1d),zeta_weights(nnodes1d))

        call gl_nodes(nnodes1d,xi_vals)
        call gl_nodes(nnodes1d,eta_vals)
        call gl_nodes(nnodes1d,zeta_vals)
        call gl_weights(nnodes1d,xi_weights)
        call gl_weights(nnodes1d,eta_weights)
        call gl_weights(nnodes1d,zeta_weights)



        !
        ! Specialize for 2D, 3D
        !
        if ( spacedim == THREE_DIM ) then
            nnodes_zeta = nnodes1d
        else if ( spacedim == TWO_DIM ) then
            nnodes_zeta  = 1
            zeta_weights = ONE
            zeta_vals    = ZERO
        end if






        !
        ! Volume node coordinates and weights
        !
        inode = 1
        do izeta = 1,nnodes_zeta        ! specialized for 2D/3D
            do ieta = 1,nnodes1d
                do ixi = 1,nnodes1d
                    xi   = xi_vals(  ixi)
                    eta  = eta_vals( ieta)
                    zeta = zeta_vals(izeta)

                    call self%nodes(inode)%set(xi,eta,zeta)
                    self%weights(inode) = xi_weights(ixi)*eta_weights(ieta)*zeta_weights(izeta)

                    inode = inode + 1
                end do
            end do
        end do






        !
        ! Initialize values and partial derivatives of each modal
        ! polynomial at each volume quadrature node
        !
        do iterm = 1,nterms
            do inode = 1,nnodes
                    node = self%nodes(inode)
                    self%val(   inode,iterm) =  polynomialVal(spacedim,nterms,iterm,node)
                    self%ddxi(  inode,iterm) = dpolynomialVal(spacedim,nterms,iterm,node,XI_DIR)
                    self%ddeta( inode,iterm) = dpolynomialVal(spacedim,nterms,iterm,node,ETA_DIR)
                    self%ddzeta(inode,iterm) = dpolynomialVal(spacedim,nterms,iterm,node,ZETA_DIR)
            end do
        end do





        !
        !   Initialize reference mass matrix
        !
        temp = transpose(self%val)
        do iterm = 1,nterms
            temp(iterm,:) = temp(iterm,:)*(self%weights)
        end do

        self%mass = matmul(temp,self%val)



        !
        ! Set reference mass diagonal
        !
        do iterm = 1,nterms
            self%dmass(iterm) = self%mass(iterm,iterm)
        end do



    end subroutine init
    !***********************************************************************************************************




    

end module type_volumeQuadrature
