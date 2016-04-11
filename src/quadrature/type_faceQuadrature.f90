module type_faceQuadrature
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                  XI_DIR,ETA_DIR,ZETA_DIR,SPACEDIM, ZERO
    use mod_polynomial,     only: polynomialVal,dpolynomialVal
    use mod_gaussLegendre,  only: gl_nodes, gl_weights
    use mod_inv,            only: inv

    use type_point,     only: point_t

    implicit none
    private



    !>  
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    type, public :: faceQuadrature_t
        integer(ik)                                 :: nnodes
        type(point_t),  dimension(:,:), allocatable :: nodes    !< (:, iface)
        real(rk),       dimension(:,:), allocatable :: weights  !< (:, iface)

        real(rk), dimension(:,:,:), allocatable :: val      !< (:,:,iface)
        real(rk), dimension(:,:,:), allocatable :: ddxi     !< ''
        real(rk), dimension(:,:,:), allocatable :: ddeta    !< ''
        real(rk), dimension(:,:,:), allocatable :: ddzeta   !< ''

    contains

        procedure :: init

    end type faceQuadrature_t
    !***************************************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine init(self,nnodes_face,nterms)
        class(faceQuadrature_t),    intent(inout)   :: self
        integer(ik),                intent(in)      :: nnodes_face
        integer(ik),                intent(in)      :: nterms

        integer(ik)                             :: inode, iterm, iface, cface
        integer(ik)                             :: nnodes1d,     nterms1d,       nterms2d
        integer(ik)                             :: ixi,          ieta,           izeta
        real(rk)                                :: xi,           eta,            zeta
        real(rk), dimension(:), allocatable     :: xi_vals,      eta_vals,       zeta_vals
        real(rk), dimension(:), allocatable     :: xi_weights,   eta_weights,    zeta_weights
        real(rk), dimension(:), allocatable     :: face_indices
        type(point_t)                           :: node

        self%nnodes = nnodes_face


        !
        ! allocate quadrature storage
        !
        allocate(self%nodes(nnodes_face,NFACES),self%weights(nnodes_face,NFACES))
        allocate(self%val(   nnodes_face,nterms,NFACES))
        allocate(self%ddxi(  nnodes_face,nterms,NFACES))
        allocate(self%ddeta( nnodes_face,nterms,NFACES))
        allocate(self%ddzeta(nnodes_face,nterms,NFACES))



        !
        ! Find number of terms in 1D variable
        !
        nterms1d = 0

        if ( SPACEDIM == 3 ) then
            do while (nterms1d*nterms1d*nterms1d /= nterms)
                nterms1d = nterms1d + 1
            end do
            if (nterms1d*nterms1d*nterms1d > nterms) stop "Error in face quadrature term count"


        else if ( SPACEDIM == 2 ) then
            do while (nterms1d*nterms1d /= nterms)
                nterms1d = nterms1d + 1
            end do
            if (nterms1d*nterms1d > nterms) stop "Error in face quadrature term count"

        else
            call chidg_signal(FATAL,"faceQuadrature%init: Invalid SPACEDIM")

        end if



        !
        ! Initialize quadrature node coordinates for face sets
        !
        ! find number nodes in 1D polynomial
        nnodes1d = 0

        if ( SPACEDIM == 3 ) then

            do while (nnodes1d*nnodes1d /= nnodes_face)
                nnodes1d = nnodes1d + 1
            end do
            if (nnodes1d*nnodes1d > nnodes_face) stop "Error in face quadrature node count"

        else if ( SPACEDIM == 2 ) then
            nnodes1d = nnodes_face
            if (nnodes1d > nnodes_face) stop "Error in face quadrature node count"

        else 
            call chidg_signal(FATAL,"faceQuadrature%init: Invalid SPACEDIM")
        end if


        allocate(xi_vals(nnodes1d),eta_vals(nnodes1d),zeta_vals(nnodes1d))
        allocate(xi_weights(nnodes1d),eta_weights(nnodes1d),zeta_weights(nnodes1d))


        call gl_nodes(nnodes1d,xi_vals)
        call gl_nodes(nnodes1d,eta_vals)
        call gl_nodes(nnodes1d,zeta_vals)
        call gl_weights(nnodes1d,xi_weights)
        call gl_weights(nnodes1d,eta_weights)
        call gl_weights(nnodes1d,zeta_weights)





        if ( SPACEDIM == 3 ) then



            !
            ! xi_min face coordinates
            !
            inode = 1
            do izeta = 1,nnodes1d
                do ieta = 1,nnodes1d
                    xi   = -1._rk
                    eta  = eta_vals(ieta)
                    zeta = zeta_vals(izeta)

                    call self%nodes(inode,XI_MIN)%set(xi,eta,zeta)
                    self%weights(inode,XI_MIN) = eta_weights(ieta)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do



            !
            ! xi_max face coordinates
            !
            inode = 1
            do izeta = 1,nnodes1d
                do ieta = 1,nnodes1d
                    xi   = 1._rk
                    eta  = eta_vals(ieta)
                    zeta = zeta_vals(izeta)

                    call self%nodes(inode,XI_MAX)%set(xi,eta,zeta)
                    self%weights(inode,XI_MAX) = eta_weights(ieta)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do


            
            !
            ! eta_min face coordinates
            !
            inode = 1
            do izeta = 1,nnodes1d
                do ixi = 1,nnodes1d
                    xi   = xi_vals(ixi)
                    eta  = -1._rk
                    zeta = zeta_vals(izeta)

                    call self%nodes(inode,ETA_MIN)%set(xi,eta,zeta)
                    self%weights(inode,ETA_MIN) = xi_weights(ixi)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do

            
            !
            ! eta_max face coordinates
            !
            inode = 1
            do izeta = 1,nnodes1d
                do ixi = 1,nnodes1d
                    xi   = xi_vals(ixi)
                    eta  = 1._rk
                    zeta = zeta_vals(izeta)

                    call self%nodes(inode,ETA_MAX)%set(xi,eta,zeta)
                    self%weights(inode,ETA_MAX) = xi_weights(ixi)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do


            !
            ! zeta_min face coordinates
            !
            inode = 1
            do ieta = 1,nnodes1d
                do ixi = 1,nnodes1d
                    xi   = xi_vals(ixi)
                    eta  = eta_vals(ieta)
                    zeta = -1._rk

                    call self%nodes(inode,ZETA_MIN)%set(xi,eta,zeta)
                    self%weights(inode,ZETA_MIN) = xi_weights(ixi)*eta_weights(ieta)
                    inode = inode + 1
                end do
            end do


            !
            ! zeta_max face coordinates
            !
            inode = 1
            do ieta = 1,nnodes1d
                do ixi = 1,nnodes1d
                    xi   = xi_vals(ixi)
                    eta  = eta_vals(ieta)
                    zeta = 1._rk

                    call self%nodes(inode,ZETA_MAX)%set(xi,eta,zeta)
                    self%weights(inode,ZETA_MAX) = xi_weights(ixi)*eta_weights(ieta)
                    inode = inode + 1
                end do
            end do






        else if ( SPACEDIM == 2 ) then


            !
            ! xi_min face coordinates
            !
            inode = 1
            do ieta = 1,nnodes1d
                xi   = -1._rk
                eta  = eta_vals(ieta)
                zeta = ZERO

                call self%nodes(inode,XI_MIN)%set(xi,eta,zeta)
                self%weights(inode,XI_MIN) = eta_weights(ieta)
                inode = inode + 1
            end do



            !
            ! xi_max face coordinates
            !
            inode = 1
            do ieta = 1,nnodes1d
                xi   = 1._rk
                eta  = eta_vals(ieta)
                zeta = ZERO

                call self%nodes(inode,XI_MAX)%set(xi,eta,zeta)
                self%weights(inode,XI_MAX) = eta_weights(ieta)
                inode = inode + 1
            end do


            
            !
            ! eta_min face coordinates
            !
            inode = 1
            do ixi = 1,nnodes1d
                xi   = xi_vals(ixi)
                eta  = -1._rk
                zeta = ZERO

                call self%nodes(inode,ETA_MIN)%set(xi,eta,zeta)
                self%weights(inode,ETA_MIN) = xi_weights(ixi)
                inode = inode + 1
            end do

            
            !
            ! eta_max face coordinates
            !
            inode = 1
            do ixi = 1,nnodes1d
                xi   = xi_vals(ixi)
                eta  = 1._rk
                zeta = ZERO

                call self%nodes(inode,ETA_MAX)%set(xi,eta,zeta)
                self%weights(inode,ETA_MAX) = xi_weights(ixi)
                inode = inode + 1
            end do


            !
            ! zeta_min face coordinates
            !
            do inode = 1,nnodes1d
                    xi   = ZERO
                    eta  = ZERO
                    zeta = ZERO

                    call self%nodes(inode,ZETA_MIN)%set(xi,eta,zeta)
                    self%weights(inode,ZETA_MIN) = ZERO
            end do


            !
            ! zeta_max face coordinates
            !
            do inode = 1,nnodes1d
                    xi   = ZERO
                    eta  = ZERO
                    zeta = ZERO

                    call self%nodes(inode,ZETA_MAX)%set(xi,eta,zeta)
                    self%weights(inode,ZETA_MAX) = ZERO
            end do





        else
            call chidg_signal(FATAL,"faceQuadrature%init: Invalid SPACEDIM")

        end if





        !===========================================================================
        ! Initialize values and partial derivatives of each modal
        ! polynomial at each face quadrature node
        !===========================================================================

        face_indices = [XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX]


        if ( SPACEDIM == 3 ) then


            do iterm = 1,nterms
                do inode = 1,nnodes_face
                    do iface = 1,size(face_indices)

                        !
                        ! Get current face
                        !                
                        cface = face_indices(iface)

                        !
                        ! on the xi_min face
                        !
                        node = self%nodes(inode,cface)
                        self%val(   inode,iterm,cface) =  polynomialVal(SPACEDIM,nterms,iterm,node)
                        self%ddxi(  inode,iterm,cface) = dpolynomialVal(SPACEDIM,nterms,iterm,node,XI_DIR)
                        self%ddeta( inode,iterm,cface) = dpolynomialVal(SPACEDIM,nterms,iterm,node,ETA_DIR)
                        self%ddzeta(inode,iterm,cface) = dpolynomialVal(SPACEDIM,nterms,iterm,node,ZETA_DIR)

                    end do
                end do
            end do


        else if ( SPACEDIM == 2 ) then


            do iterm = 1,nterms
                do inode = 1,nnodes_face
                    do iface = 1,size(face_indices)

                        !
                        ! Get current face
                        !                
                        cface = face_indices(iface)

                        !
                        ! on the xi_min face
                        !
                        node = self%nodes(inode,cface)
                        self%val(   inode,iterm,cface) =  polynomialVal(SPACEDIM,nterms,iterm,node)
                        self%ddxi(  inode,iterm,cface) = dpolynomialVal(SPACEDIM,nterms,iterm,node,XI_DIR)
                        self%ddeta( inode,iterm,cface) = dpolynomialVal(SPACEDIM,nterms,iterm,node,ETA_DIR)
                        self%ddzeta(inode,iterm,cface) = ZERO

                    end do
                end do
            end do


        else
            call chidg_signal(FATAL,"faceQuadrature%init: Invalid SPACEDIM")


        end if


    end subroutine init
    !*******************************************************************************************************
    




end module type_faceQuadrature
