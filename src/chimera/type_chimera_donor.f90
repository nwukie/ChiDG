module type_chimera_donor
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NO_ID, CARTESIAN, CYLINDRICAL, XI_DIR, ETA_DIR, ZETA_DIR, ZERO
    use type_element_info,  only: element_info_t
    use mod_polynomial,     only: polynomial_val, dpolynomial_val
    use mod_inv,            only: inv, inv_3x3, dinv_3x3
    implicit none



    !>  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !!
    !---------------------------------------------------------------
    type, public :: chimera_donor_t

        type(element_info_t)    :: elem_info

        ! Donor info
        real(rk),    allocatable    :: nodes_to_modes(:,:)  ! Donor nodes_to_modes matrix
        real(rk),    allocatable    :: modal_coords(:,:)    ! Donor modal coordinate
!        real(rk),    allocatable    :: modal_coords_def(:,:)! Donor modal deformed coordinate, NOT PASSED YET


        ! Node information
        integer(ik),    allocatable :: node_index(:)    ! index in a node set the donor is providing data for
        real(rk),       allocatable :: coords(:,:)      ! donor-local node coordinate(xi,eta,zeta)
        real(rk),       allocatable :: metric(:,:,:)    ! For each node, a matrix of metric terms
        real(rk),       allocatable :: jinv(:)          ! For each node, inverse element jacobian


        ! Interpolators
        real(rk),       allocatable :: value(:,:)
        real(rk),       allocatable :: grad1(:,:)
        real(rk),       allocatable :: grad2(:,:)
        real(rk),       allocatable :: grad3(:,:)

        real(rk),       allocatable :: dgrad1_dx(:,:,:,:)
        real(rk),       allocatable :: dgrad2_dx(:,:,:,:)
        real(rk),       allocatable :: dgrad3_dx(:,:,:,:)


    contains

        procedure   :: add_node
        procedure   :: set_properties
        procedure   :: nnodes

        procedure   :: update_interpolations_dx
        procedure   :: release_interpolations_dx

    
    end type chimera_donor_t
    !***************************************************************


    interface chimera_donor
        module procedure chimera_donor
    end interface


contains


    !>  Contructor for chimera_donor_t
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !------------------------------------------------------------------
    function chimera_donor(donor,nodes_to_modes,modal_coords) result(instance)
        type(element_info_t),                   intent(in)  :: donor
        real(rk),               allocatable,    intent(in)  :: nodes_to_modes(:,:) 
        real(rk),               allocatable,    intent(in)  :: modal_coords(:,:) 

        type(chimera_donor_t)   :: instance

        instance%elem_info      = donor
        instance%nodes_to_modes = nodes_to_modes
        instance%modal_coords   = modal_coords

    end function chimera_donor
    !******************************************************************





    !>  Add a node to the donor.
    !!
    !!  For some node set on the receiver face:
    !!
    !!      |--o---o-----o-----o---o--|
    !!         1   2     3     4   5
    !!
    !!  We are adding one of the nodes that is being provided 
    !!  by the current donor object.
    !!
    !!
    !!  For example, inode is the index of the node as it exists
    !!  in the face node set. This way, we know where to distribute
    !!  data to.
    !!
    !!  So, if the node was the 4th in the face node set, inode would
    !!  equal 4.
    !!
    !!                       inode
    !!      |--o---o-----o-----o---o--|
    !!         1   2     3     4   5
    !!
    !!
    !!  coord(3) is the [xi,eta,zeta] location of where the node exists 
    !!  within the donor element. In this way, we can evaluate quantities
    !!  on the donor element at [xi,eta,zeta] that correspond to the 
    !!  physical location on the receiver face.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !----------------------------------------------------------------
    subroutine add_node(self,inode,coord,metric,jinv)
        class(chimera_donor_t), intent(inout)   :: self
        integer(ik),            intent(in)      :: inode
        real(rk),               intent(in)      :: coord(3)
        real(rk),               intent(in)      :: metric(:,:)
        real(rk),               intent(in)      :: jinv

        integer(ik)                 :: ierr
        integer(ik),    allocatable :: tmp_nodes(:)
        real(rk),       allocatable :: tmp_coords(:,:)
        real(rk),       allocatable :: tmp_metric(:,:,:)
        real(rk),       allocatable :: tmp_jinv(:)


        
        !
        ! Extend allocations
        !
        allocate(tmp_nodes(      self%nnodes() + 1),         &
                 tmp_coords(     self%nnodes() + 1, 3),      &
                 tmp_jinv(       self%nnodes() + 1),         &
                 tmp_metric(3,3, self%nnodes() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy existing data
        !
        if (self%nnodes() > 0 ) then
            tmp_nodes(1:self%nnodes())      = self%node_index(1:self%nnodes())
            tmp_coords(1:self%nnodes(),:)   = self%coords(1:self%nnodes(),:)
            tmp_jinv(1:self%nnodes())       = self%jinv(1:self%nnodes())
            tmp_metric(:,:,1:self%nnodes()) = self%metric(:,:,1:self%nnodes())
        end if


        !
        ! Move allocation
        !
        call move_alloc(tmp_nodes,  self%node_index)
        call move_alloc(tmp_coords, self%coords)
        call move_alloc(tmp_metric, self%metric)
        call move_alloc(tmp_jinv,   self%jinv)


        !
        ! Add new data to end
        !
        self%node_index(self%nnodes()) = inode
        self%coords(self%nnodes(),:)   = coord
        self%metric(:,:,self%nnodes()) = metric
        self%jinv(self%nnodes())       = jinv


    end subroutine add_node
    !****************************************************************





    !>  Set the properties of the donor.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !!  Added nodes_to_modes and modal_coords as input
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/30/2018
    !!
    !------------------------------------------------------------------
    subroutine set_properties(self,nodes_to_modes,modal_coords)
        class(chimera_donor_t),     intent(inout)   :: self
        real(rk),   allocatable,    intent(in)      :: nodes_to_modes(:,:) 
        real(rk),   allocatable,    intent(in)      :: modal_coords(:,:) 

        self%nodes_to_modes     = nodes_to_modes
        self%modal_coords       = modal_coords

    end subroutine set_properties
    !******************************************************************





    !>  Return the number of nodes this donor is responsible for.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !----------------------------------------------------------------
    function nnodes(self) result(n)
        class(chimera_donor_t),    intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%node_index)) then
            n = size(self%node_index)
        else
            n = 0
        end if

    end function nnodes
    !****************************************************************



    !>  Compute derivatives of gradient interpolators
    !!
    !!  FIXME: Do not use with deformed coordinate frame. Fix this with Nathan
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/14/2018
    !!
    !!  TODO: tests. Already indirectly tested in test_update_space_linearadvection.pf
    !!        Might need specific test
    !!
    !----------------------------------------------------------------
    subroutine update_interpolations_dx(self,coordinate_frame,coordinate_scaling)
        class(chimera_donor_t), intent(inout)           :: self
        character(*),           intent(in), optional    :: coordinate_frame
        logical,                intent(in), optional    :: coordinate_scaling
        
        real(rk)                             :: r, dXdxi(3,3), coord(3), ddxi_s,    &
                                                ddeta_s, ddzeta_s
        real(rk), dimension(self%elem_info%nterms_c)   :: ddxi_c, ddeta_c, ddzeta_c, val_c
        integer(ik)                          :: iterm, nnodes_r, idiff_n, icoord,    &
                                                ierr, ipt, spacedim, idir
        logical                              :: scale_metric
        real(rk),       allocatable          :: modal_coords(:,:), dmodes(:,:)
        character(:),   allocatable          :: frame_selector
        real(rk),       allocatable          :: djacobian_dx(:,:,:,:),dmetric_dx(:,:,:,:)


        spacedim = 3

        ! Check coordiante system is correct
        if (self%elem_info%coordinate_system == 0) then
            call chidg_signal(FATAL,"type_chimera_donor: coordiante_system not assigned to the donor")
        end if

        ! Set default values for optional inputs
        frame_selector = 'Undeformed'
        scale_metric   = .true.


        ! Overwrite optional values if inputs are present
        if (present(coordinate_scaling)) scale_metric   = coordinate_scaling
        if (present(coordinate_frame)  ) frame_selector = coordinate_frame


        ! Retrieve modal coordinates
        select case(trim(frame_selector))
            case('Undeformed')
                modal_coords = self%modal_coords
            case('Deformed')
                !modal_coords = self%modal_coords_def
                call write_line("update_interpolations_dx: deformed not yet implemented")
        end select


        ! Allocate differential gradient interpolator 
        if (allocated(self%dgrad1_dx)) deallocate(self%dgrad1_dx,self%dgrad2_dx,self%dgrad3_dx)
        !allocate( self%dgrad1_dx(self%nnodes(),self%elem_info%nterms_s,self%nnodes_r,3), &
        !          self%dgrad2_dx(self%nnodes(),self%elem_info%nterms_s,self%nnodes_r,3), &
        !          self%dgrad3_dx(self%nnodes(),self%elem_info%nterms_s,self%nnodes_r,3)  ) 
        allocate( self%dgrad1_dx(self%nnodes(),self%elem_info%nterms_s,self%elem_info%nterms_c,3), &
                  self%dgrad2_dx(self%nnodes(),self%elem_info%nterms_s,self%elem_info%nterms_c,3), &
                  self%dgrad3_dx(self%nnodes(),self%elem_info%nterms_s,self%elem_info%nterms_c,3)  ) 

        ! Loop through each donor node and compute differential gradient interpolators
        do ipt = 1,self%nnodes()
            ! Donor coordiantes(xi,eta,zeta)
            coord = self%coords(ipt,:)

            ! Evaluate basis modes at node location
            do iterm = 1,self%elem_info%nterms_c
                val_c(iterm)    =  polynomial_val(spacedim,self%elem_info%nterms_c,iterm,coord)
                ddxi_c(iterm)   = dpolynomial_val(spacedim,self%elem_info%nterms_c,iterm,coord,XI_DIR  )
                ddeta_c(iterm)  = dpolynomial_val(spacedim,self%elem_info%nterms_c,iterm,coord,ETA_DIR )
                ddzeta_c(iterm) = dpolynomial_val(spacedim,self%elem_info%nterms_c,iterm,coord,ZETA_DIR)
            end do

            ! Evaluate mesh point from dot product of modes and polynomial values
            dXdxi(1,1) = dot_product(ddxi_c,   modal_coords(:,1))
            dXdxi(1,2) = dot_product(ddeta_c,  modal_coords(:,1))
            dXdxi(1,3) = dot_product(ddzeta_c, modal_coords(:,1))
                                             
            dXdxi(2,1) = dot_product(ddxi_c,   modal_coords(:,2))
            dXdxi(2,2) = dot_product(ddeta_c,  modal_coords(:,2))
            dXdxi(2,3) = dot_product(ddzeta_c, modal_coords(:,2))
                                             
            dXdxi(3,1) = dot_product(ddxi_c,   modal_coords(:,3))
            dXdxi(3,2) = dot_product(ddeta_c,  modal_coords(:,3))
            dXdxi(3,3) = dot_product(ddzeta_c, modal_coords(:,3))


            !
            ! Compute coordinate derivatives of jacobian matrix wrt grid ndoes 
            ! at interpolation nodes
            !
            !nnodes_r = self%nnodes_r
            nnodes_r = self%elem_info%nterms_c
            dmodes   = self%nodes_to_modes
            if (allocated(djacobian_dx)) deallocate(djacobian_dx)
            allocate(djacobian_dx(3,3,nnodes_r,3), stat=ierr)
            if (ierr /= 0) call AllocationError
            
            !
            ! Initialize djacobian_dx with zeros
            !
            djacobian_dx = ZERO
            
            do idiff_n = 1,nnodes_r
                do icoord = 1,3
                    
                    djacobian_dx(icoord,1,idiff_n,icoord) = dot_product(ddxi_c,   dmodes(:,idiff_n))
                    djacobian_dx(icoord,2,idiff_n,icoord) = dot_product(ddeta_c,  dmodes(:,idiff_n))
                    djacobian_dx(icoord,3,idiff_n,icoord) = dot_product(ddzeta_c, dmodes(:,idiff_n))
                    
                end do
            end do
            
            
            !
            ! Apply transformation to second row of dJacobian/d1
            !
            if (self%elem_info%coordinate_system == CYLINDRICAL) then
                do idiff_n = 1,nnodes_r
                    djacobian_dx(2,1,idiff_n,1) = dot_product(val_c,dmodes(:,idiff_n)) * dXdxi(2,1)
                    djacobian_dx(2,2,idiff_n,1) = dot_product(val_c,dmodes(:,idiff_n)) * dXdxi(2,2)
                    djacobian_dx(2,3,idiff_n,1) = dot_product(val_c,dmodes(:,idiff_n)) * dXdxi(2,3)
                end do
            end if


            !
            ! Apply coordinate system scaling
            !
            if (scale_metric) then
                select case(self%elem_info%coordinate_system)
                    case(CYLINDRICAL) 
                        select case(frame_selector)
                            case('Undeformed')
                                r = dot_product(self%modal_coords(:,1),val_c)
                            case('Deformed')
                                !r = dot_product(self%modal_coords_def(:,1),val_c)
                                call write_line("update_interpolations_dx: deformed not yet implemented")
                        end select
                        dXdxi(2,:) = dXdxi(2,:) * r
                        do idiff_n = 1,nnodes_r
                            ! Apply the same transformation to djacobian_dx
                            djacobian_dx(2,1,idiff_n,2) = djacobian_dx(2,1,idiff_n,2)*r
                            djacobian_dx(2,2,idiff_n,2) = djacobian_dx(2,2,idiff_n,2)*r
                            djacobian_dx(2,3,idiff_n,2) = djacobian_dx(2,3,idiff_n,2)*r
                        end do
                end select
            end if



            !
            ! Invert jacobian matrix at interpolation node
            !
            if (allocated(dmetric_dx)) deallocate(dmetric_dx)
            allocate(dmetric_dx(3,3,nnodes_r,3), stat=ierr)
            dmetric_dx = ZERO
            do idiff_n = 1,nnodes_r
                do icoord = 1,3
                    dmetric_dx(:,:,idiff_n,icoord) = &
                    dinv_3x3(dXdxi,djacobian_dx(:,:,idiff_n,icoord))
                end do !icoord
            end do !idiff_n


            


            !
            ! Compute differential gradient interpolators
            !
            do iterm = 1,self%elem_info%nterms_s


                !
                ! Compute interpolators
                !
                ddxi_s   = dpolynomial_val(spacedim,self%elem_info%nterms_s,iterm,coord,XI_DIR  )
                ddeta_s  = dpolynomial_val(spacedim,self%elem_info%nterms_s,iterm,coord,ETA_DIR )
                ddzeta_s = dpolynomial_val(spacedim,self%elem_info%nterms_s,iterm,coord,ZETA_DIR)

                !do idiff_n = 1,self%nnodes_r
                do idiff_n = 1,self%elem_info%nterms_c
                    do idir = 1,3
                        
                        ! Compute physical derivative interpolators for gq node
                        self%dgrad1_dx(ipt,iterm,idiff_n,idir) = dmetric_dx(1,1,idiff_n,idir) * ddxi_s   + &
                                                                 dmetric_dx(2,1,idiff_n,idir) * ddeta_s  + &
                                                                 dmetric_dx(3,1,idiff_n,idir) * ddzeta_s 
                        self%dgrad2_dx(ipt,iterm,idiff_n,idir) = dmetric_dx(1,2,idiff_n,idir) * ddxi_s   + &
                                                                 dmetric_dx(2,2,idiff_n,idir) * ddeta_s  + &
                                                                 dmetric_dx(3,2,idiff_n,idir) * ddzeta_s 
                        self%dgrad3_dx(ipt,iterm,idiff_n,idir) = dmetric_dx(1,3,idiff_n,idir) * ddxi_s   + &
                                                                 dmetric_dx(2,3,idiff_n,idir) * ddeta_s  + &
                                                                 dmetric_dx(3,3,idiff_n,idir) * ddzeta_s 
                    
                    end do !idir
                end do !idiff_n

            end do !iterm

        end do ! ipt

    end subroutine update_interpolations_dx
    !****************************************************************







    !>  Release allocated memory for derivatives of gradient interpolators
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/14/2018
    !!
    !----------------------------------------------------------------
    subroutine release_interpolations_dx(self)
        class(chimera_donor_t), intent(inout) :: self
        
        !
        ! Deallocate differential gradient interpolator 
        !
        if (allocated(self%dgrad1_dx)) then
            deallocate(self%dgrad1_dx, self%dgrad2_dx, self%dgrad3_dx)
        end if

    end subroutine release_interpolations_dx
    !****************************************************************














end module type_chimera_donor
