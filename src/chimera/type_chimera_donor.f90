module type_chimera_donor
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NO_ID
    implicit none



    !>  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !!
    !---------------------------------------------------------------
    type, public :: chimera_donor_t

        ! Donor location
        integer(ik)                 :: idomain_g
        integer(ik)                 :: idomain_l
        integer(ik)                 :: ielement_g
        integer(ik)                 :: ielement_l
        integer(ik)                 :: iproc            ! donor processor rank

        ! Donor properties
        integer(ik)                 :: nfields   = 0        ! Number of equations in donor element
        integer(ik)                 :: ntime     = 0        ! Number of equations in donor element
        integer(ik)                 :: nterms_s  = 0        ! Number of terms in donor expansion
        integer(ik)                 :: nterms_c  = 0        ! Number of terms in donor expansion
        integer(ik)                 :: dof_start = 0        ! Starting gobal dof index for overset donor.
        integer(ik)                 :: dof_local_start = 0  ! Starting local dof index for overset donor
        integer(ik)                 :: eqn_ID    = NO_ID    ! Equation set identifier


        ! Parallel access information
        integer(ik)                 :: pelem_ID     = NO_ID ! ID in mesh%parallel_elements(pelem_ID) (only if off-processor)
        integer(ik)                 :: recv_comm    = NO_ID ! location of donor solution in native storage
        integer(ik)                 :: recv_domain  = NO_ID ! location of donor solution in native storage
        integer(ik)                 :: recv_element = NO_ID ! location of donor solution in native storage
        integer(ik)                 :: recv_dof     = NO_ID ! location of donor solution in petsc storage


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


    contains

        procedure   :: set_properties

        procedure   :: add_node
        procedure   :: nnodes

    
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
    function chimera_donor(idomain_g, idomain_l, ielement_g, ielement_l, iproc) result(instance)
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_g
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iproc

        type(chimera_donor_t)   :: instance

        instance%idomain_g  = idomain_g
        instance%idomain_l  = idomain_l
        instance%ielement_g = ielement_g
        instance%ielement_l = ielement_l
        instance%iproc      = iproc

    end function chimera_donor
    !******************************************************************




    !>  Set the properties of the donor.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !------------------------------------------------------------------
    subroutine set_properties(self,nterms_c,nterms_s,ntime,nfields,eqn_ID,dof_start,dof_local_start)
        class(chimera_donor_t), intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms_c
        integer(ik),            intent(in)      :: nterms_s
        integer(ik),            intent(in)      :: ntime
        integer(ik),            intent(in)      :: nfields
        integer(ik),            intent(in)      :: eqn_ID
        integer(ik),            intent(in)      :: dof_start
        integer(ik),            intent(in)      :: dof_local_start

        self%nterms_c        = nterms_c
        self%nterms_s        = nterms_s
        self%ntime           = ntime
        self%nfields         = nfields
        self%eqn_ID          = eqn_ID
        self%dof_start       = dof_start
        self%dof_local_start = dof_local_start

    end subroutine set_properties
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




end module type_chimera_donor
