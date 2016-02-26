module type_chidg_data
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use type_point,                     only: point_t
    use type_dict,                      only: dict_t

    ! Primary chidg_data_t components
    use type_domaininfo,                only: domaininfo_t
    use type_mesh,                      only: mesh_t
    use type_bcset,                     only: bcset_t
    use type_equationset_wrapper,       only: equationset_wrapper_t
    use type_solverdata,                only: solverdata_t

    !
    use type_equationset_function_data, only: equationset_function_data_t
    use type_bcset_coupling,            only: bcset_coupling_t

    ! Factory methods
    use mod_equations,                  only: create_equationset
    use mod_bc,                         only: create_bc

    ! Classes
    use type_bc,                        only: bc_t

    implicit none




    !> Container for ChiDG data.
    !!
    !!  The format here is to have arrays of mesh_t, bcset_t, and eqnset_t components. The 
    !!  index of those arrays corresponds to a domain in the local ChiDG environment. A 
    !!  solverdata_t component holds chidgVector_t and chidgMatrix_t components that are
    !!  initialized from the domain components and are informed of the number of domains
    !!  in addition to their dependencies on each other.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------------------------
    type, public  :: chidg_data_t

        logical                                     :: solverInitialized = .false.
        integer(ik),        private                 :: ndomains_ = 0
        type(domaininfo_t),             allocatable :: info(:)      !< General container for domain information

        
        type(mesh_t),                   allocatable :: mesh(:)      !< Array of mesh instances. One for each domain.
        type(bcset_t),                  allocatable :: bcset(:)     !< Array of boundary condition set instances. One for each domain.
        type(equationset_wrapper_t),    allocatable :: eqnset(:)    !< Array of equation set instances. One for each domain.
        type(solverdata_t)                          :: sdata        !< Solver data container for solution vectors and matrices

    contains

        ! Initialization procedure for solution data. Execute after all domains are added.
        procedure   :: init_sdata

        ! Modifiers for adding domains and boundary conditions
        procedure   :: add_domain
        procedure   :: add_bc

        ! Accessors
        procedure   :: get_domain_index     !< Given a domain name, return domain index
        procedure   :: ndomains             !< Return number of domains in chidg instance

    end type chidg_data_t
    !*************************************************************************************************************




contains


    !> Initialize solution data storage structures. Needs to be called before accessing any 
    !! solution storage containers so they are allocated and initialized. 
    !!
    !! All domains must be added before calling this procedure, since the array of mesh 
    !! structures are used for the initialization routine.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    ! subroutine initialize_solver_data
    subroutine init_sdata(self)
        class(chidg_data_t),     intent(inout)   :: self

        integer(ik) :: idom, ndom, ierr

        type(equationset_function_data_t),  allocatable :: function_data(:)
        type(bcset_coupling_t),             allocatable :: bcset_coupling(:)


        !
        ! Assemble array of function_data from the eqnset array to pass to the solver data structure for 
        ! initialization
        !
        ndom = self%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%ndomains()
            function_data(idom) = self%eqnset(idom)%item%function_data
        end do


        !
        ! Assemble boundary condition coupling information to pass to sdata initialization for LHS storage
        !
        ndom = self%ndomains()
        allocate(bcset_coupling(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%ndomains()
            bcset_coupling(idom) = self%bcset(idom)%get_bcset_coupling()
        end do



        !
        ! Initialize solver data 
        !
        call self%sdata%init(self%mesh, bcset_coupling, function_data)

    end subroutine init_sdata
    !*************************************************************************************************************








    !> Add a domain to ChiDG. Calls initialization routines for components which define a 
    !! domain. That is, a mesh_t, an equationset_t, and the number of terms in their 
    !! polynomial expansions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  points      point_t matrix defining the mesh
    !!  @param[in]  nterms_c    Integer defining the number of terms in the element coordinate expansion
    !!  @param[in]  eqnset      Character string defining the equationset_t for the domain
    !!  @param[in]  nterms_s    Integer defining the number of terms in the solution expansion
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine add_domain(self,name,points,nterms_c,eqnset,nterms_s)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: name
        type(point_t),          intent(in)      :: points(:,:,:)
        integer(ik),            intent(in)      :: nterms_c
        character(*),           intent(in)      :: eqnset
        integer(ik),            intent(in)      :: nterms_s

        integer(ik) :: idom, ierr


        type(domaininfo_t),             allocatable :: temp_info(:)
        type(mesh_t),                   allocatable :: temp_mesh(:)
        type(bcset_t),                  allocatable :: temp_bcset(:)
        type(equationset_wrapper_t),    allocatable :: temp_eqnset(:)


        !
        ! Increment number of domains by one
        !
        self%ndomains_ = self%ndomains_ + 1
        idom = self%ndomains_



        !
        ! Resize array storage
        !
        ! Allocate new storage arrays
        allocate(temp_info(self%ndomains_),   &
                 temp_mesh(self%ndomains_),   &
                 temp_bcset(self%ndomains_),  &
                 temp_eqnset(self%ndomains_), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Copy previously initialized instances to new array. Be careful about pointers components here!
        ! For example, a pointer from a face to an element would no longer be valid in the new array.
        if (self%ndomains_ > 1) then
            !temp_mesh(   1:size(self%mesh))    = self%mesh     ! ifort segfaults on this for cases with sevaral domains
            !temp_bcset(  1:size(self%bcset))   = self%bcset
            !temp_eqnset( 1:size(self%eqnset))  = self%eqnset
            temp_info(   1:size(self%info))    = self%info(1:size(self%mesh))
            temp_mesh(   1:size(self%mesh))    = self%mesh(1:size(self%mesh))
            temp_bcset(  1:size(self%bcset))   = self%bcset(1:size(self%mesh))
            temp_eqnset( 1:size(self%eqnset))  = self%eqnset(1:size(self%mesh))
        end if


        !
        ! Set domain info
        !
        temp_info(idom)%name = name


        !
        ! Initialize new mesh
        !
        call temp_mesh(idom)%init_geom(idom,nterms_c,points)


        !
        ! Allocate equation set
        !
        call create_equationset(eqnset,temp_eqnset(idom)%item)


        !
        ! Initialize mesh numerics based on equation set and polynomial expansion order
        !
        call temp_mesh(idom)%init_sol(temp_eqnset(idom)%item%neqns,nterms_s)


        !
        ! Move rezied temp allocation back to chidg_data container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_info,self%info)
        call move_alloc(temp_mesh,self%mesh)
        call move_alloc(temp_bcset,self%bcset)
        call move_alloc(temp_eqnset,self%eqnset)

    end subroutine add_domain
    !***************************************************************************************************************













    !> Add a boundary condition to a ChiDG domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  domain      Character string of the selected domain
    !!  @param[in]  bc          Character string indicating the boundary condition to add
    !!  @param[in]  face        Integer of the block face to which the boundary condition will be applied
    !!  @param[in]  options     Boundary condition options dictionary
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine add_bc(self,domain,bc,face)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: domain
        class(bc_t),            intent(inout)   :: bc
        integer(ik),            intent(in)      :: face

        integer(ik)                 :: idom, ierr
        class(bc_t), allocatable    :: bc_copy


        !
        ! Get domain index from domain string
        !
        idom = self%get_domain_index(domain)


        !
        ! Create a copy of the incoming boundary condition. The copy can then be initialized with a domain geometry
        !
        allocate(bc_copy, source=bc, stat=ierr)
        if ( ierr /= 0 ) call AllocationError


        !
        ! Initialize new boundary condition from mesh data and face index
        !
        call bc_copy%init(self%mesh(idom),face)


        !
        ! Add initialized boundary condition to bcset_t for domain 'idom'
        !
        call self%bcset(idom)%add(bc_copy)


    end subroutine add_bc
    !**********************************************************************************************************











    !> Given a character string corresponding to the name of a given domain,
    !! find and return the index of that domain in the ChiDG_data instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  dname           String associated with a given domain
    !!  @return     domain_index    Integer index of the associated domain
    !!
    !----------------------------------------------------------------------------------------------------------
    function get_domain_index(self,dname) result(domain_index)
        class(chidg_data_t),    intent(in)      :: self
        character(*),           intent(in)      :: dname

        integer(ik)  :: idom
        integer(ik)  :: domain_index
        
        domain_index = 0

        do idom = 1,self%ndomains_

            !
            ! Test name
            !
            if ( trim(dname) == trim(self%info(idom)%name) ) then
                domain_index = idom
                exit
            end if

        end do


        if (domain_index == 0) call chidg_signal(FATAL,"chidg_data%get_domain_index :: no domain found matching given name")

    end function get_domain_index
    !**********************************************************************************************************










    !> Return the number of domains in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    function ndomains(self) result(ndom)
        class(chidg_data_t),    intent(in)      :: self

        integer :: ndom

        ndom = self%ndomains_

    end function ndomains
    !**********************************************************************************************************









end module type_chidg_data
