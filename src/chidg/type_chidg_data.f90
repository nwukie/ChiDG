module type_chidg_data
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use type_point,                     only: point_t
    use type_dict,                      only: dict_t
    use type_domain_connectivity,       only: domain_connectivity_t
    use type_boundary_connectivity,     only: boundary_connectivity_t

    ! Primary chidg_data_t components
    use type_domain_info,               only: domain_info_t
    use type_mesh,                      only: mesh_t
    use type_bcset,                     only: bcset_t
    use type_bc,                        only: bc_t
    use type_bc_state,                  only: bc_state_t
    use type_bcvector,                  only: bcvector_t
    use type_bc_group,                  only: bc_group_t
    use type_svector,                   only: svector_t
    use mod_string,                     only: string_t
    use type_equation_set,              only: equation_set_t
    use type_solverdata,                only: solverdata_t
    use type_time_manager,              only: time_manager_t

    use type_equationset_function_data, only: equationset_function_data_t
    use type_bcset_coupling,            only: bcset_coupling_t

    ! Factory methods
    use mod_equations,                  only: equation_builder_factory


    implicit none




    !>  Container for ChiDG data.
    !!
    !!  The format here is to have arrays of mesh_t, bcset_t, and eqnset_t components. The 
    !!  index of those arrays corresponds to a domain in the local ChiDG environment. A 
    !!  solverdata_t component holds chidg_vector_t and chidg_matrix_t components that are
    !!  initialized from the domain components and are informed of the number of domains
    !!  in addition to their dependencies on each other.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public  :: chidg_data_t

        logical                                     :: solverInitialized = .false.
        integer(ik),        private                 :: ndomains_ = 0
        !integer(ik),        private                 :: ntime_    = 1   !now part of time_manager
        integer(ik),        private                 :: spacedim_ = 3    !< Default 3D 

        
        ! For each domain: info, a mesh, a boundary condition set, and an equation set
        type(domain_info_t),            allocatable :: info(:)     
        type(mesh_t),                   allocatable :: mesh(:)     
        type(bcset_t),                  allocatable :: bcset(:)    
        type(equation_set_t),           allocatable :: eqnset(:)   

        ! An object containing matrix and vector storage
        type(solverdata_t)                          :: sdata

        ! An object containing time information
        type(time_manager_t)                        :: time_manager


    contains

        ! Modifiers for adding domains and boundary conditions
        procedure   :: add_domain
        procedure   :: add_bc

        ! Initialization procedure for solution data. Execute after all domains are added.
        procedure   :: initialize_solution_domains
        procedure   :: initialize_solution_solver

        ! Accessors
        procedure   :: get_domain_index             !< Given a domain name, return domain index
        procedure   :: ndomains                     !< Return number of domains in chidg instance
        procedure   :: ntime
        procedure   :: get_dimensionality
        procedure   :: get_auxiliary_field_names    !< Return required auxiliary fields

        procedure   :: report

    end type chidg_data_t
    !*******************************************************************************************




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
    !------------------------------------------------------------------------------------------
    subroutine initialize_solution_solver(self)
        class(chidg_data_t),     intent(inout)   :: self

        integer(ik) :: idom, ndom, ierr

        type(equationset_function_data_t),  allocatable :: function_data(:)
        type(bcset_coupling_t),             allocatable :: bcset_coupling(:)


        !
        ! Assemble array of function_data from the eqnset array to pass to the solver data 
        ! structure for initialization
        !
        ndom = self%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%ndomains()
            function_data(idom) = self%eqnset(idom)%function_data
        end do


        !
        ! Assemble boundary condition coupling information to pass to sdata initialization 
        ! for LHS storage
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

    end subroutine initialize_solution_solver
    !******************************************************************************************








    !> Add a domain to ChiDG. Calls initialization routines for components which define a 
    !! domain. That is, a mesh_t, an equationset_t, and the number of terms in their 
    !! polynomial expansions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  points      point_t matrix defining the mesh
    !!  @param[in]  nterms_c    Integer defining the number of terms in the element coordinate 
    !!                          expansion
    !!  @param[in]  eqnset      Character string defining the equationset_t for the domain
    !!  @param[in]  nterms_s    Integer defining the number of terms in the solution expansion
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_domain(self,name,nodes,connectivity,spacedim,nterms_c,eqnset,coord_system)
        class(chidg_data_t),            intent(inout)   :: self
        character(*),                   intent(in)      :: name
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: spacedim
        integer(ik),                    intent(in)      :: nterms_c
        character(*),                   intent(in)      :: eqnset
        character(*),                   intent(in)      :: coord_system

        integer(ik)                 :: idomain_l, ierr, idom
        character(:),   allocatable :: user_msg


        type(domain_info_t),            allocatable :: temp_info(:)
        type(mesh_t),                   allocatable :: temp_mesh(:)
        type(bcset_t),                  allocatable :: temp_bcset(:)
        type(equation_set_t),           allocatable :: temp_eqnset(:)



        !
        ! Increment number of domains by one
        !
        self%ndomains_ = self%ndomains_ + 1
        idomain_l      = self%ndomains_


        !
        ! Resize array storage
        !
        allocate( &
                 temp_info(self%ndomains_),   &
                 temp_mesh(self%ndomains_),   &
                 temp_bcset(self%ndomains_),  &
                 temp_eqnset(self%ndomains_), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%ndomains_ > 1) then
            temp_info(   1:size(self%info))    = self%info(1:size(self%mesh))
            temp_mesh(   1:size(self%mesh))    = self%mesh(1:size(self%mesh))
            temp_bcset(  1:size(self%bcset))   = self%bcset(1:size(self%mesh))
            temp_eqnset( 1:size(self%eqnset))  = self%eqnset(1:size(self%mesh))
        end if


        !
        ! Set domain info
        !
        temp_info(idomain_l)%name = name


        !
        ! Initialize new mesh
        !
        call temp_mesh(idomain_l)%init_geom(idomain_l,spacedim,nterms_c,nodes,connectivity,coord_system)


        !
        ! Check that a domain with the same global index wasn't already added. For example, if 
        ! a block got split and put on the same processor. Some of the MPI communication assumes 
        ! one unique global domain index for each domain on the processor.
        !
        user_msg = "chidg_data%add_domain: Two domains have the same global index. MPI &
                    communication assumes this does not happen."
        if (self%ndomains_ > 1) then
            do idom = 1,size(self%mesh)
                if (self%mesh(idom)%idomain_g == temp_mesh(idomain_l)%idomain_g) call chidg_signal(FATAL,user_msg)
            end do !idom
        end if


        !
        ! Allocate equation set
        !
        temp_eqnset(idomain_l) = equation_builder_factory%produce(eqnset,'default')



        !
        ! Move resized temp allocation back to chidg_data container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_info,self%info)
        call move_alloc(temp_mesh,self%mesh)
        call move_alloc(temp_bcset,self%bcset)
        call move_alloc(temp_eqnset,self%eqnset)


    end subroutine add_domain
    !*******************************************************************************************











    !>  For a ChiDG domain, add a boundary condition patche and associate it with a boundary 
    !!  condition group.
    !!
    !!
    !!  Boundary condition groups hold sets of state functions that are used to compute an 
    !!  exterior state on the boundary. The boundary condition groups are defined for the 
    !!  global problem. Here, the individual patches of a given domain are being set to a 
    !!  specific group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  domain          Character string of the selected domain.
    !!  @param[in]  bc_connectivity Face connectivities defining the boundary condition patch.
    !!  @param[in]  bc_group        Name of the boundary condition group to associate with the 
    !!                              patch.
    !!  @param[in]  bc_groups       bc_group_t's for the global problem that can be searched 
    !!                              through and used to initialize.
    !!
    !!  To force a particular bc_state on a boundary condition, one can pass a bc_state_t in 
    !!  as an option for bc_wall, bc_inlet, bc_outlet, bc_symmetry
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc(self,domain,bc_connectivity,bc_group,bc_groups,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)
        class(chidg_data_t),            intent(inout)           :: self
        character(*),                   intent(in)              :: domain
        type(boundary_connectivity_t),  intent(in)              :: bc_connectivity
        character(*),                   intent(in)              :: bc_group
        type(bc_group_t),               intent(in)              :: bc_groups(:)
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic


        character(:),       allocatable     :: user_msg
        class(bc_state_t),  allocatable     :: bc_state
        type(bc_t)                          :: bc
        integer(ik)                         :: idom, BC_ID, istate, igroup, ierr
        logical                             :: group_found, group_set


        !
        ! Get domain index from domain string
        !
        idom = self%get_domain_index(domain)


        !
        ! Add a new boundary condition and get ID
        !
        BC_ID = self%bcset(idom)%add(bc)



        !
        ! Initialize new boundary condition from mesh data and connectivity information.
        ! NOTE: init_bc needs called after the boundary condition has been added to the 
        !       set so it can inform the mesh about it's BC_ID.
        !
        call self%bcset(idom)%bcs(BC_ID)%init_bc(self%mesh(idom),bc_connectivity,bc_group,bc_groups, &
                                                                                 bc_wall,            &
                                                                                 bc_inlet,           &
                                                                                 bc_outlet,          &
                                                                                 bc_symmetry,        &
                                                                                 bc_farfield,        &
                                                                                 bc_periodic)


    end subroutine add_bc
    !******************************************************************************************








    !>  For each domain, call solution initialization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine initialize_solution_domains(self,nterms_s)
        class(chidg_data_t),    intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms_s

        integer(ik) :: idomain, neqns

        !self%ntime_ = ntime    !we do not need to pass ntime in, since ntime is contained in data%time_manager%ntime

        ! Initialize mesh numerics based on equation set and polynomial expansion order
        do idomain = 1,self%ndomains()
            neqns = self%eqnset(idomain)%prop%nprimary_fields()
            call self%mesh(idomain)%init_sol(neqns,nterms_s,self%time_manager%ntime)
        end do


    end subroutine initialize_solution_domains
    !******************************************************************************************









    !> Given a character string corresponding to the name of a given domain,
    !! find and return the index of that domain in the ChiDG_data instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  domain_name     String associated with a given domain
    !!  @return     domain_index    Integer index of the associated domain
    !!
    !-------------------------------------------------------------------------------------------
    function get_domain_index(self,domain_name) result(domain_index)
        class(chidg_data_t),    intent(in)      :: self
        character(*),           intent(in)      :: domain_name

        character(:),   allocatable :: user_msg
        integer(ik)  :: idom
        integer(ik)  :: domain_index
        
        domain_index = 0

        do idom = 1,self%ndomains_
            if ( trim(domain_name) == trim(self%info(idom)%name) ) then
                domain_index = idom
                exit
            end if
        end do


        user_msg = "chidg_data%get_domain_index: No domain was found that had a name &
                    that matched the incoming string"
        if (domain_index == 0) call chidg_signal_one(FATAL,user_msg,domain_name)

    end function get_domain_index
    !*******************************************************************************************




    !> Return the number of domains in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function ndomains(self) result(ndom)
        class(chidg_data_t),    intent(in)      :: self

        integer :: ndom

        ndom = self%ndomains_

    end function ndomains
    !*******************************************************************************************







    !>  Return the dimensionality of the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !-------------------------------------------------------------------------------------------
    function get_dimensionality(self) result(dimensionality)
        class(chidg_data_t),    intent(in)      :: self

        integer :: dimensionality

        dimensionality = self%spacedim_

    end function get_dimensionality
    !*******************************************************************************************












    !>  Return a vector of auxiliary fields that are required.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !-------------------------------------------------------------------------------------------
    function get_auxiliary_field_names(self) result(field_names)
        class(chidg_data_t),    intent(in)  :: self

        integer(ik)                 :: idom, ifield
        type(svector_t)             :: field_names
        character(:),   allocatable :: field_name



        do idom = 1,self%ndomains()
            do ifield = 1,self%eqnset(idom)%prop%nauxiliary_fields()

                field_name = self%eqnset(idom)%prop%get_auxiliary_field_name(ifield)
                call field_names%push_back_unique(string_t(field_name))

            end do !ifield
        end do !idom



    end function get_auxiliary_field_names
    !*******************************************************************************************







    !> Return the number of time instances in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function ntime(self) result(ndom)
        class(chidg_data_t),    intent(in)      :: self

        integer :: ndom

        ndom = self%time_manager%ntime

    end function ntime
    !*******************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/7/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine report(self,selection)
        class(chidg_data_t),    intent(in)  :: self
        character(*),           intent(in)  :: selection

        integer(ik) :: idom


        if ( trim(selection) == 'grid' ) then

            do idom = 1,self%ndomains()
                call write_line('Domain ', idom, '  :  ', self%mesh(idom)%nelem, ' Elements', io_proc=IRANK)
            end do


        else


        end if



    end subroutine report
    !*******************************************************************************************






end module type_chidg_data
