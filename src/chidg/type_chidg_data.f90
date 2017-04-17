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
    use type_bc,                        only: bc_t
    use type_bc_state,                  only: bc_state_t
    use type_bc_group,                  only: bc_group_t
    use type_bcvector,                  only: bcvector_t
    use type_svector,                   only: svector_t
    use mod_string,                     only: string_t
    use type_equation_set,              only: equation_set_t
    use type_solverdata,                only: solverdata_t
    use type_time_manager,              only: time_manager_t

    use type_equationset_function_data, only: equationset_function_data_t
    use type_bcset_coupling,            only: bcset_coupling_t

    ! Factory methods
    use mod_equations,                  only: equation_builder_factory

    !Mesh motion
    use type_prescribed_mesh_motion,               only: prescribed_mesh_motion_t
    use type_prescribed_mesh_motion_group,               only: prescribed_mesh_motion_group_t


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
        integer(ik),        private                 :: spacedim_ = 3

        
        ! For each domain: info, a mesh, and an equation set
        type(domain_info_t),            allocatable :: info(:)     
        type(mesh_t),                   allocatable :: mesh(:)     

        ! Boundary conditions are not specified per-domain. 
        ! A boundary condition for each bc_group in the file.
        type(bc_t),                     allocatable :: bc(:)    
        type(equation_set_t),           allocatable :: eqnset(:)

        ! An object containing matrix and vector storage
        type(solverdata_t)                          :: sdata

        ! An object containing time information
        type(time_manager_t)                        :: time_manager

        ! For each domain, a mesh motion type
        type(prescribed_mesh_motion_t),            allocatable :: pmm(:)

    contains

        ! Modifiers for adding domains and boundary conditions
        procedure   :: add_domain
        procedure   :: add_bc_group
        procedure   :: add_bc_patch
        procedure   :: new_bc
        procedure   :: new_equation_set

        ! Initialization procedure for solution data. Execute after all domains are added.
        procedure   :: initialize_solution_domains
        procedure   :: initialize_solution_bc
        procedure   :: initialize_solution_solver

        ! Mesh Motion
        procedure   :: new_pmm
        procedure   :: add_pmm_domain
        procedure   :: add_pmm_group


        ! Release allocated memory
        procedure   :: release

        ! Accessors
        procedure   :: get_domain_index             ! Given a domain name, return domain index
        procedure   :: ndomains                     ! Return number of allocated domains
        procedure   :: nbcs                         ! Return number of allocated boundary conditions
        procedure   :: nequation_sets               ! Return number of allocated equation sets
        procedure   :: ntime
        procedure   :: get_dimensionality
        procedure   :: get_auxiliary_field_names    ! Return required auxiliary fields

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

        integer(ik) :: idom, ndom, ierr, eqn_ID

        type(equationset_function_data_t),  allocatable :: function_data(:)
        type(bcset_coupling_t),             allocatable :: bcset_coupling(:)


        call write_line("Matrices/Vectors: allocating...", io_proc=GLOBAL_MASTER)
        !
        ! Assemble array of function_data from the eqnset array to pass to the solver data 
        ! structure for initialization
        !
        ndom = self%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%ndomains()
            eqn_ID = self%mesh(idom)%eqn_ID
            function_data(idom) = self%eqnset(eqn_ID)%function_data
        end do


!        !
!        ! Assemble boundary condition coupling information to pass to sdata initialization 
!        ! for LHS storage
!        !
!        ndom = self%ndomains()
!        allocate(bcset_coupling(ndom), stat=ierr)
!        if ( ierr /= 0 ) call AllocationError
!
!        do idom = 1,self%ndomains()
!            bcset_coupling(idom) = self%bcset(idom)%get_bcset_coupling()
!        end do


        !
        ! Initialize solver data 
        !
        !call self%sdata%init(self%mesh, bcset_coupling, function_data)
        !call self%sdata%init(self%mesh, self%bc, function_data)
        call self%sdata%init(self%mesh, function_data)

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
    subroutine add_domain(self,name,nodes,connectivity,nelements_g,spacedim,nterms_c,eqnset,coord_system)
        class(chidg_data_t),            intent(inout)   :: self
        character(*),                   intent(in)      :: name
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: nelements_g
        integer(ik),                    intent(in)      :: spacedim
        integer(ik),                    intent(in)      :: nterms_c
        character(*),                   intent(in)      :: eqnset
        character(*),                   intent(in)      :: coord_system

        integer(ik)                 :: idomain_l, ierr, idom, ieqn, eqn_ID
        logical                     :: already_added
        character(:),   allocatable :: user_msg


        type(domain_info_t),            allocatable :: temp_info(:)
        type(mesh_t),                   allocatable :: temp_mesh(:)
        type(equation_set_t),           allocatable :: temp_eqnset(:)


        !
        ! Increment location of new domain by one
        !
        idomain_l = self%ndomains() + 1




        !
        ! Check if equation set has already been added. If so, get identifier index eqn_ID
        !
        already_added = .false.
        do ieqn = 1,self%nequation_sets()
            already_added = (trim(self%eqnset(ieqn)%name) == trim(eqnset)) 
            if (already_added) eqn_ID = ieqn
            if (already_added) exit
        end do
        

        !
        ! Add new equation set if it doesn't already exist and get new eqn_ID
        !
        if (.not. already_added) then
            eqn_ID = self%new_equation_set()
            self%eqnset(eqn_ID) = equation_builder_factory%produce(eqnset, 'default')
            self%eqnset(eqn_ID)%eqn_ID = eqn_ID
        end if


        !
        ! Resize array storage
        !
        allocate( temp_info(  self%ndomains()+1),   &
                  temp_mesh(  self%ndomains()+1), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%ndomains() > 0) then
            temp_info(   1:size(self%info))    = self%info(1:size(self%info))
            temp_mesh(   1:size(self%mesh))    = self%mesh(1:size(self%mesh))
        end if


        !
        ! Set domain info
        !
        temp_info(idomain_l)%name = name


        !
        ! Initialize new mesh: geometry and equation set identifier
        !
        call temp_mesh(idomain_l)%init_geom(idomain_l,nelements_g,spacedim,nterms_c,nodes,connectivity,coord_system)
        call temp_mesh(idomain_l)%init_eqn(eqn_ID)

        !
        ! Check that a domain with the same global index wasn't already added. For example, if 
        ! a block got split and put on the same processor. Some of the MPI communication assumes 
        ! one unique global domain index for each domain on the processor.
        !
        user_msg = "chidg_data%add_domain: Two domains have the same global index. MPI &
                    communication assumes this does not happen."
        if (self%ndomains() > 1) then
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


    end subroutine add_domain
    !*******************************************************************************************











    !>  Create a new boundary condition and add the incoming
    !>  !For a ChiDG domain, add a boundary condition patch and associate it with a boundary 
    !!  !condition group.
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
    subroutine add_bc_group(self,bc_group,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic)
        class(chidg_data_t),            intent(inout)           :: self
        type(bc_group_t),               intent(in)              :: bc_group
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic


        integer(ik)                         :: bc_ID


        !
        ! Create a new boundary condition
        !
        bc_ID = self%new_bc()



        !
        ! Initialize boundary condition state functions from bc_group
        !
        call self%bc(bc_ID)%init_bc_group(bc_group, bc_wall,        &
                                                    bc_inlet,       &
                                                    bc_outlet,      &
                                                    bc_symmetry,    &
                                                    bc_farfield,    &
                                                    bc_periodic)



    end subroutine add_bc_group
    !******************************************************************************************













    !>  Add a bc_patch_t to the appropriate boundary condition.
    !!
    !!  NOTE: boundary condition groups must be added BEFORE this is called.
    !!
    !!  Searches through boundary condition groups that have already been added, self%bc,
    !!  searching for one with the correct name defined in the patch. Once found, initialize
    !!  the patch connectivity data on the boundary condition identified.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc_patch(self, domain_name, patch_bc_name, bc_connectivity)
        class(chidg_data_t),            intent(inout)   :: self
        character(*),                   intent(in)      :: domain_name
        character(*),                   intent(in)      :: patch_bc_name
        type(boundary_connectivity_t),  intent(in)      :: bc_connectivity

        character(:),   allocatable :: bc_name, user_msg
        integer(ik)                 :: bc_ID, ibc, idom
        logical                     :: found_bc

        
        !
        ! Find the correct boundary condition to add bc_patch to
        !
        do ibc = 1,size(self%bc)

            bc_name = self%bc(ibc)%get_name()
            found_bc = (trim(patch_bc_name) == trim(bc_name))

            if (found_bc) bc_ID = ibc
            if (found_bc) exit

        end do




        !
        ! Once bc is found, initialize bc_patch on bc
        !
        if (found_bc) then
            
            !
            ! Find domain index in mesh(:) from domain_name
            !
            idom = self%get_domain_index(domain_name)

            call self%bc(bc_ID)%init_bc_patch(self%mesh(idom), bc_connectivity)

        else


            user_msg = "chidg_data%add_bc_patch: It looks like we didn't find a boundary state &
                        group that matches with the string indicated in a boundary patch. Make &
                        sure that a boundary state group with the correct name exists. Also make &
                        sure that the name set on the boundary patch corresponds to one of the &
                        boundary state groups that exists."
            if ( (trim(patch_bc_name) /= 'empty') .and. &
                 (trim(patch_bc_name) /= 'Empty') ) &
                call chidg_signal_one(FATAL,user_msg,trim(patch_bc_name))

        end if


    end subroutine add_bc_patch
    !*******************************************************************************************








    

    !>  Extend the self%bc array to include another instance. Return the ID of the new
    !!  boundary condition where it can be found in the array as self%bc(bc_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_bc(self) result(bc_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(bc_t), allocatable :: temp_bcs(:)
        integer(ik)             :: bc_ID, ierr, nbc


        !
        ! Get number of boundary conditions
        !
        if (allocated(self%bc)) then
            nbc = size(self%bc)
        else
            nbc = 0
        end if


        !
        ! Increment number of boundary conditions
        !
        nbc = nbc + 1


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_bcs(nbc), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( nbc > 1) then
            temp_bcs(1:size(self%bc)) = self%bc(1:size(self%bc))
        end if


        !
        ! Set ID of new bc and store to array
        !
        bc_ID = nbc
        temp_bcs(bc_ID)%bc_ID = bc_ID


        !
        ! Attach extended allocation to chidg_data%bc
        !
        call move_alloc(temp_bcs,self%bc)



    end function new_bc
    !******************************************************************************************


    !
    !   Mesh Motion
    !------------------------------------------------------------------------------------------

    !> 
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_pmm_group(self,pmm_group)
        class(chidg_data_t),            intent(inout)           :: self
        type(prescribed_mesh_motion_group_t),               intent(inout)              :: pmm_group

        integer(ik)                         :: pmm_ID


        !
        ! Create a new boundary condition
        !
        pmm_ID = self%new_pmm()



        !
        ! Initialize pmm from pmm_group. Note that, since pmm_group contains pmm,
        ! we must pass the pmm info instead of pmm_group to avoid a circular dependency.
        !
        call self%pmm(pmm_ID)%init_pmm_group(pmm_group%pmm)


    end subroutine add_pmm_group
    !******************************************************************************************












    !>  Extend the self%eqnset array to include another instance. Return the ID of the new
    !!  equation set where it can be found in the array as self%eqnset(eqn_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_equation_set(self) result(eqn_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(equation_set_t), allocatable   :: temp_eqnset(:)
        integer(ik)                         :: eqn_ID, ierr, neqnset


        !
        ! Get number of boundary conditions
        !
        neqnset = self%nequation_sets()


        !
        ! Increment number of boundary conditions
        !
        neqnset = neqnset + 1


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_eqnset(neqnset), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( neqnset > 1) then
            temp_eqnset(1:size(self%eqnset)) = self%eqnset(1:size(self%eqnset))
        end if


        !
        ! Set ID of new bc and store to array
        !
        eqn_ID = neqnset
        temp_eqnset(eqn_ID)%eqn_ID = eqn_ID


        !
        ! Attach extended allocation to chidg_data%eqnset
        !
        call move_alloc(temp_eqnset,self%eqnset)



    end function new_equation_set
    !******************************************************************************************

    !> 
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_pmm_domain(self, domain_name, domain_pmm_name)
        class(chidg_data_t),            intent(inout)   :: self
        character(*),                   intent(in)      :: domain_name
        character(*),                   intent(in)      :: domain_pmm_name

        character(:),   allocatable :: pmm_name, user_msg
        integer(ik)                 :: pmm_ID, ipmm, idom
        logical                     :: found_pmm

        
        !
        ! Find the correct boundary condition to add bc_patch to
        !
        do ipmm = 1,size(self%pmm)

            pmm_name = self%pmm(ipmm)%get_name()
            found_pmm = (trim(domain_pmm_name) == trim(pmm_name))

            if (found_pmm) pmm_ID = ipmm
            if (found_pmm) exit

        end do




        !
        ! Once pmm is found, initialize pmm_patch on pmm
        !
        if (found_pmm) then
            
            !
            ! Find domain index in mesh(:) from domain_name
            !
            idom = self%get_domain_index(domain_name)

            call self%pmm(pmm_ID)%init_pmm_domain(self%mesh(idom))

        else


            user_msg = "chidg_data%add_pmm_patch: It looks like we didn't find a prescribed mesh motion  &
                        group that matches with the string indicated in a pmm domain. Make &
                        sure that a pmm group with the correct name exists. Also make &
                        sure that the name set on the pmm domain corresponds to one of the &
                        pmm groups that exists."
            if ( (trim(domain_pmm_name) /= 'empty') .and. &
                (trim(domain_pmm_name) /= 'Empty') ) &
                call chidg_signal_one(FATAL,user_msg,trim(domain_pmm_name))

        end if


    end subroutine add_pmm_domain
    !*******************************************************************************************
   

    !>
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_pmm(self) result(pmm_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(prescribed_mesh_motion_t), allocatable :: temp_pmms(:)
        integer(ik)             :: pmm_ID, ierr, npmm


        !
        ! Get number of boundary conditions
        !
        if (allocated(self%pmm)) then
            npmm = size(self%pmm)
        else
            npmm = 0
        end if


        !
        ! Increment number of boundary conditions
        !
        npmm = npmm + 1


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_pmms(npmm), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( npmm > 1) then
            temp_pmms(1:size(self%pmm)) = self%pmm(1:size(self%pmm))
        end if


        !
        ! Set ID of new pmm and store to array
        !
        pmm_ID = npmm
        temp_pmms(pmm_ID)%pmm_ID = pmm_ID


        !
        ! Attach extended allocation to chidg_data%pmm
        !
        call move_alloc(temp_pmms,self%pmm)



    end function new_pmm
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

        integer(ik) :: idomain, nfields, eqn_ID

        ! Initialize mesh numerics based on equation set and polynomial expansion order
        call write_line(" ", ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line("Initialize: domain equation space...", io_proc=GLOBAL_MASTER)
        do idomain = 1,self%ndomains()
            eqn_ID = self%mesh(idomain)%eqn_ID
            nfields = self%eqnset(eqn_ID)%prop%nprimary_fields()
            call self%mesh(idomain)%init_sol(nfields,nterms_s,self%time_manager%ntime)
        end do

        call write_line(" ", ltrim=.false., io_proc=GLOBAL_MASTER)

    end subroutine initialize_solution_domains
    !******************************************************************************************






    
    !>  Initialize the 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine initialize_solution_bc(self)
        class(chidg_data_t),    intent(inout)   :: self

        integer(ik) :: ibc


        call write_line("BC: initializing communication...", io_proc=GLOBAL_MASTER)
        do ibc = 1,self%nbcs()

            !
            ! Prepare boundary condition parallel communication
            !
            call self%bc(ibc)%init_bc_comm(self%mesh)

            !
            ! Call bc-specific specialized routine. Default does nothing
            !
            call self%bc(ibc)%init_bc_specialized(self%mesh)

            !
            ! Initialize boundary condition coupling. 
            !
            call self%bc(ibc)%init_bc_coupling(self%mesh)

            !
            ! Propagate boundary condition coupling. 
            !
            call self%bc(ibc)%propagate_bc_coupling(self%mesh)

        end do



    end subroutine initialize_solution_bc
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

        do idom = 1,self%ndomains()
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
    function ndomains(self) result(ndomains_)
        class(chidg_data_t),    intent(in)      :: self

        integer :: ndomains_

        if (allocated(self%mesh)) then
            ndomains_ = size(self%mesh)
        else
            ndomains_ = 0
        end if

    end function ndomains
    !*******************************************************************************************





    !> Return the number of boundary conditions in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/10/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function nbcs(self) result(nbcs_)
        class(chidg_data_t),    intent(in)      :: self

        integer :: nbcs_

        if (allocated(self%bc)) then
            nbcs_ = size(self%bc)
        else
            nbcs_ = 0
        end if

    end function nbcs
    !*******************************************************************************************




    !> Return the number of equation sets in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function nequation_sets(self) result(nsets_)
        class(chidg_data_t),    intent(in)      :: self

        integer :: nsets_

        if (allocated(self%eqnset)) then
            nsets_ = size(self%eqnset)
        else
            nsets_ = 0
        end if

    end function nequation_sets
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

        integer(ik)                 :: idom, ifield, eqn_ID
        type(svector_t)             :: field_names
        character(:),   allocatable :: field_name



        do idom = 1,self%ndomains()
            eqn_ID = self%mesh(idom)%eqn_ID
            do ifield = 1,self%eqnset(eqn_ID)%prop%nauxiliary_fields()

                field_name = self%eqnset(eqn_ID)%prop%get_auxiliary_field_name(ifield)
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







    !>  Release allocated memory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_data_t),    intent(inout)   :: self

        if (allocated(self%mesh))   deallocate(self%mesh)
        if (allocated(self%eqnset)) deallocate(self%eqnset)
        if (allocated(self%bc))     deallocate(self%bc)

        call self%sdata%release()

    end subroutine release
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
