module type_chidg_data
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use type_point,                     only: point_t
    use type_domain_connectivity,       only: domain_connectivity_t
    use type_boundary_connectivity,     only: boundary_connectivity_t

    ! Primary chidg_data_t components
    use type_mesh,                      only: mesh_t
    use type_bc,                        only: bc_t
    use type_bc_state,                  only: bc_state_t
    use type_bc_state_group,            only: bc_state_group_t
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
    !--------------------------------------------------------------------------------------
    type, public  :: chidg_data_t

        logical                                     :: solverInitialized = .false.
        integer(ik),        private                 :: spacedim_ = 3

        ! mesh
        type(mesh_t)                                :: mesh

        ! bc's and equation set's
        type(bc_t),                     allocatable :: bc(:)    
        type(bc_state_group_t),         allocatable :: bc_state_group(:)
        type(equation_set_t),           allocatable :: eqnset(:)

        ! An object containing chidg matrices/vectors
        type(solverdata_t)                          :: sdata

        ! An object containing time information
        type(time_manager_t)                        :: time_manager


    contains

        ! Boundary conditions
        procedure   :: add_bc_group
        procedure   :: add_bc_patch
        procedure   :: new_bc
        procedure   :: nbcs


        ! Equations
        procedure   :: add_equation_set
        procedure   :: new_equation_set
        procedure   :: get_equation_set_id
        procedure   :: nequation_sets


        ! Initialization procedure for solution data. Execute after all domains are added.
        procedure   :: initialize_solution_domains
        procedure   :: initialize_solution_bc
        procedure   :: initialize_solution_solver


        ! Release allocated memory
        procedure   :: release


        ! Accessors
        procedure   :: get_domain_index             ! Given domain name, return domain index
        procedure   :: ntime
        procedure   :: get_dimensionality
        procedure   :: get_auxiliary_field_names    ! Return required auxiliary fields

        procedure   :: report

    end type chidg_data_t
    !***************************************************************************************




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
    !---------------------------------------------------------------------------------------
    subroutine initialize_solution_solver(self)
        class(chidg_data_t),     intent(inout)   :: self

        integer(ik) :: idom, ndom, ierr, eqn_ID

        type(equationset_function_data_t),  allocatable :: function_data(:)
        type(bcset_coupling_t),             allocatable :: bcset_coupling(:)


        call write_line("Initialize: matrix/vector allocation...", io_proc=GLOBAL_MASTER)
        !
        ! Assemble array of function_data from the eqnset array to pass to the solver data 
        ! structure for initialization
        !
        !ndom = self%ndomains()
        ndom = self%mesh%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        !do idom = 1,self%ndomains()
        do idom = 1,self%mesh%ndomains()
            eqn_ID = self%mesh%domain(idom)%eqn_ID
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
    !***************************************************************************************










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
    !!  @param[in]  bc_connectivity Face connectivities defining the boundary 
    !!                              condition patch.
    !!  @param[in]  bc_group        Name of the boundary condition group to associate 
    !!                              with the patch.
    !!  @param[in]  bc_groups       bc_group_t's for the global problem that can be searched 
    !!                              through and used to initialize.
    !!
    !!  To force a particular bc_state on a boundary condition, one can pass a bc_state_t in 
    !!  as an option for bc_wall, bc_inlet, bc_outlet, bc_symmetry
    !!
    !---------------------------------------------------------------------------------------
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
    !***************************************************************************************













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
    !---------------------------------------------------------------------------------------
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
            ! Find domain index in mesh from domain_name
            !
            idom = self%get_domain_index(domain_name)

            call self%bc(bc_ID)%init_bc_patch(self%mesh%domain(idom), bc_connectivity)


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
    !****************************************************************************************








    

    !>  Extend the self%bc array to include another instance. Return the ID of the new
    !!  boundary condition where it can be found in the array as self%bc(bc_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
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
    !***************************************************************************************










    !>  Add a new equation set to the data instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine add_equation_set(self, eqn_name)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: eqn_name

        logical     :: already_added
        integer(ik) :: ieqn, eqn_ID


        !
        ! Check if equation set has already been added.
        !
        already_added = .false.
        do ieqn = 1,self%nequation_sets()
            already_added = (trim(self%eqnset(ieqn)%name) == trim(eqn_name)) 
            if (already_added) exit
        end do
        

        !
        ! Add new equation set if it doesn't already exist and get new eqn_ID
        !
        if (.not. already_added) then
            eqn_ID = self%new_equation_set()
            self%eqnset(eqn_ID) = equation_builder_factory%produce(eqn_name, 'default')
            self%eqnset(eqn_ID)%eqn_ID = eqn_ID
        end if


    end subroutine add_equation_set
    !***************************************************************************************















    !>  Extend the self%eqnset array to include another instance. Return the ID of the new
    !!  equation set where it can be found in the array as self%eqnset(eqn_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
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
    !***************************************************************************************








    !>  Given an equation set name, return its index identifier in chidg_data.
    !!
    !!  Returns eqn_ID such that data%eqnset(eqn_ID) is valid and corresponds to the
    !!  equation set with name eqn_name that was given.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !-------------------------------------------------------------------------------------
    function get_equation_set_id(self, eqn_name) result(eqn_ID)
        class(chidg_data_t),    intent(in)  :: self
        character(*),           intent(in)  :: eqn_name

        integer(ik)                 :: ieqn, eqn_ID
        character(:),   allocatable :: user_msg
        logical                     :: names_match


        eqn_ID = 0
        do ieqn = 1,self%nequation_sets()

            names_match = trim(self%eqnset(ieqn)%name) == trim(eqn_name)

            if (names_match) eqn_ID = ieqn
            if (names_match) exit

        end do


        user_msg = "chidg_data%get_equation_set_id: No equation set was found that had a name &
                    matching the incoming string"
        if (eqn_ID == 0) call chidg_signal_one(FATAL,user_msg,eqn_name)


    end function get_equation_set_id
    !**************************************************************************************












    !>  For each domain, call solution initialization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine initialize_solution_domains(self,nterms_s)
        class(chidg_data_t),    intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms_s

        integer(ik) :: idomain, nfields, eqn_ID

        ! Initialize mesh numerics based on equation set and polynomial expansion order
        call write_line(" ", ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line("Initialize: domain equation space...", io_proc=GLOBAL_MASTER)

        do idomain = 1,self%mesh%ndomains()
            eqn_ID = self%mesh%domain(idomain)%eqn_ID
            nfields = self%eqnset(eqn_ID)%prop%nprimary_fields()
            call self%mesh%domain(idomain)%init_sol(nfields,nterms_s,self%time_manager%ntime)
        end do


    end subroutine initialize_solution_domains
    !***************************************************************************************






    
    !>  Initialize the 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine initialize_solution_bc(self)
        class(chidg_data_t),    intent(inout)   :: self

        integer(ik) :: ibc


        call write_line("Initialize: bc communication...", io_proc=GLOBAL_MASTER)
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
    !***************************************************************************************









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
    !----------------------------------------------------------------------------------------
    function get_domain_index(self,domain_name) result(domain_index)
        class(chidg_data_t),    intent(in)      :: self
        character(*),           intent(in)      :: domain_name

!        character(:),   allocatable :: user_msg
!        integer(ik)  :: idom
        integer(ik)  :: domain_index
        
!        domain_index = 0
!
!        do idom = 1,self%mesh%ndomains()
!            !if ( trim(domain_name) == trim(self%info(idom)%name) ) then
!            if ( trim(domain_name) == trim(self%mesh%domain(idom)%name) ) then
!                domain_index = idom
!                exit
!            end if
!        end do
!
!
!        user_msg = "chidg_data%get_domain_index: No domain was found that had a name &
!                    that matched the incoming string"
!        if (domain_index == 0) call chidg_signal_one(FATAL,user_msg,domain_name)
        domain_index = self%mesh%get_domain_id(domain_name)

    end function get_domain_index
    !****************************************************************************************







    !> Return the number of boundary conditions in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/10/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function nbcs(self) result(nbcs_)
        class(chidg_data_t),    intent(in)      :: self

        integer :: nbcs_

        if (allocated(self%bc)) then
            nbcs_ = size(self%bc)
        else
            nbcs_ = 0
        end if

    end function nbcs
    !****************************************************************************************







    !> Return the number of equation sets in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function nequation_sets(self) result(nsets_)
        class(chidg_data_t),    intent(in)      :: self

        integer :: nsets_

        if (allocated(self%eqnset)) then
            nsets_ = size(self%eqnset)
        else
            nsets_ = 0
        end if

    end function nequation_sets
    !****************************************************************************************





    !>  Return the dimensionality of the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_dimensionality(self) result(dimensionality)
        class(chidg_data_t),    intent(in)      :: self

        integer :: dimensionality

        dimensionality = self%spacedim_

    end function get_dimensionality
    !****************************************************************************************












    !>  Return a vector of auxiliary fields that are required.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !----------------------------------------------------------------------------------------
    function get_auxiliary_field_names(self) result(field_names)
        class(chidg_data_t),    intent(in)  :: self

        integer(ik)                 :: idom, ifield, eqn_ID
        type(svector_t)             :: field_names
        character(:),   allocatable :: field_name



        do idom = 1,self%mesh%ndomains()
            eqn_ID = self%mesh%domain(idom)%eqn_ID
            do ifield = 1,self%eqnset(eqn_ID)%prop%nauxiliary_fields()

                field_name = self%eqnset(eqn_ID)%prop%get_auxiliary_field_name(ifield)
                call field_names%push_back_unique(string_t(field_name))

            end do !ifield
        end do !idom



    end function get_auxiliary_field_names
    !****************************************************************************************





    !> Return the number of time instances in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function ntime(self) result(ndom)
        class(chidg_data_t),    intent(in)      :: self

        integer :: ndom

        ndom = self%time_manager%ntime

    end function ntime
    !****************************************************************************************







    !>  Release allocated memory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_data_t),    intent(inout)   :: self

        if (allocated(self%eqnset)) deallocate(self%eqnset)
        if (allocated(self%bc))     deallocate(self%bc)

        call self%mesh%release()
        call self%sdata%release()

    end subroutine release
    !****************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/7/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine report(self,selection)
        class(chidg_data_t),    intent(in)  :: self
        character(*),           intent(in)  :: selection

        integer(ik) :: idom


        if ( trim(selection) == 'grid' ) then

            do idom = 1,self%mesh%ndomains()
                call write_line('Domain ', idom, '  :  ', self%mesh%domain(idom)%nelem, ' Elements', io_proc=IRANK)
            end do


        else


        end if



    end subroutine report
    !****************************************************************************************






end module type_chidg_data
