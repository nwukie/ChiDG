module type_chidg_data
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NO_ID
    use type_domain_connectivity,       only: domain_connectivity_t
    use type_boundary_connectivity,     only: boundary_connectivity_t

    ! Primary chidg_data_t components
    use type_point,                     only: point_t
    use type_mesh,                      only: mesh_t
    use type_bc_state,                  only: bc_state_t
    use type_bc_state_group,            only: bc_state_group_t
    use type_svector,                   only: svector_t
    use mod_string,                     only: string_t
    use type_equation_set,              only: equation_set_t
    use type_solverdata,                only: solverdata_t
    use type_time_manager,              only: time_manager_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_svector,                   only: svector_t
    use mod_string,                     only: string_t

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
        type(bc_state_group_t),         allocatable :: bc_state_group(:)
        type(equation_set_t),           allocatable :: eqnset(:)

        ! An object containing chidg matrices/vectors
        type(solverdata_t)                          :: sdata

        ! An object containing time information
        type(time_manager_t)                        :: time_manager


    contains

        ! Boundary conditions
        procedure   :: add_bc_state_group
        procedure   :: new_bc_state_group
        procedure   :: nbc_state_groups
        procedure   :: get_bc_state_group_id


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


        call write_line("Initialize: matrix/vector allocation...", io_proc=GLOBAL_MASTER)
        !
        ! Assemble array of function_data from the eqnset array to pass to the solver data 
        ! structure for initialization
        !
        ndom = self%mesh%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%mesh%ndomains()
            eqn_ID = self%mesh%domain(idom)%eqn_ID
            function_data(idom) = self%eqnset(eqn_ID)%function_data
        end do


        !
        ! Initialize solver data 
        !
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
    subroutine add_bc_state_group(self,bc_state_group,bc_wall,bc_inlet,bc_outlet,bc_symmetry,bc_farfield,bc_periodic, bc_scalar)
        class(chidg_data_t),            intent(inout)           :: self
        type(bc_state_group_t),         intent(inout)           :: bc_state_group
        class(bc_state_t),              intent(in), optional    :: bc_wall
        class(bc_state_t),              intent(in), optional    :: bc_inlet
        class(bc_state_t),              intent(in), optional    :: bc_outlet
        class(bc_state_t),              intent(in), optional    :: bc_symmetry
        class(bc_state_t),              intent(in), optional    :: bc_farfield
        class(bc_state_t),              intent(in), optional    :: bc_periodic
        class(bc_state_t),              intent(in), optional    :: bc_scalar


        integer(ik) :: bc_ID


        
        !
        ! Set override boundary condition states if they were passed in:
        !   if overriding:
        !       - clear state group
        !       - add overriding bc_state
        !
        if ( present(bc_wall) .and. (trim(bc_state_group%family) == 'Wall') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_wall)

        else if ( present(bc_inlet) .and. (trim(bc_state_group%family) == 'Inlet') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_inlet)

        else if ( present(bc_outlet) .and. (trim(bc_state_group%family) == 'Outlet') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_outlet)

        else if ( present(bc_symmetry) .and. (trim(bc_state_group%family) == 'Symmetry') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_symmetry)

        else if ( present(bc_farfield) .and. (trim(bc_state_group%family) == 'Farfield') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_farfield)

        else if ( present(bc_periodic) .and. (trim(bc_state_group%family) == 'Periodic') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_periodic)

        else if ( present(bc_scalar) .and. (trim(bc_state_group%family) == 'Scalar') ) then
            call bc_state_group%remove_states()
            call bc_state_group%add_bc_state(bc_scalar)

        end if



        !
        ! Assign:
        !   - Create new boundary condition state group
        !   - Set newly allocated object
        !
        bc_ID = self%new_bc_state_group()
        bc_state_group%bc_ID = bc_ID
        self%bc_state_group(bc_ID) = bc_state_group
                                                                            




    end subroutine add_bc_state_group
    !***************************************************************************************











    

    !>  Extend the self%bc array to include another instance. Return the ID of the new
    !!  boundary condition where it can be found in the array as self%bc(bc_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    function new_bc_state_group(self) result(bc_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(bc_state_group_t), allocatable     :: temp_bcs(:)
        integer(ik)                             :: bc_ID, ierr


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_bcs(self%nbc_state_groups() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( self%nbc_state_groups() > 0) then
            temp_bcs(1:size(self%bc_state_group)) = self%bc_state_group(1:size(self%bc_state_group))
        end if


        !
        ! Set ID of new bc and store to array
        !
        bc_ID = size(temp_bcs)
        temp_bcs(bc_ID)%bc_ID = bc_ID


        !
        ! Attach extended allocation to chidg_data%bc
        !
        call move_alloc(temp_bcs,self%bc_state_group)



    end function new_bc_state_group
    !*************************************************************************************







    !>  Return the number of boundary conditions state groups on the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/7/2017
    !!
    !!
    !-------------------------------------------------------------------------------------
    function nbc_state_groups(self) result(n)
        class(chidg_data_t),    intent(in)      :: self

        integer(ik) :: n

        if (allocated(self%bc_state_group)) then
            n = size(self%bc_state_group)
        else
            n = 0
        end if

    end function nbc_state_groups
    !*************************************************************************************







    !>  Given a group name for a bc_state_group, return its identifier so it can be 
    !!  located as self%bc_state_group(group_ID).
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/7/2017
    !!
    !-------------------------------------------------------------------------------------
    function get_bc_state_group_id(self,group_name_in) result(group_ID)
        class(chidg_data_t),    intent(in)  :: self
        character(*),           intent(in)  :: group_name_in

        integer(ik)                 :: igroup, group_ID
        character(:),   allocatable :: group_name
        logical                     :: found_group

        group_ID = NO_ID
        do igroup = 1,self%nbc_state_groups()

            found_group = trim(group_name_in) == trim(self%bc_state_group(igroup)%name)
            
            if (found_group) group_ID = igroup
            if (found_group) exit

        end do


    end function get_bc_state_group_id
    !*************************************************************************************








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
        integer(ik)                         :: eqn_ID, ierr


        !
        ! Allocate number of boundary conditions
        !
        allocate(temp_eqnset(self%nequation_sets() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any previously allocated boundary conditions to new array
        !
        if ( self%nequation_sets() > 0) then
            temp_eqnset(1:size(self%eqnset)) = self%eqnset(1:size(self%eqnset))
        end if


        !
        ! Set ID of new bc and store to array
        !
        eqn_ID = size(temp_eqnset)
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

            self%mesh%ntime_ = self%time_manager%ntime
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
        do ibc = 1,self%nbc_state_groups()

            !
            ! Prepare boundary condition parallel communication
            !
            call self%bc_state_group(ibc)%init_comm(self%mesh)

            !
            ! Call bc-specific specialized routine. Default does nothing
            !
            call self%bc_state_group(ibc)%init_specialized(self%mesh)

            !
            ! Initialize boundary condition coupling. 
            !
            call self%bc_state_group(ibc)%init_coupling(self%mesh)


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

        integer(ik)  :: domain_index
        
        domain_index = self%mesh%get_domain_id(domain_name)

    end function get_domain_index
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

        if (allocated(self%eqnset))         deallocate(self%eqnset)
        if (allocated(self%bc_state_group)) deallocate(self%bc_state_group)

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
