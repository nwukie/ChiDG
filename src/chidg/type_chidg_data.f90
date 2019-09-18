module type_chidg_data
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NO_ID, NO_MM_ASSIGNED, ZERO, NO_ELEMENT, &
                                              MAX_ELEMENTS_PER_NODE, NFACES, TWO, THIRD
    use mod_chidg_mpi,                  only: IRANK, NRANK
    use type_domain_connectivity,       only: domain_connectivity_t
    use type_boundary_connectivity,     only: boundary_connectivity_t

    ! Primary chidg_data_t components
    use type_point,                     only: point_t
    use type_mesh,                      only: mesh_t
    use type_bc_state,                  only: bc_state_t
    use type_bc_state_group,            only: bc_state_group_t
    use type_svector,                   only: svector_t
    use type_ivector,                   only: ivector_t
    use mod_string,                     only: string_t
    use type_equation_set,              only: equation_set_t
    use type_solverdata,                only: solverdata_t
    use type_time_manager,              only: time_manager_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_svector,                   only: svector_t
    use mod_string,                     only: string_t

    ! Factory methods
    use mod_equations,                  only: equation_set_factory

    ! Vertex-based interpolation for h-smoothing
    use mod_interpolate,                only: interpolate_from_vertices

    !Mesh motion
    use type_mesh_motion,               only: mesh_motion_t
    use type_mesh_motion_wrapper,       only: mesh_motion_wrapper_t
    use type_mesh_motion_group,         only: mesh_motion_group_t
    use type_radial_basis_function
    use mod_radial_basis_function
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

        type(mesh_t)                                :: mesh
        type(solverdata_t)                          :: sdata
        type(time_manager_t)                        :: time_manager

        type(bc_state_group_t),         allocatable :: bc_state_group(:)
        type(equation_set_t),           allocatable :: eqnset(:)
        type(mesh_motion_wrapper_t),    allocatable :: mesh_motion(:)

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
        procedure   :: get_equation_set_name
        procedure   :: nequation_sets


        ! Initialization procedure for solution data. Execute after all domains are added.
        procedure   :: initialize_solution_domains
        procedure   :: initialize_solution_bc
        procedure   :: initialize_solution_solver
        procedure   :: initialize_postcomm_bc

        ! Mesh Motion
        procedure   :: new_mm
        procedure   :: add_mm_domain
        procedure   :: add_mm_group
        procedure   :: nmm_groups


        procedure   :: update_grid


        ! RBF
        procedure   :: construct_rbf_arrays
        procedure   :: node_index_u_to_g
        procedure   :: domain_index_global_to_local
        procedure   :: register_rbf_with_elements
        procedure   :: find_who_rbfs_touch

        ! Artificial Viscosity
        procedure   :: compute_area_weighted_h
        procedure   :: record_mesh_size
        procedure   :: perform_h_smoothing 

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
    !---------------------------------------------------------------------------------------
    subroutine initialize_solution_solver(self)
        class(chidg_data_t),     intent(inout)   :: self

        integer(ik) :: idom, ndom, ierr, eqn_ID

        type(equationset_function_data_t),  allocatable :: function_data(:)


        call write_line("Initialize: matrix/vector allocation...", io_proc=GLOBAL_MASTER)
        ! Assemble array of function_data from the eqnset array to pass to the solver data 
        ! structure for initialization
        ndom = self%mesh%ndomains()
        allocate(function_data(ndom), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        do idom = 1,self%mesh%ndomains()
            ! Assume that each element has the same eqn_ID
            eqn_ID = self%mesh%domain(idom)%elems(1)%eqn_ID
            function_data(idom) = self%eqnset(eqn_ID)%function_data
        end do


        ! Initialize solver data 
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


        ! Set override boundary condition states if they were passed in:
        !   if overriding:
        !       - clear state group
        !       - add overriding bc_state
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


        ! Assign:
        !   - Create new boundary condition state group
        !   - Set newly allocated object
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
    !---------------------------------------------------------------------------------------
    function new_bc_state_group(self) result(bc_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(bc_state_group_t), allocatable     :: temp_bcs(:)
        integer(ik)                             :: bc_ID, ierr


        ! Allocate number of boundary conditions
        allocate(temp_bcs(self%nbc_state_groups() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Copy any previously allocated boundary conditions to new array
        if ( self%nbc_state_groups() > 0) then
            temp_bcs(1:size(self%bc_state_group)) = self%bc_state_group(1:size(self%bc_state_group))
        end if


        ! Set ID of new bc and store to array
        bc_ID = size(temp_bcs)
        temp_bcs(bc_ID)%bc_ID = bc_ID


        ! Attach extended allocation to chidg_data%bc
        call move_alloc(temp_bcs,self%bc_state_group)


    end function new_bc_state_group
    !*************************************************************************************





    !>  Return the number of boundary conditions state groups on the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/7/2017
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
        character(:),   allocatable :: group_name, user_msg
        logical                     :: found_group

        group_ID = NO_ID
        do igroup = 1,self%nbc_state_groups()

            found_group = trim(group_name_in) == trim(self%bc_state_group(igroup)%name)
            
            if (found_group) group_ID = igroup
            if (found_group) exit

        end do

        if ((group_ID == NO_ID) .and. &
            (trim(group_name_in) /= 'empty') .and. &
            (trim(group_name_in) /= 'Empty') .and. &
            (trim(group_name_in) /= 'EMPTY')) then
            user_msg = "chidg_data%get_bc_state_group_id: Didn't find boundary condition state group. Make sure this group exists and that the patch group assignment is spelled correctly."
            call chidg_signal_one(FATAL,user_msg,trim(group_name_in))
        end if

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

        ! Check if equation set has already been added.
        already_added = .false.
        do ieqn = 1,self%nequation_sets()
            already_added = (trim(self%eqnset(ieqn)%name) == trim(eqn_name)) 
            if (already_added) exit
        end do
        
        ! Add new equation set if it doesn't already exist and get new eqn_ID
        if (.not. already_added) then
            eqn_ID = self%new_equation_set()
            self%eqnset(eqn_ID) = equation_set_factory%produce(eqn_name, 'default')
            self%eqnset(eqn_ID)%eqn_ID = eqn_ID
        end if


    end subroutine add_equation_set
    !***************************************************************************************






    !
    !   Mesh Motion
    !------------------------------------------------------------------------------------------


    !> 
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_mm_group(self,mm_group)
        class(chidg_data_t),            intent(inout)           :: self
        type(mesh_motion_group_t),               intent(inout)              :: mm_group

        integer(ik)                         :: mm_ID

        ! Create a new boundary condition
        mm_ID = self%new_mm()

        ! Initialize mm from mm_group. Note that, since mm_group contains mm,
        ! we must pass the mm info instead of mm_group to avoid a circular dependency.
        allocate(self%mesh_motion(mm_ID)%mm, source=mm_group%mm)
        self%mesh_motion(mm_ID)%mm%mm_ID = mm_ID

    end subroutine add_mm_group
    !******************************************************************************************




    !>  Return number of prescribed mesh motion groups.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/14/2018
    !!
    !------------------------------------------------------------------------------------------
    function nmm_groups(self) result(nmm)
        class(chidg_data_t),    intent(in)  :: self

        integer(ik) :: nmm

        if (allocated(self%mesh_motion)) then
            nmm = size(self%mesh_motion)
        else
            nmm = 0
        end if

    end function nmm_groups
    !******************************************************************************************




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


        ! Allocate number of boundary conditions
        allocate(temp_eqnset(self%nequation_sets() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Copy any previously allocated boundary conditions to new array
        if ( self%nequation_sets() > 0) then
            temp_eqnset(1:size(self%eqnset)) = self%eqnset(1:size(self%eqnset))
        end if


        ! Set ID of new bc and store to array
        eqn_ID = size(temp_eqnset)
        temp_eqnset(eqn_ID)%eqn_ID = eqn_ID


        ! Attach extended allocation to chidg_data%eqnset
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




    !> 
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_mm_domain(self, domain_name, domain_mm_name)
        class(chidg_data_t),            intent(inout)   :: self
        character(*),                   intent(in)      :: domain_name
        character(*),                   intent(in)      :: domain_mm_name

        character(:),   allocatable :: mm_name, user_msg
        integer(ik)                 :: mm_ID, imm, idom
        logical                     :: found_mm

        ! Find the correct boundary condition to add bc_patch to
        do imm = 1,size(self%mesh_motion)

            mm_name = self%mesh_motion(imm)%mm%get_name()
            found_mm = (trim(domain_mm_name) == trim(mm_name))

            if (found_mm) mm_ID = imm
            if (found_mm) exit

        end do

        ! Once mm is found, initialize mm_patch on mm
        if (found_mm) then
            ! Find domain index in mesh(:) from domain_name
            idom = self%get_domain_index(domain_name)
            call self%mesh_motion(mm_ID)%mm%init_mm_domain(self%mesh%domain(idom))
        else

            user_msg = "chidg_data%add_mm_domain: It looks like we didn't find a mesh motion  &
                        group that matches with the string indicated in a mm domain. Make &
                        sure that a mm group with the correct name exists. Also make &
                        sure that the name set on the mm domain corresponds to one of the &
                        mm groups that exists."
            if ( (trim(domain_mm_name) /= 'empty') .and. &
                (trim(domain_mm_name) /= 'Empty') ) &
                call chidg_signal_one(FATAL,user_msg,trim(domain_mm_name))

        end if


    end subroutine add_mm_domain
    !*******************************************************************************************



    !>
    !!  @author Eric Wolf
    !!  @date   4/3/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_mm(self) result(mm_ID)
        class(chidg_data_t),    intent(inout)   :: self

        type(mesh_motion_wrapper_t), allocatable :: temp_mms(:)
        integer(ik)             :: mm_ID, ierr, nmm

        ! Get number of boundary conditions
        if (allocated(self%mesh_motion)) then
            nmm = size(self%mesh_motion)
        else
            nmm = 0
        end if

        ! Increment number of boundary conditions
        nmm = nmm + 1

        ! Allocate number of boundary conditions
        allocate(temp_mms(nmm), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Copy any previously allocated boundary conditions to new array
        if ( nmm > 1) then
            temp_mms(1:size(self%mesh_motion)) = self%mesh_motion(1:size(self%mesh_motion))
        end if


        ! Set ID of new mm and store to array
        mm_ID = nmm

        ! Attach extended allocation to chidg_data%mm
        call move_alloc(temp_mms,self%mesh_motion)

    end function new_mm
    !******************************************************************************************




    !>  Spatial loop through domains, elements, and faces. Functions get called for each element/face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Improved layout, added computation of diffusion terms.
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine update_grid(self,timing,info)
        class(chidg_data_t), intent(inout)   :: self 
        real(rk),           optional        :: timing
        integer(ik),        optional        :: info

        integer(ik) :: inode, mm_ID, idom, ierr, imm


        call write_line('Updating to mesh motion...', io_proc=GLOBAL_MASTER)

        ! Update mesh motion objects
        !do imm = 1, size(self%mesh_motion) 
        do imm = 1, self%nmm_groups()
            call self%mesh_motion(imm)%mm%update(self%mesh, self%time_manager%t)
            call self%mesh_motion(imm)%mm%apply(self%mesh, self%time_manager%t)
        end do

        ! Loop through domains and apply mesh motions
        do idom = 1,self%mesh%ndomains()
            mm_ID = self%mesh%domain(idom)%mm_ID
            if (mm_ID /= NO_MM_ASSIGNED) then
                call self%mesh%domain(idom)%set_displacements_velocities(self%mesh%domain(idom)%dnodes, self%mesh%domain(idom)%vnodes)
                call self%mesh%domain(idom)%update_interpolations_ale()
            end if
        end do  ! idom


        call self%mesh%comm_send()
        call self%mesh%comm_recv()
        call self%mesh%comm_wait()


    end subroutine update_grid
    !******************************************************************************************************************

    

   

    !>  Given an equation set name, return its index identifier in chidg_data.
    !!
    !!  Returns eqn_ID such that data%eqnset(eqn_ID) is valid and corresponds to the
    !!  equation set with name eqn_name that was given.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !-------------------------------------------------------------------------------------
    function get_equation_set_name(self, eqn_ID) result(eqn_name)
        class(chidg_data_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: eqn_ID

        character(:),   allocatable :: eqn_name, user_msg

        ! Check if eqn_ID is within bounds
        if ( eqn_ID > self%nequation_sets() ) call chidg_signal(FATAL,"chidg_data%get_equation_set_name: eqn_ID is out of bounds.")
        if ( eqn_ID < 1 )                     call chidg_signal(FATAL,"chidg_data%get_equation_set_name: eqn_ID is out of bounds.")

        ! Get name
        eqn_name = self%eqnset(eqn_ID)%name

    end function get_equation_set_name
    !**************************************************************************************





    !>  For each domain, call solution initialization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine initialize_solution_domains(self,interpolation,level,nterms_s)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: interpolation
        integer(ik),            intent(in)      :: level
        integer(ik),            intent(in)      :: nterms_s

        integer(ik) :: idomain, nfields, ntime, eqn_ID, domain_dof_start, domain_dof_local_start

        ! Initialize mesh numerics based on equation set and polynomial expansion order
        call write_line(" ", ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line("Initialize: domain equation space...", io_proc=GLOBAL_MASTER)

        ! Set interpolation node set
        self%mesh%interpolation = interpolation

        ! Initialize domains
        do idomain = 1,self%mesh%ndomains()

            ! Initialize mesh_dof_start
            eqn_ID = self%mesh%domain(idomain)%elems(1)%eqn_ID
            ntime   = self%time_manager%ntime
            nfields = self%eqnset(eqn_ID)%prop%nprimary_fields()
            self%mesh%mesh_dof_start = sum(self%mesh%nelements_per_proc(1:IRANK))*nterms_s*nfields*ntime + 1


            ! Get the starting dof index for the domain
            if (idomain==1) then
                domain_dof_start       = self%mesh%mesh_dof_start
                domain_dof_local_start = 1
            else
                domain_dof_start       = self%mesh%domain(idomain-1)%get_dof_end() + 1
                domain_dof_local_start = self%mesh%domain(idomain-1)%get_dof_local_end() + 1
            end if

            ! Call initialization
            self%mesh%ntime_ = self%time_manager%ntime
            call self%mesh%domain(idomain)%init_sol(interpolation,level,nterms_s,nfields,self%time_manager%ntime,domain_dof_start,domain_dof_local_start)
            call self%mesh%domain(idomain)%update_interpolations_ale()
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

        integer(ik) :: ibc, ierr

        call write_line("Initialize: bc communication...", io_proc=GLOBAL_MASTER)
        do ibc = 1,self%nbc_state_groups()

            ! Prepare boundary condition parallel communication
            call self%bc_state_group(ibc)%init_comm(self%mesh)

            ! Call bc-specific specialized routine. Default does nothing
            call self%bc_state_group(ibc)%init_precomm(self%mesh)

            ! Initialize boundary condition coupling. 
            call self%bc_state_group(ibc)%init_coupling(self%mesh)

        end do

    end subroutine initialize_solution_bc
    !***************************************************************************************





    !>  Initialize the 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine initialize_postcomm_bc(self)
        class(chidg_data_t),    intent(inout)   :: self

        integer(ik) :: ibc

        call write_line("Initialize: bc specializations...", io_proc=GLOBAL_MASTER)
        do ibc = 1,self%nbc_state_groups()
            ! Call bc-specific specialized routine. Default does nothing
            call self%bc_state_group(ibc)%init_postcomm(self%mesh)
        end do

    end subroutine initialize_postcomm_bc
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
            ! Assume each element has the same eqn_ID
            eqn_ID = self%mesh%domain(idom)%elems(1)%eqn_ID
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




    !>
    !! 
    !! @author  Eric M. Wolf
    !! @date    03/05/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_area_weighted_h(self)
        class(chidg_data_t),    intent(inout) :: self

        integer(ik) :: idom

        do idom = 1, self%mesh%ndomains()
           call self%mesh%domain(idom)%compute_area_weighted_h() 
        end do

    end subroutine compute_area_weighted_h 
    !********************************************************************************










    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    07/16/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine construct_rbf_arrays(self)
        class(chidg_data_t),    intent(inout)      :: self

        integer(ik)     :: idom, ielem, rbf_index, idir

        self%sdata%rbf_center = ZERO
        self%sdata%rbf_radius = ZERO
        do idom = 1, self%mesh%ndomains()
            do ielem = 1, self%mesh%domain(idom)%nelem
                
                rbf_index = sum(self%sdata%nelems_per_domain(1:self%mesh%domain(idom)%elems(ielem)%idomain_g-1)) + &
                            self%mesh%domain(idom)%elems(ielem)%ielement_g

                do idir = 1,3
                    self%sdata%rbf_center(rbf_index,idir) = self%mesh%domain(idom)%elems(ielem)%bb_centroid(idir)
                    self%sdata%rbf_radius(rbf_index,idir) = self%mesh%domain(idom)%elems(ielem)%h(idir)
                end do

            end do
        end do

    end subroutine construct_rbf_arrays
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/13/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine node_index_u_to_g(self, inode_u, idomain_g, inode_g) 
        class(chidg_data_t),    intent(in)      :: self
        integer(ik),            intent(in)      :: inode_u
        integer(ik),            intent(inout)   :: idomain_g
        integer(ik),            intent(inout)   :: inode_g

        integer(ik) :: nnodes_old, nnodes, idom, ndomains_g

        ndomains_g = size(self%sdata%nnodes_per_domain)

        nnodes_old = 0
        nnodes = 0
        do idom = 1,ndomains_g
            nnodes = nnodes_old + self%sdata%nnodes_per_domain(idom)

            if ((nnodes_old<inode_u) .and. (inode_u<=nnodes)) then
                ! Node located in the present domain
                idomain_g = idom
                inode_g = inode_u - nnodes_old
                exit

            end if
            nnodes_old = nnodes

        end do

    end subroutine node_index_u_to_g
    !****************************************************************************************




    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    09/13/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine domain_index_global_to_local(self,idomain_g,idomain_l,domain_found) 
        class(chidg_data_t),    intent(in)      :: self
        integer(ik),            intent(in)      :: idomain_g
        integer(ik),            intent(inout)   :: idomain_l
        logical,                intent(inout)   :: domain_found

        integer(ik) :: idom, ndomains

        domain_found = .false.
        ndomains = self%mesh%ndomains()
        idomain_l = -1

        do idom = 1,ndomains
            if (idomain_g == self%mesh%domain(idom)%idomain_g) then
                idomain_l = idom
                domain_found = .true.
                exit
            end if
        end do

    end subroutine domain_index_global_to_local
    !********************************************************************************





    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    03/20/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine find_who_rbfs_touch(self)
        class(chidg_data_t),    intent(inout) :: self

        integer(ik) :: irbf, nrbf
        real(rk)    :: center(3), radius(3), alpha(3)

        ! Loop over the element-centered RBFs and register them with each element
        alpha = TWO
        nrbf = size(self%sdata%rbf_center(:,1))
        do irbf = 1,nrbf
            center = self%sdata%rbf_center(irbf,:)
            radius = alpha*self%sdata%rbf_radius(irbf,:)
            call self%register_rbf_with_elements(center, radius, irbf, 1)
        end do

    end subroutine find_who_rbfs_touch
    !********************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/13/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine register_rbf_with_elements(self, center, radius, rbf_ID, rbf_set_ID)
        class(chidg_data_t),                intent(inout) :: self
        real(rk),                           intent(in)      :: center(3), radius(3)
        integer(ik),                        intent(in)      :: rbf_ID, rbf_set_ID

        integer(ik)     :: inode, inode_u, idomain_g, inode_g, idomain_l, ielem, nelems, test_elem
        logical         :: domain_found
        type(ivector_t) :: hit_list

        ! Perform radius neighbor search
        call self%mesh%octree%radius_search(self%sdata%global_nodes,self%mesh%octree%root_box_ID,center, radius, hit_list)

        do inode = 1, hit_list%size()
            inode_u = hit_list%at(inode)
            ! The inode_u index corresponds to the node's index in the global_nodes array,
            ! which contains all nodes in the entire mesh.
            ! The subroutine node_index_u_to_g gets the corresponding idomain_g and inode_g index.
            call self%node_index_u_to_g(inode_u, idomain_g, inode_g) 

            ! The global domain idomain_g may or may not exist on the local processor.
            ! If it does, we get the idomain_l index.
            call self%domain_index_global_to_local(idomain_g, idomain_l, domain_found)

            if (domain_found) then
                nelems = 0
                ! Loop over all of the elements containing this node (up to MAX_ELEMENTS_PER_NODE=8 elements).
                ! There may be none, if this node is off-processor.
                do ielem = 1, MAX_ELEMENTS_PER_NODE
                    test_elem = self%mesh%domain(idomain_l)%nodes_elems(inode_g, ielem)

                    ! If inode_g is off-processor, then nodes_elems array should contain only NO_ELEMENT
                    ! Otherwise, it will contain at least one valid ielement_l index.
                    if (test_elem /= NO_ELEMENT) then
                        call self%mesh%domain(idomain_l)%elems(test_elem)%register_rbf(rbf_set_ID, rbf_ID)
                    end if

                end do
            end if

        end do

    end subroutine register_rbf_with_elements
    !********************************************************************************




    !>  Record mesh size
    !!
    !!  @author Eric M. Wolf 
    !!  @date   09/18/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine record_mesh_size(self)
        class(chidg_data_t),    intent(inout)   :: self

        integer(ik) :: idomain, ielem, idomain_g, idomain_l, ielement_g, &
                       ielement_l, ii, inode_loc, inode_global, idir
        real(rk)    :: h_loc(3)

        self%sdata%mesh_size_elem       = ZERO
        self%sdata%mesh_size_vertex     = ZERO
        self%sdata%min_mesh_size_vertex = ZERO
        self%sdata%avg_mesh_size_vertex = ZERO
        self%sdata%sum_mesh_size_vertex = ZERO
        self%sdata%num_elements_touching_vertex = 0

        do idomain = 1,self%mesh%ndomains()
            do ielem = 1, self%mesh%domain(idomain)%nelem

                idomain_g  = self%mesh%domain(idomain)%elems(ielem)%idomain_g
                idomain_l  = self%mesh%domain(idomain)%elems(ielem)%idomain_l
                ielement_g = self%mesh%domain(idomain)%elems(ielem)%ielement_g
                ielement_l = self%mesh%domain(idomain)%elems(ielem)%ielement_l
                ii = sum(self%sdata%nelems_per_domain(1:idomain_g-1))+ielement_g
                
                ! OPTION 1: Get element h from bounding box
                h_loc = self%mesh%domain(idomain)%elems(ielem)%h

                ! OPTION 2: Get element h from area weighted h
                !h_loc = self%mesh%domain(idomain)%elems(ielem)%area_weighted_h

                ! Record the present elements mesh size
                self%sdata%mesh_size_elem(ii, :)  = h_loc
                self%sdata%min_mesh_size_elem(ii) = minval(h_loc)
                self%sdata%area_weighted_h(ii,:)  = self%mesh%domain(idomain)%elems(ielem)%area_weighted_h

                ! Loop over the nodes belonging to the element and overwrite if the new value is larger than the previous value
                do inode_loc = 1, size(self%mesh%domain(idomain)%elems(ielem)%connectivity)

                    ! Get the nodes index
                    inode_global = sum(self%sdata%nnodes_per_domain(1:idomain_g-1))+self%mesh%domain(idomain)%elems(ielem)%connectivity(inode_loc)
                    if (inode_global > sum(self%sdata%nnodes_per_domain)) print *, 'av node out of bounds'

                    ! Compare with the current nodal value, overwrite if larger
                    do idir = 1, 3
                        if (self%sdata%mesh_size_vertex(inode_global, idir) < h_loc(idir)) self%sdata%mesh_size_vertex(inode_global, idir) = h_loc(idir) 
                        if (self%sdata%min_mesh_size_vertex(inode_global) < h_loc(idir)) self%sdata%min_mesh_size_vertex(inode_global) = h_loc(idir) 
                    end do

                    ! Record nodal information
                    ! NOTE: this statement fills only a small portion of the allocated arrays since
                    ! they are representative of the global problem. Later, via MPI Reductions, their 
                    ! data will be synchronized across all participating ranks (See type_chidg%read_mesh).
                    self%sdata%sum_mesh_size_vertex(inode_global)         = self%sdata%sum_mesh_size_vertex(inode_global) + (h_loc(1)*h_loc(2)*h_loc(3))**(THIRD)
                    self%sdata%sum_mesh_h_vertex(inode_global,:)          = self%sdata%sum_mesh_h_vertex(inode_global,:) + h_loc 
                    self%sdata%num_elements_touching_vertex(inode_global) = self%sdata%num_elements_touching_vertex(inode_global) + 1

                end do !inode_loc

            end do !ielem
        end do !idomain

    end subroutine record_mesh_size
    !***************************************************************************************





    !> This subroutine creates a continous h-field, used to generate a smooth(er) artificial viscosity.
    !!
    !! @author  Eric M. Wolf
    !! @date    03/01/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine perform_h_smoothing(self)
        class(chidg_data_t),        intent(inout) :: self

        real(rk)        :: h_loc(3), vert_vals_h(8,3), vert_vals_size(8)
        integer(ik)     :: idomain, ielement, ii, ivertex, inode, inode_g, iface, idomain_g, idomain_l, ielement_g, ielement_l, &
                           nnodes_elem, nnodes_face, idir, ierr

        real(rk),   allocatable :: h_field_elem(:,:), h_field_face(:,:), size_field_elem(:), size_field_face(:), nodes_elem(:,:), nodes_face(:,:)
        real(rk),   allocatable :: normalizer_elem(:), normalizer_face(:), node1(:), node2(:), node3(:)
        real(rk)                :: alpha(3), h_field_coeff(3), size_field_coeff, center(3), radius(3), evalnode(3)

        integer(ik), allocatable    :: rbf_index_list(:)
        integer(ik)                 :: irbf, nrbf
        class(radial_basis_function_t), allocatable  :: rbf

        integer                     :: unit, msg
        logical                     :: file_exists
        character(len=100)          :: rbf_type, smoothing_method

        namelist /h_smoothing/ smoothing_method 
        namelist /rbf_av/      rbf_type

        alpha = TWO

        ! Look in models.nml to get alpha (RBF radius expansion factors)
        inquire(file='artificial_viscosity.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='artificial_viscosity.nml')
            read(unit,nml=rbf_av,iostat=msg)
            if (msg /= 0) print *, 'rbf_av NOT read, msg = ', msg
            if (msg /= 0) rbf_type = 'wc2_rbf'
            close(unit)
        else
            rbf_type     = 'wc2_rbf'
        end if       


        if (file_exists) then
            open(newunit=unit,form='formatted',file='artificial_viscosity.nml')
            read(unit,nml=h_smoothing,iostat=msg)
            if (msg /= 0) print *, 'h_smoothing NOT read, msg = ', msg
            if (msg /= 0) smoothing_method = 'rbf'
            close(unit)
        else
            smoothing_method = 'rbf'
        end if       



        select case (trim(smoothing_method))
            ! Vertex-based interpolation
            case ('vbi')
                do idomain = 1,self%mesh%ndomains()
                    do ielement = 1, self%mesh%domain(idomain)%nelements()

                        idomain_g  = self%mesh%domain(idomain)%elems(ielement)%idomain_g
                        idomain_l  = self%mesh%domain(idomain)%elems(ielement)%idomain_l
                        ielement_g = self%mesh%domain(idomain)%elems(ielement)%ielement_g
                        ielement_l = self%mesh%domain(idomain)%elems(ielement)%ielement_l

                        ! Collect vertex h values
                        do ivertex = 1, 8
                            inode_g = sum(self%sdata%nnodes_per_domain(1:idomain_g-1)) + self%mesh%domain(idomain)%elems(ielement)%vertex_indices(ivertex)

                            vert_vals_size(ivertex) = self%sdata%avg_mesh_size_vertex(inode_g)
                            do idir = 1, 3
                                vert_vals_h(ivertex,idir) = self%sdata%avg_mesh_h_vertex(inode_g,idir)
                            end do
                        end do

                        ! Interpolate to volume quadrature nodes
                        nodes_elem = self%mesh%domain(idomain)%elems(ielement)%basis_s%nodes_elem_
                        nnodes_elem = size(nodes_elem(:,1))

                        if (allocated(size_field_elem)) deallocate(size_field_elem)
                        if (allocated(h_field_elem))    deallocate(h_field_elem)
                        allocate(size_field_elem(nnodes_elem), &
                                 h_field_elem(nnodes_elem,3), stat=ierr)
                        if (ierr /= 0) call AllocationError

                        do inode = 1,nnodes_elem 
                            size_field_elem(inode) = interpolate_from_vertices(vert_vals_size, nodes_elem(inode,:))
                            do idir=1,3
                                h_field_elem(inode, idir) = interpolate_from_vertices(vert_vals_h(:,idir), nodes_elem(inode,:))
                            end do
                        end do
                        self%mesh%domain(idomain)%elems(ielement)%h_smooth    = h_field_elem
                        self%mesh%domain(idomain)%elems(ielement)%size_smooth = size_field_elem

                        ! Interpolate to face quadrature nodes
                        do iface = 1, NFACES
                            nodes_face = self%mesh%domain(idomain)%elems(ielement)%basis_s%nodes_face_(:,:, iface)
                            nnodes_face = size(nodes_face(:,1))
                            if (allocated(size_field_face)) deallocate(size_field_face)
                            if (allocated(h_field_face))    deallocate(h_field_face)
                            allocate(size_field_face(nnodes_face), &
                                     h_field_face(nnodes_face, 3), stat=ierr)
                            if (ierr /= 0) call AllocationError
                            do inode = 1, nnodes_face 
                                size_field_face(inode) = interpolate_from_vertices(vert_vals_size, nodes_face(inode,:))
                                do idir = 1,3
                                    h_field_face(inode, idir) = interpolate_from_vertices(vert_vals_h(:,idir), nodes_face(inode,:))
                                end do
                            end do
                            self%mesh%domain(idomain)%faces(ielement, iface)%h_smooth = h_field_face
                            self%mesh%domain(idomain)%faces(ielement, iface)%size_smooth = size_field_face

                        end do !iface

                    end do !ielement
                end do !idomain

            ! Smoothing via Radial Basis Function Interpolation
            case ('rbf')

                call create_radial_basis_function(rbf,trim(rbf_type))
                do idomain = 1,self%mesh%ndomains()
                    do ielement = 1, self%mesh%domain(idomain)%nelements()

                        idomain_g  = self%mesh%domain(idomain)%elems(ielement)%idomain_g
                        idomain_l  = self%mesh%domain(idomain)%elems(ielement)%idomain_l
                        ielement_g = self%mesh%domain(idomain)%elems(ielement)%ielement_g
                        ielement_l = self%mesh%domain(idomain)%elems(ielement)%ielement_l

                        if (allocated(rbf_index_list)) deallocate(rbf_index_list)
                        rbf_index_list = self%mesh%domain(idomain)%elems(ielement)%get_rbf_indices(1)
                        nrbf = size(rbf_index_list)
                        
                        if (nrbf==0) print *, 'no rbfs found for h smoothing'

                        node1 = self%mesh%domain(idomain_l)%elems(ielement_l)%interp_coords_def(:,1)
                        node2 = self%mesh%domain(idomain_l)%elems(ielement_l)%interp_coords_def(:,2)
                        node3 = self%mesh%domain(idomain_l)%elems(ielement_l)%interp_coords_def(:,3)

                        nodes_elem = self%mesh%domain(idomain)%elems(ielement)%basis_s%nodes_elem_
                        nnodes_elem = size(nodes_elem(:,1))

                        if (allocated(size_field_elem)) deallocate(size_field_elem)
                        if (allocated(h_field_elem))    deallocate(h_field_elem)
                        if (allocated(normalizer_elem)) deallocate(normalizer_elem)

                        allocate(size_field_elem(nnodes_elem), &
                                 h_field_elem(nnodes_elem,3),  &
                                 normalizer_elem(nnodes_elem), stat=ierr)
                        if (ierr /= 0) call AllocationError
                        
                        size_field_elem = ZERO
                        h_field_elem    = ZERO
                        normalizer_elem = ZERO


                        do irbf = 1, nrbf

                            ii     = rbf_index_list(irbf)
                            center = self%sdata%rbf_center(ii,:)
                            radius = alpha*self%sdata%rbf_radius(ii,:)

                            h_field_coeff       = self%sdata%rbf_radius(ii,:)
                            !h_field_coeff       = self%sdata%area_weighted_h(ii,:)
                            size_field_coeff    = sum(h_field_coeff)/3.0_rk
        
                            do inode = 1,nnodes_elem 
                                evalnode = [node1(inode), node2(inode), node3(inode)]
                                                           
                                normalizer_elem(inode)  = normalizer_elem(inode)    + rbf%compute(evalnode, center, radius)
                                size_field_elem(inode)  = size_field_elem(inode)    + size_field_coeff*rbf%compute(evalnode, center, radius)        
                                h_field_elem(inode,:)   = h_field_elem(inode,:)     + h_field_coeff*rbf%compute(evalnode, center, radius)        
                            end do
                        end do

                        do idir = 1,3
                            self%mesh%domain(idomain)%elems(ielement)%h_smooth(:,idir) = h_field_elem(:,idir)/normalizer_elem
                        end do
                        self%mesh%domain(idomain)%elems(ielement)%size_smooth = size_field_elem/normalizer_elem

                        ! Interpolate to face quadrature nodes
                        do iface = 1, NFACES

                            nodes_face = self%mesh%domain(idomain)%elems(ielement)%basis_s%nodes_face_(:,:, iface)
                            nnodes_face = size(nodes_face(:,1))
                            if (allocated(size_field_face)) deallocate(size_field_face)
                            if (allocated(h_field_face))    deallocate(h_field_face)
                            if (allocated(normalizer_face)) deallocate(normalizer_face)
                            allocate(size_field_face(nnodes_face), &
                                     h_field_face(nnodes_face, 3), &
                                     normalizer_face(nnodes_face), stat=ierr)
                            if (ierr /= 0) call AllocationError

                            node1 = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,1)
                            node2 = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,2)
                            node3 = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,3)

                            size_field_face = ZERO
                            h_field_face    = ZERO
                            normalizer_face = ZERO

                            do irbf = 1,nrbf
                                ii     = rbf_index_list(irbf)
                                center = self%sdata%rbf_center(ii,:)
                                radius = alpha*self%sdata%rbf_radius(ii,:)

                                h_field_coeff    = self%sdata%rbf_radius(ii,:)
                                size_field_coeff = sum(h_field_coeff)/3.0_rk

                                do inode = 1,nnodes_face
                                    evalnode = [node1(inode), node2(inode), node3(inode)]
                             
                                    normalizer_face(inode) = normalizer_face(inode) + rbf%compute(evalnode, center, radius)
                                    size_field_face(inode) = size_field_face(inode) + size_field_coeff*rbf%compute(evalnode, center, radius)        
                                    h_field_face(inode,:)  = h_field_face(inode,:)  + h_field_coeff*rbf%compute(evalnode, center, radius)        
                                end do
                            end do

                            do idir = 1,3
                                self%mesh%domain(idomain)%faces(ielement,iface)%h_smooth(:,idir) = h_field_face(:,idir)/normalizer_face
                            end do
                            self%mesh%domain(idomain)%faces(ielement, iface)%size_smooth = size_field_face/normalizer_face
                        end do ! iface

                    end do ! ielem
                end do ! idomain

            case default
                call chidg_signal(FATAL,'chidg_data%perform_h_smoothing: invalid parameter for h-smoothing interpolation.')
        end select

    end subroutine perform_h_smoothing 
    !********************************************************************************
















    !>  Release allocated memory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
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
