module type_solverdata
#include <messenger.h>
    use mod_kinds,                          only: rk,ik
    use mod_constants,                      only: NFACES, ZERO, NO_ID
    use mod_string,                         only: string_t
    use mod_io,                             only: backend
    use type_chidg_vector,                  only: chidg_vector_t, chidg_vector
    use type_chidg_matrix,                  only: chidg_matrix_t, chidg_matrix
    use type_mesh,                          only: mesh_t
    use type_function_status,               only: function_status_t
    use type_equationset_function_data,     only: equationset_function_data_t
    use type_element_info,                  only: element_info_t
    use type_face_info,                     only: face_info_t
    use type_function_info,                 only: function_info_t
    implicit none



    !> Container for solver data.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   3/15/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public  :: solverdata_t

        ! Base solver data
        type(chidg_vector_t)            :: q           ! Solution vector
        type(chidg_vector_t)            :: dq          ! Change in solution vector
        type(chidg_vector_t)            :: rhs         ! Residual of the spatial scheme
        type(chidg_matrix_t)            :: lhs         ! Linearization of the spatial scheme

 
        ! Container for reading data
        type(chidg_vector_t)            :: q_in     ! For reading data
        type(chidg_vector_t)            :: q_out    ! For post-processing

        
        ! Auxiliary fields
        type(string_t),         allocatable :: auxiliary_field_name(:)
        type(chidg_vector_t),   allocatable :: auxiliary_field(:)

        ! Time information
        real(rk),       allocatable :: dt(:,:)         ! Element-local time-step, (ndomains,maxelems)


        ! Mesh size information
        real(rk),       allocatable :: mesh_size_elem(:,:), mesh_size_vertex(:,:),       &
                                       min_mesh_size_elem(:), min_mesh_size_vertex(:),   &
                                       avg_mesh_size_vertex(:), sum_mesh_size_vertex(:), &
                                       avg_mesh_h_vertex(:,:), sum_mesh_h_vertex(:,:),   &
                                       area_weighted_h(:,:)
        integer(ik),    allocatable :: num_elements_touching_vertex(:)


        ! RBF-related information
        integer(ik),    allocatable :: nelems_per_domain(:)
        real(rk),       allocatable :: rbf_center(:,:), rbf_radius(:,:)

        ! Vertex-based smoothing information
        integer(ik),    allocatable :: nnodes_per_domain(:)

        ! Global nodes array - used for octree operation
        real(rk),       allocatable :: global_nodes(:,:)

        ! Function registration
        type(function_status_t)         :: function_status ! Status of function residuals and linearizations


        logical                         :: solverInitialized = .false.


    contains

        generic, public     :: init => init_base
        procedure, private  :: init_base


        procedure           :: nauxiliary_fields
        procedure           :: add_auxiliary_field
        procedure           :: new_auxiliary_field
        procedure           :: get_auxiliary_field_index
        procedure           :: get_auxiliary_field_name

        procedure           :: set_nelems_per_domain
        procedure           :: set_nnodes_per_domain
        procedure           :: set_global_nodes


        procedure           :: release

    end type solverdata_t
    !*******************************************************************************************







contains







    !>  Initialize solver base data structures
    !!      - allocate and initialize q, dq, rhs, and linearization.
    !!      - Should be called by specialized 'init' procedure for derived solvers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh            Array of mesh_t instances which define storage requirements.
    !!  @param[in]  bcset_coupling  Array of bcset_coupling instances which describe the 
    !!                              coupling of elements in bcs.
    !!  @param[in]  function_data   Array of containers that hold information on number of 
    !!                              each function in eqnset.
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init_base(self,mesh,function_data)
        class(solverdata_t),                intent(inout)   :: self
        type(mesh_t),                       intent(inout)   :: mesh
        type(equationset_function_data_t),  intent(in)      :: function_data(:)
        
        integer(ik) :: ierr, ndom, maxelems, idom, iaux, aux_ID
        logical     :: increase_maxelems = .false.


        ! Create vector/matrix containers
        self%q     = chidg_vector(trim(backend))
        self%dq    = chidg_vector(trim(backend))
        self%rhs   = chidg_vector(trim(backend))
        self%q_in  = chidg_vector(trim(backend))
        self%q_out = chidg_vector(trim(backend))

        self%lhs   = chidg_matrix(trim(backend))

        ! Initialize vectors
        call self%q%init(    mesh,mesh%ntime_)
        call self%dq%init(   mesh,mesh%ntime_)
        call self%rhs%init(  mesh,mesh%ntime_)
        call self%q_in%init( mesh,mesh%ntime_)
        call self%q_out%init(mesh,mesh%ntime_)

        ! Initialize matrix and parallel recv data
        call self%lhs%init(mesh,'full')
        call self%lhs%init_recv(self%rhs)

        ! By default, create 5 auxiliary field vectors. Each initialized with 'empty' field string.
        do iaux = 1,5
            aux_ID = self%new_auxiliary_field()
            call self%auxiliary_field(aux_ID)%init(mesh,mesh%ntime_)
        end do

        ! Find maximum number of elements in any domain
        ndom = mesh%ndomains()
        maxelems = 0
        do idom = 1,ndom
            increase_maxelems = ( mesh%domain(idom)%nelem > maxelems )
            if (increase_maxelems) then
                maxelems = mesh%domain(idom)%nelem
            end if
        end do


        ! Allocate timestep storage
        if (allocated(self%dt)) deallocate(self%dt)
        allocate(self%dt(ndom,maxelems),stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initialize storage on flux and linearization registration
        call self%function_status%init( mesh, function_data)

        ! Confirm solver initialization
        self%solverInitialized = .true.


    end subroutine init_base
    !******************************************************************************************









    !>  Add chidg_vector for storing an auxiliary field for the problem.
    !!
    !!  One could call this as:
    !!      call solverdata%add_auxiliary_field('my field')
    !!
    !!  which would just add an empty chidg_vector for storing the auxiliary field
    !!
    !!  One could also call this as:
    !!      call solverdata%add_auxiliary_field('my field', my_vector)
    !!
    !!  which would create space for a new auxiliary vector and assign the incoming 
    !!  chidg_vector, my_vector, to the field.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function add_auxiliary_field(self,fieldname,auxiliary_vector) result(aux_ID)
        class(solverdata_t),    intent(inout)           :: self
        character(*),           intent(in)              :: fieldname
        type(chidg_vector_t),   intent(in), optional    :: auxiliary_vector

        integer(ik) :: aux_ID

        ! Try and find 'empty' auxiliary vector
        aux_ID = self%get_auxiliary_field_index('empty')
        if (aux_ID == NO_ID) aux_ID = self%new_auxiliary_field()

        ! Set field name
        self%auxiliary_field_name(aux_ID) = string_t(trim(fieldname))

        ! Store incoming vector if present
        if (present(auxiliary_vector)) self%auxiliary_field(aux_ID) = auxiliary_vector

    end function add_auxiliary_field
    !*****************************************************************************************








    !>  Allocate storage for a new auxiliary vector. Return index of new vector 
    !!  in self%auxiliary_vectors(:)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/24/2019
    !!
    !!
    !------------------------------------------------------------------------------------------
    function new_auxiliary_field(self) result(aux_ID)
        class(solverdata_t),    intent(inout)   :: self

        integer(ik) :: naux_vectors, ierr, aux_ID, iaux

        type(string_t),         allocatable :: temp_names(:)
        type(chidg_vector_t),   allocatable :: temp_vectors(:)


        ! Get new size for self%auxiliary_field(:)
        if (allocated(self%auxiliary_field)) then
            naux_vectors = size(self%auxiliary_field) + 1
        else
            naux_vectors = 1
        end if
        ! New auxiliary ID is last entry be default.
        aux_ID = naux_vectors


        ! Allocate temp storage for new vectors and names 
        allocate(temp_names(naux_vectors), &
                 temp_vectors(naux_vectors), stat=ierr)
        if (ierr /= 0) call AllocationError


        do iaux = 1,size(temp_vectors)
            temp_vectors(iaux) = chidg_vector(trim(backend))
        end do


        ! Copy previously added auxiliary fields to new array
        if (naux_vectors > 1) then
            temp_names(1:size(self%auxiliary_field_name)) = self%auxiliary_field_name(1:size(self%auxiliary_field_name))

            ! TODO: Not working correctly with petsc, maybe due to elemental assignment not defined?
            !temp_vectors(1:size(self%auxiliary_field))    = self%auxiliary_field(1:size(self%auxiliary_field))
            ! Fix:
            do iaux = 1,size(self%auxiliary_field)
                temp_vectors(iaux) = self%auxiliary_field(iaux)
            end do
        end if


        ! Set field name: default = 'empty'
        temp_names(aux_ID) = string_t('empty')
        temp_vectors(aux_ID) = chidg_vector(trim(backend))


        ! Move resized temp allocation back to solverdata_t container.
        call move_alloc(temp_names,self%auxiliary_field_name)
        call move_alloc(temp_vectors,self%auxiliary_field)


    end function new_auxiliary_field
    !******************************************************************************************








    !>  Given the name of an auxiliary field, return the index of the vector containing 
    !!  the field in self%auxiliary_field. 
    !!
    !!  The field represented as a chidg_vector would then be accessed as:
    !!      solverdata%auxiliary_field(index)
    !!
    !!  If not found, returns 0.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/1/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_auxiliary_field_index(self,fieldname) result(field_index)
        class(solverdata_t),    intent(in)  :: self
        character(*),           intent(in)  :: fieldname

        integer(ik)                 :: field_index, ifield
        
        ! Loop through names to try and find field
        if (allocated(self%auxiliary_field_name)) then
            field_index = NO_ID
            do ifield = 1,size(self%auxiliary_field_name)
                if (self%auxiliary_field_name(ifield)%get() == trim(fieldname)) then
                    field_index = ifield
                    exit
                end if
            end do
        else
            field_index = NO_ID
        end if


    end function get_auxiliary_field_index
    !******************************************************************************************








    !>  Given the index of an auxiliary field, return the name of the field the auxiliary
    !!  vector is representing.
    !!
    !!  The field represented as a chidg_vector would then be accessed as:
    !!      solverdata%auxiliary_field(index)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/21/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_auxiliary_field_name(self,field_index) result(field_name)
        class(solverdata_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: field_name, user_msg
        
        ! Check bounds
        user_msg = "solverdata%get_auxiliary_field_name: Index is out of bounds."
        if (field_index > self%nauxiliary_fields()) call chidg_signal(FATAL,user_msg)


        ! Return the appropriate field name
        field_name = trim(self%auxiliary_field_name(field_index)%get())


    end function get_auxiliary_field_name
    !******************************************************************************************







    !>  Return the number of auxiliary fields that have been added to the solverdata_t
    !!  container.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/21/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function nauxiliary_fields(self) result(nfields)
        class(solverdata_t),    intent(in)  :: self

        integer(ik) :: nfields


        if (allocated(self%auxiliary_field)) then
            nfields = size(self%auxiliary_field)
        else
            nfields = 0
        end if


    end function nauxiliary_fields
    !*****************************************************************************************






    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    07/13/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine set_nelems_per_domain(self,nelems_per_domain)
        class(solverdata_t),    intent(inout)   :: self
        integer(ik),            intent(in)      :: nelems_per_domain(:)

        integer(ik) :: nelements_g, ierr

        self%nelems_per_domain = nelems_per_domain
        nelements_g = sum(nelems_per_domain)
        
        if (allocated(self%rbf_center))         deallocate(self%rbf_center)
        if (allocated(self%rbf_radius))         deallocate(self%rbf_radius)
        if (allocated(self%area_weighted_h))    deallocate(self%area_weighted_h)
        if (allocated(self%mesh_size_elem))     deallocate(self%mesh_size_elem)
        if (allocated(self%min_mesh_size_elem)) deallocate(self%min_mesh_size_elem)
        allocate(self%area_weighted_h(nelements_g, 3), &
                 self%mesh_size_elem(nelements_g, 3),  &
                 self%min_mesh_size_elem(nelements_g), &
                 self%rbf_center(nelements_g, 3),      &
                 self%rbf_radius(nelements_g, 3), stat=ierr)
        if (ierr /= 0) call AllocationError

        self%rbf_center         = ZERO
        self%rbf_radius         = ZERO
        self%area_weighted_h    = ZERO
        self%mesh_size_elem     = ZERO
        self%min_mesh_size_elem = ZERO

    end subroutine set_nelems_per_domain
    !*****************************************************************************************





    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    07/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine set_nnodes_per_domain(self,nnodes_per_domain)
        class(solverdata_t),    intent(inout)   :: self
        integer(ik),            intent(in)      :: nnodes_per_domain(:)

        integer(ik) :: nnodes, ierr

        self%nnodes_per_domain = nnodes_per_domain
        nnodes = sum(nnodes_per_domain)

        if (allocated(self%mesh_size_vertex))             deallocate(self%mesh_size_vertex)
        if (allocated(self%min_mesh_size_vertex))         deallocate(self%min_mesh_size_vertex)
        if (allocated(self%avg_mesh_size_vertex))         deallocate(self%avg_mesh_size_vertex)
        if (allocated(self%sum_mesh_size_vertex))         deallocate(self%sum_mesh_size_vertex)
        if (allocated(self%avg_mesh_h_vertex))            deallocate(self%avg_mesh_h_vertex)
        if (allocated(self%sum_mesh_h_vertex))            deallocate(self%sum_mesh_h_vertex)
        if (allocated(self%num_elements_touching_vertex)) deallocate(self%num_elements_touching_vertex)
        allocate(self%mesh_size_vertex(nnodes, 3),  &
                 self%min_mesh_size_vertex(nnodes), &
                 self%avg_mesh_size_vertex(nnodes), &
                 self%sum_mesh_size_vertex(nnodes), &
                 self%avg_mesh_h_vertex(nnodes,3),  &
                 self%sum_mesh_h_vertex(nnodes,3),  &
                 self%num_elements_touching_vertex(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError

        self%mesh_size_vertex     = ZERO
        self%min_mesh_size_vertex = ZERO
        self%avg_mesh_size_vertex = ZERO
        self%sum_mesh_size_vertex = ZERO
        self%avg_mesh_h_vertex    = ZERO
        self%sum_mesh_h_vertex    = ZERO
        self%num_elements_touching_vertex = 0

    end subroutine set_nnodes_per_domain
    !*****************************************************************************************






    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/30/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine set_global_nodes(self,global_nodes)
        class(solverdata_t),    intent(inout)   :: self
        real(rk),               intent(in)      :: global_nodes(:,:)

        self%global_nodes = global_nodes

    end subroutine set_global_nodes
    !*****************************************************************************************








    !>  Release allocated data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine release(self)
        class(solverdata_t), intent(inout)   :: self 

        integer(ik) :: iaux

        ! Release chidg_vector data
        call self%q%release()
        call self%dq%release()
        call self%q_in%release()
        call self%q_out%release()
        call self%rhs%release()

        ! Release chidg_matrix data
        call self%lhs%release()

        ! Release auxiliary_field vectors
        do iaux = 1,self%nauxiliary_fields()
            call self%auxiliary_field(iaux)%release()
        end do
        if (allocated(self%auxiliary_field)) deallocate(self%auxiliary_field)

    end subroutine release
    !****************************************************************************************










end module type_solverdata
