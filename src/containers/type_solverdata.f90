module type_solverdata
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NFACES
    use mod_string,                     only: string_t
    use type_chidg_vector,               only: chidg_vector_t
    use type_chidg_matrix,               only: chidg_matrix_t
    use type_mesh,                      only: mesh_t
    use type_function_status,           only: function_status_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_element_info,              only: element_info_t
    use type_face_info,                 only: face_info_t
    use type_function_info,             only: function_info_t
    implicit none



    !> Container for solver data.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   3/15/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public  :: solverdata_t

        !
        ! Base solver data
        !
        type(chidg_vector_t)             :: q              !< Solution vector
        type(chidg_vector_t)             :: dq             !< Change in solution vector
        type(chidg_vector_t)             :: rhs            !< Residual of the spatial scheme
        type(chidg_matrix_t)             :: lhs            !< Linearization of the spatial scheme


        !
        ! Container for reading data
        !
        type(chidg_vector_t)            :: q_in

        !
        ! Auxiliary fields
        !
        type(string_t),         allocatable :: auxiliary_field_name(:)
        type(chidg_vector_t),   allocatable :: auxiliary_field(:)

        !
        ! Time information
        !
        real(rk)                        :: t               !< Global time
        real(rk),   allocatable         :: dt(:,:)         !< Element-local time-step, (ndomains,maxelems)

        !
        ! Function registration
        !
        type(function_status_t)         :: function_status !< Status of function residuals and linearizations


        logical                         :: solverInitialized = .false.




        ! NOTE: if one wanted to add specialized data, instead of deriving from chidgData, 
        ! maybe you could add a chidgExtension class that could be specialized further which 
        ! could contain non-standard data class(chidgExtension_t)

    contains

        generic, public     :: init => init_base
        procedure, private  :: init_base

        procedure           :: add_auxiliary_field

        procedure           :: nauxiliary_fields
        procedure           :: get_auxiliary_field_index
        procedure           :: get_auxiliary_field_name

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
    !subroutine init_base(self,mesh,bcset_coupling,function_data)
    !subroutine init_base(self,mesh,bc,function_data)
    subroutine init_base(self,mesh,function_data)
        class(solverdata_t),                intent(inout)           :: self
        type(mesh_t),                       intent(inout)           :: mesh(:)
!        type(bc_t),                         intent(inout)           :: bc(:)
!        type(bcset_coupling_t),             intent(in)              :: bcset_coupling(:)
        type(equationset_function_data_t),  intent(in)              :: function_data(:)
        

        integer(ik) :: ierr, ndom, maxelems, idom
        logical     :: increase_maxelems = .false.


        ! Initialize and allocate storage
        call self%q%init(  mesh)
        call self%dq%init( mesh)
        call self%rhs%init(mesh)
!        call self%lhs%init(mesh,bcset_coupling,'full')
        call self%lhs%init(mesh,'full')
        call self%q_in%init(mesh)

        ! Initialize matrix parallel recv data
        call self%lhs%init_recv(self%rhs)




    
        !
        ! Find maximum number of elements in any domain
        !
        ndom = size(mesh)
        maxelems = 0
        do idom = 1,ndom

            increase_maxelems = ( mesh(idom)%nelem > maxelems )

            if (increase_maxelems) then
                maxelems = mesh(idom)%nelem
            end if

        end do



        !
        ! Allocate timestep storage
        !
        if (allocated(self%dt)) deallocate(self%dt)
        allocate(self%dt(ndom,maxelems),stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize storage on flux and linearization registration
        !
        call self%function_status%init( mesh, function_data)

        

        !
        ! Confirm solver initialization
        !
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
    subroutine add_auxiliary_field(self,fieldname,auxiliary_vector)
        class(solverdata_t),    intent(inout)           :: self
        character(*),           intent(in)              :: fieldname
        type(chidg_vector_t),   intent(in), optional    :: auxiliary_vector

        integer(ik) :: naux_vectors, ierr

        type(string_t),         allocatable :: temp_names(:)
        type(chidg_vector_t),   allocatable :: temp_vectors(:)


        !
        ! Resize array storage
        !
        if (allocated(self%auxiliary_field)) then
            naux_vectors = size(self%auxiliary_field) + 1
        else
            naux_vectors = 1
        end if

        
        allocate(temp_names(naux_vectors), &
                 temp_vectors(naux_vectors), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Copy previously added auxiliary fields to new array
        !
        if (naux_vectors > 1) then
            temp_names(1:size(self%auxiliary_field_name)) = self%auxiliary_field_name(1:size(self%auxiliary_field_name))
            temp_vectors(1:size(self%auxiliary_field))    = self%auxiliary_field(1:size(self%auxiliary_field))
        end if



        !
        ! Set field name
        !
        temp_names(naux_vectors) = string_t(trim(fieldname))



        !
        ! Move resized temp allocation back to solverdata_t container.
        !
        call move_alloc(temp_names,self%auxiliary_field_name)
        call move_alloc(temp_vectors,self%auxiliary_field)


        !
        ! Store incoming vector if present
        !
        if (present(auxiliary_vector)) then
            self%auxiliary_field(naux_vectors) = auxiliary_vector
        end if

    end subroutine add_auxiliary_field
    !*****************************************************************************************









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
        

        !
        ! Loop through names to try and find field
        !
        if (allocated(self%auxiliary_field_name)) then

            field_index = 0
            do ifield = 1,size(self%auxiliary_field_name)
                if (self%auxiliary_field_name(ifield)%get() == trim(fieldname)) then
                    field_index = ifield
                    exit
                end if
            end do

        else

            field_index = 0

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













    !>  Release allocated data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine release(self)
       class(solverdata_t), intent(inout)   :: self 

        ! Release chidg_vector data
        call self%q%release()
        call self%dq%release()
        call self%q_in%release()
        call self%rhs%release()

        ! Release chidg_matrix data
        call self%lhs%release()


    end subroutine release
    !****************************************************************************************










end module type_solverdata
