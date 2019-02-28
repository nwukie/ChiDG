module type_chidg_vector
#include <messenger.h>
#include "petsc/finclude/petscvec.h"
    use petscvec,                   only: PETSC_DETERMINE, VecCreate, VecSetType, VecSetSizes, VecSetUp,                &
                                          VecSetValues, tVec, tVecScatter, ADD_VALUES, INSERT_VALUES, VecCopy,          &
                                          VecAssemblyBegin, VecAssemblyEnd, VecDuplicate, NORM_2,VecGetArrayF90,        &
                                          VecRestoreArrayF90, VecNorm, VecScale, VecWAXPY, VecReciprocal, VecDestroy,   &
                                          VecScatterCreateToAll, VecScatterBegin, VecScatterEnd, SCATTER_FORWARD, VecScatterDestroy

    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: ZERO, TWO, ONE, NO_ID, NO_DATA
    use mod_chidg_mpi,              only: GROUP_MASTER, ChiDG_COMM, IRANK
    use type_mesh,                  only: mesh_t
    use type_face_info,             only: face_info_t
    use type_element_info,          only: element_info_t
    use type_function,              only: function_t
    use type_chidg_vector_send,     only: chidg_vector_send_t
    use type_chidg_vector_recv,     only: chidg_vector_recv_t
    use type_domain_vector
    use mpi_f08,                    only: MPI_AllReduce, MPI_Reduce, MPI_COMM, MPI_REAL8,    &
                                          MPI_SUM, MPI_STATUS_IGNORE, MPI_Recv, MPI_Request, &
                                          MPI_STATUSES_IGNORE, MPI_INTEGER4
    implicit none





    !>  High-level ChiDG vector container.
    !! 
    !!  Container stores a domain_vector_t for each domain_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chidg_vector_t

        ! PETSC
        Vec         :: petsc_vector
        Vec         :: petsc_vector_recv
        VecScatter  :: petsc_scatter
        logical     :: petsc_vector_created = .false.


        ! ChiDG
        type(domain_vector_t),    allocatable   :: dom(:)       ! Local block vector storage

        type(chidg_vector_send_t)               :: send         ! What to send to other processors
        type(chidg_vector_recv_t)               :: recv         ! Receive data from other processors



        ! backend dynamic procedures
        procedure(vector_init_interface),      pointer, pass   :: init             => chidg_init_vector
        procedure(vector_store_interface),     pointer, pass   :: set_field        => chidg_set_field
        procedure(vector_store_interface),     pointer, pass   :: add_field        => chidg_add_field
        procedure(vector_getfield_interface),  pointer, pass   :: select_get_field => chidg_get_field
        procedure(vector_self_interface),      pointer, pass   :: clear            => chidg_clear_vector
        procedure(vector_self_interface),      pointer, pass   :: assemble         => chidg_assemble_vector
        procedure(vector_assign_interface),    pointer, nopass :: assign_vector    => chidg_assign_vector

        integer(ik),    private                 :: ntime_       ! No. of time instances stored

    contains


        procedure,  public  :: project                          ! Project function to basis

        procedure, public :: get_field

        generic,    public  :: norm => norm_local, norm_comm    ! Compute L2 vector norm
        procedure,  public  :: norm_local                       ! proc-local L2 vector norm
        procedure,  public  :: norm_comm                        ! MPI group L2 vector norm
        generic,    public  :: norm_fields => norm_fields_comm  ! L2 norm of independent fields
        procedure,  public  :: norm_fields_comm                 ! MPI group L2 field norms

        procedure,  public  :: sumsqr                           ! Sum squared proc-local entries 
        procedure,  public  :: sumsqr_fields
        procedure,  public  :: dump

        procedure,  public  :: comm_send                        ! Nonblocking send to comm procs
        procedure,  public  :: comm_recv                        ! Blocking recv incomming data
        procedure,  public  :: comm_wait                        ! Wait to finish send data

        procedure,  public  :: release                          ! Release allocated resources
        procedure,  public  :: get_ntime                        ! Return ntime associated with
        procedure,  public  :: set_ntime                        ! Set ntime in the associated
                                                                ! densevectors
        procedure,  public  :: ndomains

        procedure,  public  :: restrict
        procedure,  public  :: prolong
                                                                    

        final               :: destroy

        procedure, public  :: assign_vector_public
        generic :: assignment(=) => assign_vector_public
        

    end type chidg_vector_t
    !*****************************************************************************************










    !------------------------       OPERATORS       --------------------------------------

    public operator (*)
    interface operator(*)
        module procedure mult_real_chidg_vector          ! real * chidg_vector
        module procedure mult_chidg_vector_real          ! chidg_vector * real
    end interface


    public operator (/)
    interface operator (/)
        module procedure div_real_chidg_vector           ! real / chidg_vector
        module procedure div_chidg_vector_real           ! chidg_vector / real
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_chidg_vector_chidg_vector    ! chidg_vector - chidg_vector
        module procedure minus_chidg_vector               ! - chidg_vector
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_chidg_vector_chidg_vector    ! chidg_vector + chidg_vector
    end interface





    interface 
        subroutine vector_init_interface(self,mesh,ntime)
            import chidg_vector_t
            import mesh_t
            import ik
            class(chidg_vector_t),  intent(inout)   :: self
            type(mesh_t),           intent(inout)   :: mesh
            integer(ik),            intent(in)      :: ntime
        end subroutine vector_init_interface
    end interface


    interface 
        subroutine vector_store_interface(self,values,element_info,ifield,itime)
            import chidg_vector_t
            import element_info_t
            import rk
            import ik
            class(chidg_vector_t),  intent(inout)   :: self
            real(rk),               intent(in)      :: values(:)
            type(element_info_t),   intent(in)      :: element_info
            integer(ik),            intent(in)      :: ifield
            integer(ik),            intent(in)      :: itime
        end subroutine vector_store_interface
    end interface


    interface 
        subroutine vector_getfield_interface(self,element_info,ifield,itime,values)
            import chidg_vector_t
            import element_info_t
            import rk
            import ik
            class(chidg_vector_t),  intent(inout),  target  :: self
            type(element_info_t),   intent(in)              :: element_info
            integer(ik),            intent(in)              :: ifield
            integer(ik),            intent(in)              :: itime
            real(rk), allocatable,  intent(inout)           :: values(:)
        end subroutine vector_getfield_interface
    end interface






    interface 
        subroutine vector_self_interface(self)
            import chidg_vector_t
            class(chidg_vector_t),  intent(inout)   :: self
        end subroutine vector_self_interface
    end interface


    interface 
        subroutine vector_assign_interface(vec_out,vec_in)
            import chidg_vector_t
            class(chidg_vector_t),  intent(inout)   :: vec_out
            class(chidg_vector_t),  intent(in)      :: vec_in
        end subroutine vector_assign_interface
    end interface



    interface chidg_vector
        module procedure new_chidg_vector
    end interface







contains


    !>
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    function new_chidg_vector(storage) result(vec)
        character(*),   intent(in)  :: storage

        type(chidg_vector_t)    :: vec

        select case(storage)
            case('native')
                call vector_assign_pointers_chidg(vec)
            case('petsc')
                call vector_assign_pointers_petsc(vec)
            case default
                call chidg_signal_one(FATAL,"new_chidg_vector: invalid parameter for 'storage'.",trim(storage))
        end select

    end function new_chidg_vector
    !***********************************************************************************




    !>  Allocate and initialize chidg_vector_t storage and data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances used to initialize each 
    !!                      domain_vector_t subcomponent.
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_init_vector(self,mesh,ntime)
        class(chidg_vector_t),  intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: ntime

        integer(ik) :: ierr, ndomains, idom


        ! Set ntime_ for the chidg_vector
        self%ntime_ = ntime

        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        if (allocated(self%dom)) deallocate(self%dom)


        ! Allocate domain_vector_t for each mesh
        ndomains = mesh%ndomains()
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Call initialization procedure for each domain_vector_t
        do idom = 1,ndomains
            call self%dom(idom)%init(mesh%domain(idom))
        end do


        ! Call initialization for determining what data to send and where
        call self%send%init(mesh)

        ! Call initialization for determining what data to receive and allocate storage for it
        call self%recv%init(mesh)

        ! Wait on outstanding mpi_reqests initiated during the send%init(mesh) call
        call self%send%init_wait()

    end subroutine chidg_init_vector
    !******************************************************************************************



    !>  Allocate and initialize chidg_vector_t storage and data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances used to initialize each 
    !!                      domain_vector_t subcomponent.
    !!
    !------------------------------------------------------------------------------------------
    subroutine petsc_init_vector(self,mesh,ntime)
        class(chidg_vector_t),  intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: ntime

        integer(ik)     :: ndomains, idom, ielem
        PetscErrorCode  ierr
        PetscInt        nlocal_rows, nglobal_rows

        ! If previously allocated, destroy and reinitialize
        if (self%petsc_vector_created) then
            call VecDestroy(self%petsc_vector,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecDestroy.')
            call VecDestroy(self%petsc_vector_recv,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecDestroy.')
            call VecScatterDestroy(self%petsc_scatter,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecScatterDestroy.')
        end if

        ! Set ntime_ for the chidg_vector
        self%ntime_ = ntime

        ! Create vector object
        call VecCreate(ChiDG_COMM%mpi_val, self%petsc_vector, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init_vector: error creating PETSc vector.')

        ! Set vector type
        call VecSetType(self%petsc_vector, 'standard', ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init_vector: error calling VecSetType.')

        ! Set vector size
        ! Compute proc-local degress-of-freedom
        nlocal_rows = 0
        do idom = 1,mesh%ndomains()
            do ielem = 1,mesh%domain(idom)%nelements()
                nlocal_rows = nlocal_rows + mesh%domain(idom)%elems(ielem)%nterms_s * mesh%domain(idom)%elems(ielem)%neqns
            end do !ielem
        end do !idom

        ! Compute global degrees-of-freedom via reduction
        call MPI_AllReduce(nlocal_rows,nglobal_rows,1,MPI_INTEGER4,MPI_SUM,ChiDG_COMM,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error reducing global degrees-of-freedom.')

        call VecSetSizes(self%petsc_vector,nlocal_rows,nglobal_rows,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecSetSizes.')

        ! Set up vector
        call VecSetUp(self%petsc_vector,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecSetUp.')

        ! Initialize parallel scatter to all
        call VecScatterCreateToAll(self%petsc_vector, self%petsc_scatter, self%petsc_vector_recv, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecScatterCreateToAll.')

        ! Clear vector storage
        call self%clear()

        ! Indicate petsc vector was created
        self%petsc_vector_created = .true.


    end subroutine petsc_init_vector
    !*************************************************************************************


    


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/30/2019
    !!
    !-------------------------------------------------------------------------------------
    subroutine chidg_set_field(self,values,element_info,ifield,itime)
        class(chidg_vector_t),  intent(inout)   :: self
        real(rk),               intent(in)      :: values(:)
        type(element_info_t),   intent(in)      :: element_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime

        call self%dom(element_info%idomain_l)%vecs(element_info%ielement_l)%setvar(ifield,itime,values)

    end subroutine chidg_set_field
    !**************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/30/2019
    !!
    !-------------------------------------------------------------------------------------
    !subroutine chidg_add_field(self,values,face_info,nterms_s,nfields,ifield,itime)
    subroutine chidg_add_field(self,values,element_info,ifield,itime)
        class(chidg_vector_t),  intent(inout)   :: self
        real(rk),               intent(in)      :: values(:)
        type(element_info_t),   intent(in)      :: element_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime

        real(rk),   allocatable :: vals(:)

        ! Access current field and add new contribution
        vals = self%dom(element_info%idomain_l)%vecs(element_info%ielement_l)%getvar(ifield,itime) + values

        ! Store updated values
        call self%dom(element_info%idomain_l)%vecs(element_info%ielement_l)%setvar(ifield,itime,vals)

    end subroutine chidg_add_field
    !**************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/29/2019
    !!
    !-----------------------------------------------------------------------------------
    subroutine petsc_set_field(self,values,element_info,ifield,itime)
        class(chidg_vector_t),  intent(inout)   :: self
        real(rk),               intent(in)      :: values(:)
        type(element_info_t),   intent(in)      :: element_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime

        PetscErrorCode          :: ierr, i
        PetscInt                :: istart
        PetscInt, allocatable   :: indices(:)

        istart = element_info%dof_start + (ifield-1)*element_info%nterms_s + (element_info%nfields*element_info%nterms_s)*(itime-1)
        indices = [(i, i=istart,(istart+element_info%nterms_s-1),1)]

        ! Decrement by 1 for 0-based indexing
        indices = indices - 1

        call VecSetValues(self%petsc_vector,size(values),indices,values,INSERT_VALUES,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_set_field: error calling VecSetValues.')

    end subroutine petsc_set_field
    !***********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/29/2019
    !!
    !-----------------------------------------------------------------------------------
    subroutine petsc_add_field(self,values,element_info,ifield,itime)
        class(chidg_vector_t),  intent(inout)   :: self
        real(rk),               intent(in)      :: values(:)
        type(element_info_t),   intent(in)      :: element_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime

        PetscErrorCode          :: ierr, i
        PetscInt                :: istart
        PetscInt, allocatable   :: indices(:)

        istart = element_info%dof_start + (ifield-1)*element_info%nterms_s + (element_info%nfields*element_info%nterms_s)*(itime-1)
        indices = [(i, i=istart,(istart+element_info%nterms_s-1),1)]

        ! Decrement by 1 for 0-based indexing
        indices = indices - 1

        call VecSetValues(self%petsc_vector,size(values),indices,values,ADD_VALUES,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_add_vector: error calling VecSetValues.')

    end subroutine petsc_add_field
    !***********************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/30/2019
    !!
    !-----------------------------------------------------------------------------------
    subroutine chidg_get_field(self,element_info,ifield,itime,values) 
        class(chidg_vector_t),  intent(inout), target   :: self
        type(element_info_t),   intent(in)              :: element_info
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime
        real(rk), allocatable,  intent(inout)           :: values(:)

        if (element_info%iproc /= IRANK) then
            values = self%recv%comm(element_info%recv_comm)%dom(element_info%recv_domain)%vecs(element_info%recv_element)%getvar(ifield,itime)
        else
            values = self%dom(element_info%idomain_l)%vecs(element_info%ielement_l)%getvar(ifield,itime)
        end if

    end subroutine chidg_get_field
    !***********************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/30/2019
    !!
    !-----------------------------------------------------------------------------------
    subroutine petsc_get_field(self,element_info,ifield,itime,values)
        class(chidg_vector_t),  intent(inout), target   :: self
        type(element_info_t),   intent(in)              :: element_info
        integer(ik),            intent(in)              :: ifield
        integer(ik),            intent(in)              :: itime
        real(rk), allocatable,  intent(inout)           :: values(:)

        integer(ik)             :: istart, iend

        PetscScalar, pointer :: array(:) => null()
        PetscErrorCode  :: ierr

        ! Get petsc array pointer
        if (self%petsc_vector_created) then
            !call VecGetArrayF90(self%petsc_vector,array,ierr)
            call VecGetArrayF90(self%petsc_vector_recv,array,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_get_field: error calling VecGetArrayF90.')
        else
            call chidg_signal(FATAL,'chidg_vector%petsc_get_field: petsc vector not created.')
        end if

        ! Compute start and end indices for accessing modes of a variable
        istart = element_info%dof_start + (ifield-1)*element_info%nterms_s + (element_info%nfields*element_info%nterms_s)*(itime-1)
        iend = istart + (element_info%nterms_s-1)

        ! Access modes
        values = array(istart:iend)

        ! Restore petsc array
        call VecRestoreArrayF90(self%petsc_vector_recv,array,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_get_field: error calling VecGetArrayF90.')


    end subroutine petsc_get_field
    !***********************************************************************************


    function get_field(self,element_info,ifield,itime) result(values)
        class(chidg_vector_t),  intent(inout)   :: self
        type(element_info_t),   intent(in)      :: element_info
        integer(ik),            intent(in)      :: ifield
        integer(ik),            intent(in)      :: itime

        real(rk),   allocatable :: values(:)

        call self%select_get_field(element_info,ifield,itime,values)

    end function get_field





    !>  Project a function onto the global solution basis.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/25/2016
    !!
    !!  @author Mayank Sharma
    !!  @date   11/12/2016
    !!
    !!  TODO: Should itime be an input parameter here?
    !!
    !------------------------------------------------------------------------------------------
    subroutine project(self,mesh,fcn,ifield)
        class(chidg_vector_t),  intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh
        class(function_t),      intent(inout)   :: fcn
        integer(ik),            intent(in)      :: ifield

        integer(ik)                 :: idom, ielem, ierr
        integer(ik)                 :: itime
        real(rk),       allocatable :: modes(:)
        character(:),   allocatable :: user_msg
        type(element_info_t)    :: element_info


        ! Loop through elements in mesh and call function projection
        do idom = 1,mesh%ndomains()

            ! Check that variable index 'ifield' is valid
            user_msg = 'project: variable index ifield exceeds the number of equations.'
            if (ifield > mesh%domain(idom)%neqns ) call chidg_signal(FATAL,user_msg)

            do ielem = 1,mesh%domain(idom)%nelem

                element_info = element_info_t(idomain_g  = mesh%domain(idom)%elems(ielem)%idomain_g,    &
                                              idomain_l  = mesh%domain(idom)%elems(ielem)%idomain_l,    &
                                              ielement_g = mesh%domain(idom)%elems(ielem)%ielement_g,   &
                                              ielement_l = mesh%domain(idom)%elems(ielem)%ielement_l,   &
                                              iproc      = mesh%domain(idom)%elems(ielem)%iproc,        &
                                              pelem_ID   = NO_ID,                                       &
                                              eqn_ID     = mesh%domain(idom)%elems(ielem)%eqn_ID,       &
                                              nfields    = mesh%domain(idom)%elems(ielem)%neqns,        &
                                              nterms_s   = mesh%domain(idom)%elems(ielem)%nterms_s,     &
                                              nterms_c   = NO_DATA,                                     &
                                              dof_start  = mesh%domain(idom)%elems(ielem)%dof_start)


                do itime = 1,mesh%domain(idom)%ntime
                    ! Call function projection
                    modes = mesh%domain(idom)%elems(ielem)%project(fcn)

                    ! Store the projected modes to the solution expansion
                    call self%set_field(modes,element_info,ifield,itime)

                    call self%assemble()

                    modes  = self%get_field(element_info,ifield,itime)

                end do ! itime
            end do ! ielem
        end do ! idomain

        
        ! Assemble vector 
        call self%assemble()


    end subroutine project
    !******************************************************************************************






    !>  Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_clear_vector(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer :: idom

        ! Call clear procedure for each domain_vector_t
        if (allocated(self%dom)) then
            do idom = 1,size(self%dom)
                call self%dom(idom)%clear()
            end do
        end if

        ! Call clear on recv storage
        call self%recv%clear()

    end subroutine chidg_clear_vector
    !******************************************************************************************





    !>  Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/29/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine petsc_clear_vector(self)
        class(chidg_vector_t),   intent(inout)   :: self

        PetscErrorCode :: perr

        call VecSet(self%petsc_vector,ZERO,perr)
        if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_clear_vector: error calling VecSet.')

    end subroutine petsc_clear_vector
    !******************************************************************************************



    !>  Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/29/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_assemble_vector(self)
        class(chidg_vector_t),   intent(inout)   :: self

        call self%comm_send()
        call self%comm_recv()
        call self%comm_wait()

    end subroutine chidg_assemble_vector
    !******************************************************************************************



    !>  Set all floating-point vector entries to zero.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/29/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine petsc_assemble_vector(self)
        class(chidg_vector_t),   intent(inout)   :: self

        PetscErrorCode ierr

        call VecAssemblyBegin(self%petsc_vector,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assmble_vector: error calling VecAssemblyBegin.')

        call VecAssemblyEnd(self%petsc_vector,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assmble_vector: error calling VecAssemblyEnd.')


        ! Scatter
        call VecSet(self%petsc_vector_recv,ZERO,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assemble_vector: error calling VecSet.')

        call VecScatterBegin(self%petsc_scatter, self%petsc_vector, self%petsc_vector_recv, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assemble_vector: error calling VecScatterBegin.')

        call VecScatterEnd(self%petsc_scatter, self%petsc_vector, self%petsc_vector_recv, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assemble_vector: error calling VecScatterBegin.')


    end subroutine petsc_assemble_vector
    !******************************************************************************************





    !>  Compute the process-local L2-Norm of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_local(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        res = self%sumsqr()

        ! Take the square root of the result
        res = sqrt(res)

    end function norm_local
    !******************************************************************************************









    !>  Compute the L2-Norm of the vector within the space of processors given by the MPI 
    !!  communicator.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/23/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_comm(self,comm) result(norm)
        class(chidg_vector_t),   intent(in)  :: self
        type(mpi_comm),         intent(in)  :: comm

        real(rk)    :: sumsqr, norm
        integer     :: ierr

        PetscReal      :: petsc_norm
        PetscErrorCode :: perr


        if (self%petsc_vector_created) then

            call VecNorm(self%petsc_vector,NORM_2,petsc_norm,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%norm_comm: error calling VecNorm.')
            norm = real(petsc_norm,rk)

        else

            norm = ZERO
            ! Compute sum of the squared elements of the processor-local vector
            sumsqr = self%sumsqr()
            ! Reduce sumsqr across all procs, distribute result back to all
            call MPI_AllReduce(sumsqr,norm,1,MPI_REAL8,MPI_SUM,comm,ierr)
            norm = sqrt(norm)
        end if


    end function norm_comm
    !******************************************************************************************











    !>  Compute the L2-Norm of the vector within the space of processors given by the MPI 
    !!  communicator.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   6/23/2016
    !!
    !!  @return res     L2-norm of the vector
    !!
    !------------------------------------------------------------------------------------------
    function norm_fields_comm(self,comm) result(norm)
        class(chidg_vector_t),  intent(in)  :: self
        type(mpi_comm),         intent(in)  :: comm

        real(rk), allocatable   :: sumsqr(:), norm(:)
        integer     :: ierr

        PetscReal      :: petsc_norm
        PetscErrorCode :: perr

        if (self%petsc_vector_created) then

            allocate(norm(1))
            call VecNorm(self%petsc_vector,NORM_2,petsc_norm,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%norm_fields_comm: error calling VecNorm.')
            norm(1) = real(petsc_norm,rk)

        else

            ! Compute sum of the squared elements of the processor-local vector
            sumsqr = self%sumsqr_fields()

            ! Alloate norm
            norm = sumsqr
            norm = ZERO

            ! Reduce sumsqr across all procs, distribute result back to all
            call MPI_AllReduce(sumsqr,norm,size(sumsqr),MPI_REAL8,MPI_SUM,comm,ierr)

            norm = sqrt(norm)
        end if


    end function norm_fields_comm
    !******************************************************************************************










    !< Return the sum of the squared chidg_vector entries 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/23/2016
    !!
    !!  @return res     sum of the squared chidg_vector entries
    !!
    !------------------------------------------------------------------------------------------
    function sumsqr(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem

        res = ZERO

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr()
        end do ! idom

    end function sumsqr
    !******************************************************************************************








    !< Return the sum of the squared chidg_vector entries for each field independently.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2017
    !!
    !!  @return res     sum of the squared chidg_vector entries
    !!
    !------------------------------------------------------------------------------------------
    function sumsqr_fields(self) result(res)
        class(chidg_vector_t),   intent(in)   :: self

        real(rk),   allocatable :: res(:)
        integer(ik) :: idom, ielem


        ! Allocate size of res based on assumption of same equation set across domains.
        res = self%dom(1)%sumsqr_fields()
        res = ZERO

        ! Loop through domain vectors and compute contribution to vector sum of the squared elements
        do idom = 1,size(self%dom)
            res = res + self%dom(idom)%sumsqr_fields()
        end do ! idom


    end function sumsqr_fields
    !******************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_send(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer(ik)         :: icomm, iproc_send, idom_send, idom, ielem_send, &
                               ielem, ierr, data_size, isend
        type(mpi_request)   :: isend_handle


        ! Loop through comms to send
        isend = 1
        do icomm = 1,size(self%send%comm)

            ! Get processor rank we are sending to
            iproc_send = self%send%comm(icomm)%proc

            ! Loop through domains/elements to send
            do idom_send = 1,self%send%comm(icomm)%dom_send%size()
                idom = self%send%comm(icomm)%dom_send%at(idom_send)
                do ielem_send = 1,self%send%comm(icomm)%elems_send(idom_send)%size()
                    ielem = self%send%comm(icomm)%elems_send(idom_send)%at(ielem_send)

                    ! Post non-blocking send message for the vector data
                    data_size = size(self%dom(idom)%vecs(ielem)%vec)
                    call MPI_ISend(self%dom(idom)%vecs(ielem)%vec, data_size, MPI_REAL8, iproc_send, 0, ChiDG_COMM, isend_handle, ierr)

                    ! Add non-blocking send handle to list of things to wait on
                    self%send%isend_handles(isend) = isend_handle

                    ! Increment send counter
                    isend = isend + 1

                end do !ielem_send
            end do !idom_send

        end do ! icomm


    end subroutine comm_send
    !*****************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(chidg_vector_t),   intent(inout)   :: self

        integer(ik) :: icomm, idom_recv, ielem_recv, proc_recv, data_size, &
                       ierr, dparent_g, eparent_g

        real(rk), allocatable   :: test(:)

        ! Receive data from each communicating processor
        do icomm = 1,size(self%recv%comm)

            ! Get process we are receiving from
            proc_recv = self%recv%comm(icomm)%proc
            
            ! Recv each element chunk
            do idom_recv = 1,size(self%recv%comm(icomm)%dom)
                do ielem_recv = 1,size(self%recv%comm(icomm)%dom(idom_recv)%vecs)

                    data_size = size(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec)
                    call MPI_Recv(self%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%vec, data_size, MPI_REAL8, proc_recv, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                end do ! ielem_recv
            end do ! idom_recv

        end do ! icomm

    end subroutine comm_recv
    !*****************************************************************************************








    !>  Wait for all non-blocking sends to complete.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine comm_wait(self)
        class(chidg_vector_t),   intent(in)  :: self

        integer(ik) :: nwait, ierr

        nwait = size(self%send%isend_handles)
        call MPI_Waitall(nwait, self%send%isend_handles, MPI_STATUSES_IGNORE, ierr)

    end subroutine comm_wait
    !*****************************************************************************************













    !> Dump contents of the vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine dump(self)
        class(chidg_vector_t),   intent(in)   :: self

        integer(ik) :: idom

        ! Loop through domain vectors and compute contribution to vecotr L2-Norm
        do idom = 1,size(self%dom)
            call self%dom(idom)%dump()
        end do ! idom

    end subroutine dump
    !****************************************************************************************







    !>  Release allocated resources.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_vector_t),  intent(inout)   :: self

        PetscErrorCode :: ierr

        if (allocated(self%dom)) deallocate(self%dom)

        if (self%petsc_vector_created) then
            call VecDestroy(self%petsc_vector,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%release: error calling VecDestroy.')
            call VecDestroy(self%petsc_vector_recv,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%release: error calling VecDestroy.')
            call VecScatterDestroy(self%petsc_scatter,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%release: error calling VecScatterDestroy.')
        end if

    end subroutine release
    !****************************************************************************************







    !>  Return ntime
    !!
    !!  @author Mayank Sharma
    !!  @date   3/9/2017
    !!
    !----------------------------------------------------------------------------------------
    function get_ntime(self) result(ntime_out)
        class(chidg_vector_t),  intent(inout)   :: self

        integer(ik)     :: ntime_out

        ! Get ntime 
        ntime_out = self%ntime_

    end function get_ntime
    !****************************************************************************************







    !>  Set ntime in the densevectors
    !!
    !!  @author Mayank Sharma
    !!  @date   3/9/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_ntime(self,ntime)
        class(chidg_vector_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: ntime

        integer(ik)     :: idom, ielem

        self%ntime_ = ntime

        ! Set ntime
        if (.not. self%petsc_vector_created) then
            do idom = 1,size(self%dom)
                do ielem = 1,size(self%dom(idom)%vecs)
                    call self%dom(idom)%vecs(ielem)%set_ntime(ntime)
                end do
            end do
        end if

    end subroutine set_ntime
    !****************************************************************************************




    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !---------------------------------------------------------------------------------------
    function ndomains(self) result(ndomains_)
        class(chidg_vector_t),  intent(in)  :: self

        integer(ik) :: ndomains_

        if (allocated(self%dom)) then
            ndomains_ = size(self%dom)
        else
            ndomains_ = 0
        end if

    end function ndomains
    !***************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !---------------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(chidg_vector_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms_r

        type(chidg_vector_t)    :: restricted
        integer(ik)             :: idom, ierr
        

        restricted%send = self%send                     ! Copy self%send directly
        restricted%recv = self%recv%restrict(nterms_r)  ! Get restricted copy of self%recv


        ! Allocate storage for each domain
        if (allocated(restricted%dom)) deallocate(restricted%dom)
        allocate(restricted%dom(self%ndomains()), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Return restricted domain_vector objects for each domain
        do idom = 1,self%ndomains()
            restricted%dom(idom) = self%dom(idom)%restrict(nterms_r)
        end do !idom

        ! Set ntime
        restricted%ntime_ = self%ntime_

    end function restrict
    !***************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/21/2017
    !!
    !---------------------------------------------------------------------------------------
    function prolong(self,nterms_p) result(prolonged)
        class(chidg_vector_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: nterms_p

        type(chidg_vector_t)    :: prolonged
        integer(ik)             :: idom, ierr
        

        prolonged%send = self%send                     ! Copy self%send directly
        prolonged%recv = self%recv%prolong(nterms_p)   ! Get prolonged copy of self%recv


        ! Allocate storage for each domain
        allocate(prolonged%dom(self%ndomains()), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Return restricted domain_vector objects for each domain
        do idom = 1,self%ndomains()
            prolonged%dom(idom) = self%dom(idom)%prolong(nterms_p)
        end do !idom

        ! Set ntime
        prolonged%ntime_ = self%ntime_

    end function prolong
    !***************************************************************************************













    !-----------------------------------------------------------------------------------------
    !
    !
    !                              Operator Implementations
    !
    !-----------------------------------------------------------------------------------------



    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------------
    function mult_real_chidg_vector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidg_vector_t),   intent(in)  :: right

        type(chidg_vector_t)    :: res
        integer(ik)             :: idom, ndom
        PetscErrorCode          :: perr

        res%ntime_ = right%ntime_

        if (right%petsc_vector_created) then

            ! Copy
            res = right

            ! Scale by inverse
            call VecScale(res%petsc_vector,left,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%mult_chidg_vector_real: error calling VecScale.')

        else


            ndom = size(right%dom)
            allocate(res%dom(ndom))

            do idom = 1,size(right%dom)
                res%dom(idom) = left * right%dom(idom)
            end do

            res%send = right%send
            res%recv = right%recv

        end if

    end function mult_real_chidg_vector
    !*****************************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function mult_chidg_vector_real(left,right) result(res)
        type(chidg_vector_t),   intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidg_vector_t)    :: res
        integer(ik)             :: idom, ndom
        PetscErrorCode          :: perr

        res%ntime_ = left%ntime_

        if (left%petsc_vector_created) then

            ! Copy
            res = left

            ! Scale by inverse
            call VecScale(res%petsc_vector,right,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%mult_chidg_vector_real: error calling VecScale.')

        else


            ndom = size(left%dom)
            allocate(res%dom(ndom))

            do idom = 1,size(left%dom)
                res%dom(idom) = left%dom(idom) * right
            end do

            res%send = left%send
            res%recv = left%recv

        end if

    end function mult_chidg_vector_real
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function div_real_chidg_vector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom
        PetscErrorCode :: perr


        res%ntime_ = right%ntime_

        if (right%petsc_vector_created) then

            ! Copy
            res = right

            ! Get vector reciprocal
            call VecReciprocal(res%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%div_real_chidg_vector: error calling VecReciprocal.')

            ! Scale by real
            call VecScale(res%petsc_vector,left,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%div_real_chidg_vector: error calling VecScale.')

        else

            ndom = size(right%dom)
            allocate(res%dom(ndom))

            do idom = 1,size(right%dom)
                res%dom(idom) = left / right%dom(idom)
            end do


            res%send = right%send
            res%recv = right%recv

        end if


    end function div_real_chidg_vector
    !****************************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function div_chidg_vector_real(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        PetscErrorCode :: perr

        res%ntime_ = left%ntime_

        if (left%petsc_vector_created) then

            ! Copy
            res = left

            ! Scale by inverse
            call VecScale(res%petsc_vector,ONE/right,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%minus_chidg_vector: error calling VecScale.')

        else

            ndom = size(left%dom)

            allocate(res%dom(ndom))

            do idom = 1,size(left%dom)
                res%dom(idom) = left%dom(idom) / right
            end do

            res%send = left%send
            res%recv = left%recv
        end if

    end function div_chidg_vector_real
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function add_chidg_vector_chidg_vector(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        PetscErrorCode :: perr

        res%ntime_ = right%ntime_

        if (left%petsc_vector_created) then

            ! Copy
            res = left
            ! Operation
            call VecWAXPY(res%petsc_vector,ONE,right%petsc_vector,left%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%sub_chidg_vector_chidg_vector: error calling VecAYPX.')

        else

            ndom = size(right%dom)

            allocate(res%dom(ndom))

            do idom = 1,size(left%dom)
                res%dom(idom) = left%dom(idom) + right%dom(idom)
            end do


            res%send = right%send
            res%recv = right%recv

        end if

    end function add_chidg_vector_chidg_vector
    !****************************************************************************************





    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function sub_chidg_vector_chidg_vector(left,right) result(res)
        type(chidg_vector_t),    intent(in)  :: left
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        PetscErrorCode :: perr

        res%ntime_ = right%ntime_

        if (left%petsc_vector_created) then

            ! Copy
            res = left
            ! Operation
            call VecWAXPY(res%petsc_vector,-ONE,right%petsc_vector,left%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%sub_chidg_vector_chidg_vector: error calling VecWAYPX.')

        else

            ndom = size(right%dom)
            allocate(res%dom(ndom))

            do idom = 1,size(left%dom)
                res%dom(idom) = left%dom(idom) - right%dom(idom)
            end do

            res%send = right%send
            res%recv = right%recv

        end if

    end function sub_chidg_vector_chidg_vector
    !****************************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function minus_chidg_vector(right) result(res)
        type(chidg_vector_t),    intent(in)  :: right

        type(chidg_vector_t) :: res
        integer(ik)         :: idom, ndom

        PetscErrorCode :: perr

        res%ntime_ = right%ntime_

        if (right%petsc_vector_created) then

            ! Copy
            res = right
            ! Negate
            call VecScale(res%petsc_vector,-ONE,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_vector%minus_chidg_vector: error calling VecScale.')

        else

            ! ChiDG
            ndom = size(right%dom)
            allocate(res%dom(ndom))
            do idom = 1,size(right%dom)
                res%dom(idom) = (-ONE)*right%dom(idom)
            end do
            res%send   = right%send
            res%recv   = right%recv

        end if

    end function minus_chidg_vector
    !****************************************************************************************




    subroutine chidg_assign_vector(vec_out,vec_in)
        class(chidg_vector_t),  intent(inout)   :: vec_out
        class(chidg_vector_t),  intent(in)      :: vec_in

        if (allocated(vec_in%dom)) vec_out%dom = vec_in%dom
        vec_out%send   = vec_in%send
        vec_out%recv   = vec_in%recv
        vec_out%ntime_ = vec_in%ntime_

        ! Update procedure pointers
        call vector_assign_pointers_chidg(vec_out)

    end subroutine chidg_assign_vector


    subroutine petsc_assign_vector(vec_out,vec_in)
        class(chidg_vector_t),  intent(inout)   :: vec_out
        class(chidg_vector_t),  intent(in)      :: vec_in

        PetscErrorCode ierr


        if (vec_in%petsc_vector_created) then

            ! If already created, only copy
            if (vec_out%petsc_vector_created) then

                call VecCopy(vec_in%petsc_vector,vec_out%petsc_vector,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecCopy.')
                call VecCopy(vec_in%petsc_vector_recv,vec_out%petsc_vector_recv,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecCopy.')

            ! If not already created, duplicate storage, then copy.
            else

                call VecDuplicate(vec_in%petsc_vector,vec_out%petsc_vector,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecDuplicate.')
                call VecDuplicate(vec_in%petsc_vector_recv,vec_out%petsc_vector_recv,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecDuplicate.')

                call VecCopy(vec_in%petsc_vector,vec_out%petsc_vector,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecCopy.')
                call VecCopy(vec_in%petsc_vector_recv,vec_out%petsc_vector_recv,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_assign_vector: error in VecCopy.')

            end if

            vec_out%petsc_vector_created = .true.
        end if

        vec_out%ntime_ = vec_in%ntime_

        ! Update procedure pointers
        call vector_assign_pointers_petsc(vec_out)

    end subroutine petsc_assign_vector


    subroutine assign_vector_public(vec_out,vec_in)
        class(chidg_vector_t),  intent(inout) :: vec_out
        class(chidg_vector_t),  intent(in)  :: vec_in

        call vec_in%assign_vector(vec_out,vec_in)

    end subroutine assign_vector_public




    subroutine vector_assign_pointers_chidg(vec)
        type(chidg_vector_t),   intent(inout)   :: vec

        vec%init             => chidg_init_vector
        vec%add_field        => chidg_add_field
        vec%set_field        => chidg_set_field
        vec%select_get_field => chidg_get_field
        vec%clear            => chidg_clear_vector
        vec%assemble         => chidg_assemble_vector
        vec%assign_vector    => chidg_assign_vector

    end subroutine vector_assign_pointers_chidg

    subroutine vector_assign_pointers_petsc(vec)
        type(chidg_vector_t),   intent(inout)   :: vec

        vec%init             => petsc_init_vector
        vec%clear            => petsc_clear_vector
        vec%set_field        => petsc_set_field
        vec%add_field        => petsc_add_field
        vec%select_get_field => petsc_get_field
        vec%assemble         => petsc_assemble_vector
        vec%assign_vector    => petsc_assign_vector

    end subroutine vector_assign_pointers_petsc






    !>
    !!
    !!
    !!
    !-------------------------------------------------
    subroutine destroy(self)
        type(chidg_vector_t),   intent(inout)   :: self

        PetscErrorCode :: ierr

        if (self%petsc_vector_created) then
            call VecDestroy(self%petsc_vector,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecDestroy.')
            call VecDestroy(self%petsc_vector_recv,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%petsc_init: error calling VecDestroy.')
            call VecScatterDestroy(self%petsc_scatter,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_vector%release: error calling VecScatterDestroy.')
        end if


    end subroutine destroy
    !*************************************************






end module type_chidg_vector
