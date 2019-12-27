module type_geometry_cache
#include<messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, NO_ID, dQ_DIFF, dX_DIFF, NO_DIFF, dQ_DIFF
    use mod_io,                 only: backend
    use type_mesh,              only: mesh_t
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
    use type_ivector,           only: ivector_t
    use type_svector,           only: svector_t
    use type_integral_cache,    only: integral_cache_t
    use type_function_info,     onlY: function_info_t

    use mod_chidg_mpi,         only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mpi_f08,               only: MPI_Barrier, MPI_Bcast, MPI_AllReduce, MPI_Gather,   &
                                     MPI_INTEGER4, MPI_REAL8, MPI_LOGICAL, MPI_CHARACTER, &
                                     MPI_SUM, MPI_LOR
    use DNAD_D
    implicit none


    !>  Container to store integral computed on a geometry. 
    !!
    !!  'reference' or 'auxiliary'
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/23/2018
    !!
    !---------------------------------------------------------------------------------------------------
    type, public :: geometry_cache_t
        
        ! Storage of each integral to compute on the geometry 
        type(integral_cache_t), allocatable :: integral_cache(:)
        
        ! Lists of elements/faces to compute the geometry on
        type(ivector_t)                     :: idomain_g
        type(ivector_t)                     :: idomain_l
        type(ivector_t)                     :: ielement_g
        type(ivector_t)                     :: ielement_l
        type(ivector_t)                     :: iface

        ! Number of entities and integrals
        integer(ik)                         :: nentities = 0
        integer(ik)                         :: nintegrals = 0
        integer(ik)                         :: integral_cache_size = 5

    contains
        
        ! Initialization procedures 
        procedure,  public  :: initialize => init   
        procedure,  private :: push_back            

        ! Memory access
        procedure,  public  :: set_value
        procedure,  public  :: get_value
        procedure,  public  :: update_global
        procedure,  public  :: get_real
        procedure,  public  :: get_deriv
        
        ! Parallel communication
        procedure,  public  :: comm

        procedure,  public  :: get_id

    end type geometry_cache_t
    !***************************************************************************************************


contains




    !>  Initialize geometry_cache 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/23/2018
    !!
    !!  @param[inout]      mesh             chidg mesh
    !!  @param[in]         geometries       geometries to compute integrals on 
    !!  @param[in]         integral_type    type of integral to compute 
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init(self,mesh,geometries,integral_type,differentiate)
        class(geometry_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(svector_t),            intent(in)      :: geometries
        character(*),               intent(in)      :: integral_type
        integer(ik),                intent(in)      :: differentiate

        integer(ik)                 :: ngeoms, ierr, domain_ID, nelems, group_ID, igeom, ielem,   &
                                       idom_g, idom_l, ielem_g, ielem_l, iface, ipatch, npatches, &
                                       nfaces, patch_iface
        logical                     :: volume_integral = .false.
        logical                     :: surface_integral = .false.
        character(:),   allocatable :: domain_name, group_name


        ! Overall number of geometries 
        ngeoms = geometries%size()


        ! Choose the type of integral
        volume_integral   = ( integral_type == 'VOLUME INTEGRAL')
        surface_integral  = ( integral_type == 'FACE INTEGRAL')


        ! If volume integral, initialize geometry_cache
        if (volume_integral) then
            
            ! Each geometry is a mesh domain
            do igeom = 1,ngeoms
                domain_name = geometries%data_(igeom)%get()
                domain_ID   = mesh%get_domain_id(domain_name)
                if (domain_ID == NO_ID) cycle

                nelems      = mesh%domain(domain_ID)%get_nelements_local()
                do ielem = 1,nelems
                    associate( elem => mesh%domain(domain_ID)%elems(ielem) )
                        idom_g  = elem%idomain_g
                        idom_l  = elem%idomain_l
                        ielem_g = elem%ielement_g
                        ielem_l = elem%ielement_l
                        iface   = NO_ID 
                        call self%push_back(idom_g,idom_l,ielem_g,ielem_l,iface)
                    end associate
                end do !ielem
            end do !igeom
        
        end if



        ! If surface integral, initialize geometry_cache
        if (surface_integral) then

            ! Each geometry is a patch group
            do igeom = 1,ngeoms
                group_name = geometries%data_(igeom)%get()
                group_ID   = mesh%get_bc_patch_group_id(group_name)
                if (group_ID == NO_ID) cycle
                npatches = mesh%bc_patch_group(group_ID)%npatches()
                do ipatch = 1,npatches

                    nfaces = mesh%bc_patch_group(group_ID)%patch(ipatch)%nfaces()
                    do patch_iface = 1,nfaces
                        associate ( face => mesh%bc_patch_group(group_ID)%patch(ipatch) ) 
                            idom_g  = face%idomain_g()
                            idom_l  = face%idomain_l()
                            ielem_g = face%ielement_g(patch_iface)
                            ielem_l = face%ielement_l(patch_iface)
                            iface   = face%iface(patch_iface) 
                            call self%push_back(idom_g,idom_l,ielem_g,ielem_l,iface)
                        end associate
                    end do !patch_iface
                end do !ipatch
            end do !igeom

        end if


        ! Start allocating storage for 5 integrals, if more integral will be needed this will be reallocated
        if ( allocated(self%integral_cache) ) deallocate (self%integral_cache)
        allocate ( self%integral_cache(self%integral_cache_size), stat=ierr )
        if (ierr /= 0) call AllocationError


    end subroutine init
    !***************************************************************************************************






    !>  Push back indicies for new geometry entity and check if it has already been added 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/25/2018
    !!
    !!  @param[inout]      idom_g             
    !!  @param[inout]      idom_l           
    !!  @param[inout]      ielem_g         
    !!  @param[inout]      ielem_l        
    !!  @param[inout]      iface         
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine push_back(self,idom_g,idom_l,ielem_g,ielem_l,iface)
        class(geometry_cache_t),  intent(inout)   :: self
        integer(ik),                intent(inout)   :: idom_g      
        integer(ik),                intent(inout)   :: idom_l      
        integer(ik),                intent(inout)   :: ielem_g      
        integer(ik),                intent(inout)   :: ielem_l      
        integer(ik),                intent(inout)   :: iface 
  
        logical             :: duplicate 
        integer(ik)         :: idomain_g, idomain_l, ielement_g, ielement_l, iface_, &
                               item

        duplicate = .false.

        ! Check if the same instance has already been addded
        do item = 1,self%nentities
            
            idomain_g  = self%idomain_g%at(item) 
            idomain_l  = self%idomain_l%at(item) 
            ielement_g = self%ielement_g%at(item) 
            ielement_g = self%ielement_l%at(item) 
            iface_     = self%iface%at(item) 
            
            duplicate = ( idomain_g  == idom_g  .and. &
                          idomain_l  == idom_l  .and. &
                          ielement_g == ielem_g .and. &
                          ielement_l == ielem_l .and. &
                          iface_     == iface )
            
            if (duplicate) exit

        end do
        
        ! If there is not a duplicate, push back the new indeces
        if (.not. duplicate) then

            call self%idomain_g%push_back(idom_g)
            call self%idomain_l%push_back(idom_l)
            call self%ielement_g%push_back(ielem_g)
            call self%ielement_l%push_back(ielem_l)
            call self%iface%push_back(iface)
        
            ! Increment number of entities in the geometry (auxiliary or reference)
            self%nentities = self%nentities + 1
        
        end if

    end subroutine push_back
    !***************************************************************************************************






    !>  Set Value of integral evaluated either on the single geometry entity or on the overall geometry.
    !!  If it is a local integral (computed on a face or element) then the real value will be added
    !!  toward the overall summmation of all face/element integral belonging to the geometry, whereas
    !!  the derivatives wrt to the current element will be simply stored.
    !!  If it is a overall integral, then the real value and derivatives will be easily stored
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_value(self,mesh,int_name,value_add,vec_model,dtype,fcn_info)
        class(geometry_cache_t),               intent(inout)   :: self
        type(mesh_t),                          intent(in)      :: mesh
        character(*),                          intent(in)      :: int_name
        type(AD_D),                            intent(in)      :: value_add
        type(chidg_vector_t),                  intent(in)      :: vec_model
        integer(ik),                           intent(in)      :: dtype
        type(function_info_t),      optional,  intent(in)      :: fcn_info

        integer(ik)     :: int_ID
        
        ! Find integral 
        int_ID = self%get_id(int_name)


        ! If the integral is not initialized yet, initialize it
        if (int_ID == 0) then
            self%nintegrals = self%nintegrals + 1

            ! Integral cache size is hardocded, check if the max limit is reached
            if (self%nintegrals > self%integral_cache_size) then
                call chidg_signal(FATAL,"type_geometry_cache%set_value: reached maximum number of integral.")
            end if

            int_ID = self%nintegrals
            call self%integral_cache(int_ID)%init(trim(int_name))
        end if

        
        ! Set 
        if ( present(fcn_info) ) then
            call self%integral_cache(int_ID)%set_value(mesh,value_add,vec_model,fcn_info)
        else
            call self%integral_cache(int_ID)%set_value(mesh,value_add,vec_model,dtype)
        end if


    end subroutine set_value
    !***************************************************************************************************







    !>  
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_value(self,mesh,int_name,vec_model,dtype,fcn_info) result(integral)
        class(geometry_cache_t),               intent(inout)   :: self
        type(mesh_t),                          intent(in)      :: mesh
        character(*),                          intent(in)      :: int_name
        type(chidg_vector_t),                  intent(in)      :: vec_model
        integer(ik),                           intent(in)      :: dtype
        type(function_info_t),      optional,  intent(in)      :: fcn_info

        integer(ik)     :: int_ID
        type(AD_D)      :: integral
        
        ! Find integral 
        int_ID = self%get_id(int_name)


        ! If the integral is not initialized yet, initialize it
        if (int_ID == 0) then
            call chidg_signal_one(FATAL,"type_geometry_cache%get_value: integral has not been computed yet.", trim(int_name))
        end if


        ! Set 
        if ( present(fcn_info) ) then
            integral = self%integral_cache(int_ID)%get_value(mesh,fcn_info)
        else
            integral = self%integral_cache(int_ID)%get_value(mesh,vec_model,dtype)
        end if

    end function get_value
    !***************************************************************************************************






    !>  This is called to update the overall integral real value with the "reduced all" value comingg from
    !!  proc communication 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine update_global(self,int_name,value_rk,vec_model,dtype)
        class(geometry_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: int_name
        real(rk),                   intent(in)      :: value_rk
        type(chidg_vector_t),       intent(in)      :: vec_model
        integer(ik),                intent(in)      :: dtype

        integer(ik)     :: int_ID
        
        ! Loop through integrals to check it already exists
        int_ID = self%get_id(int_name)
        

        ! If the integral is not initialized yet, initialize it and set derivatives to zero.
        if (int_ID == 0) then
            self%nintegrals = self%nintegrals + 1

            ! Integral cache size is hardcoded, check if the max limit is reached
            if (self%nintegrals > self%integral_cache_size) then
                call chidg_signal(FATAL,"type_geometry_cache%set_value: reached maximum number of integral.")
            end if

            int_ID = self%nintegrals
            call self%integral_cache(int_ID)%init(trim(int_name))

            ! If dtype /= NO_DIFF, initialize derivatives and set them to zero.
            ! This is need if a proc has no entities from this geometry
            if (dtype /= NO_DIFF) then
                self%integral_cache(int_ID)%integral_deriv = vec_model
                call self%integral_cache(int_ID)%integral_deriv%assemble()
                self%integral_cache(int_ID)%derivatives_initialized = .true.
            end if

        end if
        
        
        ! Update to global value 
        self%integral_cache(int_ID)%integral_value = value_rk


    end subroutine update_global
    !***************************************************************************************************





    !> Get integral ID by name
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_id(self,integral_name) result(int_ID)
        class(geometry_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: integral_name

        integer(ik)     :: int_ID, iint

        int_ID = 0
        do iint = 1,self%nintegrals
            if (self%integral_cache(iint)%name == trim(integral_name)) then
                int_ID = iint
                exit
            end if
        end do

    end function get_id
    !***************************************************************************************************







    !>  Get real value of the integral by name 
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_real(self,integral_name) result(integral)
        class(geometry_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: integral_name

        real(rk)     :: integral
        integer(ik)  :: int_ID

        ! Get integral ID 
        int_ID = self%get_id(integral_name)

        ! If the integral is not initialized, call error
        if (int_ID == 0) then
            call chidg_signal_one(FATAL,"type_geometry_cache%get_real: integral has not been computed yet.",trim(integral_name))
        end if 

        integral = self%integral_cache(int_ID)%integral_value

    end function get_real
    !***************************************************************************************************







    !>  Get real value of the integral by name 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_deriv(self,integral_name) result(integral)
        class(geometry_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: integral_name

        type(chidg_vector_t)    :: integral
        integer(ik)             :: int_ID

        ! Get integral ID 
        int_ID = self%get_id(integral_name)


        ! If the integral is not initialized, call error
        if (int_ID == 0) then
            call chidg_signal_one(FATAL,"type_geometry_cache%get_deriv: integral has not been computed yet.",trim(integral_name))
        end if

        ! If the integral is not initialized, call error
        if ( .not. self%integral_cache(int_ID)%derivatives_initialized ) then
            call chidg_signal_one(FATAL,"type_geometry_cache%get_deriv: integral has no derivatives.",trim(integral_name))
        end if

        ! Assemble derivatives
        call self%integral_cache(int_ID)%integral_deriv%assemble()
        
        integral = chidg_vector(trim(backend))
        integral = self%integral_cache(int_ID)%integral_deriv

        call integral%assemble()

    end function get_deriv
    !***************************************************************************************************






    
    !>  Parallel communication of processors to share and sum up the functional computed locally
    !!  TODO: using collective blocking MPI procedures, maybe using non-blocking will make it faster
    !!  TODO: can be simplified and improved!!!!
    !!
    !!  Only the real part of the integral is shared!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine comm(self,vec_model,dtype)
        class(geometry_cache_t),    intent(inout)   :: self
        type(chidg_vector_t),       intent(in)      :: vec_model
        integer(ik),                intent(in)      :: dtype


        integer(ik)                 :: local_entities, ientity, leader_rank, global_nint, &
                                       i_int, ierr, str_length
        logical                     :: local_cache_found, global_store
        logical,      allocatable   :: receive_info(:)
        real(rk)                    :: local_value, global_value 
        integer(ik)                 :: integral_id
        character(:), allocatable   :: integral_name
        character(50)               :: string_to_send

        ! Synchronize
        call MPI_Barrier(ChiDG_COMM,ierr)        

        
        ! Check if at least one processor has computed the functional
        ! NOTE: It might be that the entire functional geometry is handled by one single 
        !       processor. This means that all other processors do not carry out any functional
        !       computation
        local_entities = self%nentities
        local_cache_found = (local_entities /= 0)
        call MPI_AllReduce(local_cache_found,global_store,1,MPI_LOGICAL,MPI_LOR,ChiDG_COMM,ierr)


        ! If at least a processor have computed something, share information to store
        ! If the auxiliary geometry is not needed, global_store will be false
        if (global_store) then
        
            ! Master processor allocate number of flags. Flags = .true. means that the irank has a cache
            ! that can be used to share information regarding the functional
            if (IRANK==GLOBAL_MASTER) then
                allocate(receive_info(NRANK), stat=ierr)    
                if (ierr /= 0) call AllocationError
            end if

            
            ! The master processor receive from all the ranks  
            call MPI_Gather(local_cache_found,1,MPI_LOGICAL,receive_info,1,MPI_LOGICAL,GLOBAL_MASTER,ChiDG_COMM,ierr)

            
            ! The root processor loops throught receive_info buffer and return the rank ID of the first
            ! processor that contains functional info (leader_rank)
            ! NOTE: here we want to make sure that the leader_rank is one that contains integrals to store
            if (IRANK==GLOBAL_MASTER) then
                do leader_rank = 0,NRANK-1
                    if (receive_info(leader_rank+1)) then
                        exit
                    end if
                end do
            end if
            call MPI_Bcast(leader_rank,1,MPI_INTEGER4,GLOBAL_MASTER,ChiDG_COMM,ierr)

            
            ! The leader_rank broadcast the number of integrals necessary for this functional
            if (IRANK == leader_rank) then
                global_nint = self%nintegrals
            end if
            call MPI_Bcast(global_nint,1,MPI_INTEGER4,leader_rank,ChiDG_COMM,ierr)
           
           
            ! Loop through the integrals, share name and reduce_all
            do i_int = 1,global_nint
                
                ! Send integral name
                if (IRANK == leader_rank) then
                    string_to_send = self%integral_cache(i_int)%name
                end if
                call MPI_Bcast(string_to_send,50,MPI_CHARACTER,leader_rank,ChiDG_COMM,ierr)

           
                ! Find integral value by name 
                integral_name = trim(string_to_send)
                integral_id   = self%get_id(integral_name)
                if (integral_id /= 0) then
                    local_value   = self%integral_cache(integral_id)%integral_value
                else
                    local_value   = ZERO 
                end if

                ! Sum reduction and broadcast data to all processors
                call MPI_AllReduce(local_value,global_value,1,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)

                ! Update to global value
                call self%update_global(integral_name,real(global_value,rk),vec_model,dtype)

            end do

        end if 
           
    end subroutine comm
    !***************************************************************************************************





end module type_geometry_cache
