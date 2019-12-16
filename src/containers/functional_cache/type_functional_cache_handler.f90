module type_functional_cache_handler
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use type_chidg_data,            only: chidg_data_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_cache_handler,         only: cache_handler_t
    use type_element_info,          only: element_info_t
    use mod_constants,              only: DIAG, ZERO, dQ_DIFF, dX_DIFF, dBC_DIFF, NO_DIFF
    use type_functional_cache,      only: functional_cache_t
    use type_evaluator,             only: evaluator_t
    use type_integral_cache,        only: integral_cache_t
    use type_geometry_cache,        only: geometry_cache_t
    use mod_DNAD_tools,             only: face_compute_seed, element_compute_seed

    use mod_chidg_mpi,              only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                     only: verbosity
    use mpi_f08,                    only: MPI_Barrier, MPI_Bcast, MPI_AllReduce, MPI_Gather,   &
                                          MPI_INTEGER4, MPI_REAL8, MPI_LOGICAL, MPI_CHARACTER, &
                                          MPI_SUM, MPI_LOR
    use DNAD_D
    implicit none


    !>  An object to handle functional computation operations. 
    !!
    !!  This allows to temporarily store and handle the data coming from the functional calculation.
    !!  It is meant to avoid the time consuming multiple loops over the elements/faces belonging
    !!  to the functional geometries.
    !!
    !!  It is also useful for easy processor communication in parallel mode.
    !!
    !!  The functional real and derivatives values are store in a sistematic way, independently whether
    !!  it is a surface or a volume type of functional.
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/23/2018
    !!
    !!  Restructured
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    type, public :: functional_cache_handler_t
        
        type(functional_cache_t),   pointer      :: cache            ! Point to the functional cache to be used for storing info 
        type(chidg_data_t),         pointer      :: data             ! ChiDG data
        character(:),               allocatable  :: geom             ! Geometry we are computing integrals on (reference/auxiliary)
        character(:),               allocatable  :: component        ! Component element/face_interior for worker cache update
        integer(ik)                              :: dtype            ! Type of differentiation requested
        integer(ik)                              :: ifunc            ! Functional currently computed

    contains
        
        procedure,  public   :: update       ! Functional computation update
       
        ! Private procedures (made public for tests) 
        procedure,  public   :: init         ! Initialize functional cache handler
        procedure,  public   :: compute      ! Compute functional
        procedure,  public   :: store        ! Share functional integral with other processors
        procedure,  public   :: destructor   ! Nullify all pointers

    end type functional_cache_handler_t
    !***************************************************************************************************



contains



    !>  This procedure begins the functional computation either on primary or auxiliary geometries
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/23/2018
    !!
    !!  restructured
    !!
    !!  @author matteo ugolotti
    !!  @date   3/7/2019
    !!
    !!  @param[inout]      worker           chidg worker
    !!  @param[inout]      data             chidg data
    !!  @param[inout]      fcl_cache        functional storage
    !!  @param[in]         ifunc            functional ID 
    !!  @param[in]         geom_feature     geometry on which the functional needs to 
    !!                                      be compute: reference/auxiliary
    !!  @param[in]         differentiate    flag for functional derivative computation
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine update(self,worker,fcl_cache,data,ifunc,geom_feature,differentiate)
        class(functional_cache_handler_t),  intent(inout)       :: self
        type(chidg_worker_t),               intent(inout)       :: worker
        type(functional_cache_t),           intent(inout)       :: fcl_cache
        type(chidg_data_t),                 intent(inout)       :: data
        integer(ik),                        intent(in)          :: ifunc
        character(*),                       intent(in)          :: geom_feature
        integer(ik),                        intent(in)          :: differentiate

        ! Initialize functional_cache_handler 
        call self%init(data,fcl_cache,ifunc,geom_feature,differentiate)
        
        ! Compute functional
        call self%compute(worker)

        ! Store functional
        call self%store()

        ! Destructor
        call self%destructor()

    end subroutine update
    !***************************************************************************************************






    !>  Here the handler is initialize so that it knows functional and geometry we are currently working on.
    !!  Pointers are created to facilitate the access to the main objects in the handler, that are:
    !!      - ChiDG data (for storage)
    !!      - Functional cache (for integral storage during functional computation)
    !!      - Functional (containing the procedures to actually compute the specific functional)
    !!
    !!  Furthermore, type of differentiation, geometry of current computation (auxiliary/reference) and 
    !!  cache component are saved to be easily accessible during the computation.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  @param[inout]      worker           chidg worker
    !!  @param[inout]      data             chidg data
    !!  @param[inout]      fcl_cache        functional storage
    !!  @param[in]         ifunc            functional ID 
    !!  @param[in]         geom_feature     geometry on which the functional needs to 
    !!                                      be compute: reference/auxiliary
    !!  @param[in]         differentiate    flag for functional derivative computation
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init(self,data,fcl_cache,ifunc,geom_feature,differentiate)
        class(functional_cache_handler_t),  intent(inout)           :: self
        type(chidg_data_t),                 intent(inout),  target  :: data
        type(functional_cache_t),           intent(inout),  target  :: fcl_cache
        integer(ik),                        intent(in)              :: ifunc
        character(*),                       intent(in)              :: geom_feature
        integer(ik),                        intent(in)              :: differentiate


        associate (func => data%functional_group%fcl_entities(ifunc)%func )

            ! Assign active geometry (reference/auxiliary),data, cache, initialize active_cache
            select case (trim(geom_feature))
                case("auxiliary")

                    ! Currently working on auxiliary geometry
                    self%geom  =  trim(geom_feature)
                    
                    ! Create pointers
                    self%cache => fcl_cache
                    self%data  => data

                    ! Initialize fuctional auxiliary cache
                    call self%cache%init(self%geom,data%mesh,func%auxiliary_geom,func%get_int_type(),differentiate)


                case("reference")
                    
                    ! Currently working on reference geometry
                    self%geom  =  trim(geom_feature)
                    
                    ! Create pointers
                    self%cache => fcl_cache
                    self%data  => data
                    
                    ! Initialize fuctional auxiliary cache
                    call self%cache%init(self%geom,data%mesh,func%reference_geom,func%get_int_type(),differentiate)
                    
                case default
                    call chidg_signal_two(FATAL,"functional_cache_handler_t%init: wrong geometry string provided to the functional cache handler", trim(geom_feature),".") 
            end select


            ! Define type of integral 
            ! NOTE: this can be modified to handle different type of integrals on auxiliary and reference geometries
            select case (func%get_int_type())
                case('VOLUME INTEGRAL')
                    self%component = 'element'
                case('FACE INTEGRAL')
                    self%component = 'interior faces'
                case default
                    call chidg_signal_two(FATAL,"functional_cache_handler_t%init: wrong type of integral provided to the worker",func%get_int_type(),".") 
            end select

            
            ! Save locally type of differentiation
            self%dtype = differentiate

            ! Functiona ID
            self%ifunc = ifunc

        end associate


    end subroutine init
    !***************************************************************************************************






    !>  TASKS:
    !!      1) Compute the functional at each element/face by looping through the entities in the cache and 
    !!         store them in the functional cache 
    !!      2) Carry out the communication among the processor of the integrals values (functional can be made
    !!         off multiple integrals)
    !!      3) Finalize the functional on either the auxiliary or reference geometry (operations done using overall
    !!         integrals.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/25/2018
    !!
    !!  Added dX linearization capability
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/5/2018
    !!
    !!  restructured
    !!
    !!  @author matteo ugolotti
    !!  @date   3/7/2019
    !!
    !!  @param[inout]      worker           chidg worker
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(functional_cache_handler_t),  intent(inout)       :: self
        type(chidg_worker_t),               intent(inout)       :: worker
   
        type(element_info_t)        :: elem_info
        type(cache_handler_t)       :: worker_cache_handler
        integer(ik)                 :: idepend, ientity, idiff, iface, ierr, idomain_l, ielement_l
        logical                     :: compute_auxiliary, compute_reference

        associate (func => self%data%functional_group%fcl_entities(self%ifunc)%func )

            ! Select the functional computation 
            compute_auxiliary = ( self%geom == 'auxiliary' )
            compute_reference = ( self%geom == 'reference' )

            ! Loop through all the entities initialized (elements/faces) 
            do ientity = 1,self%cache%nentities(self%geom)

                ! Update the worker
                if (compute_auxiliary) then
                    idomain_l  = self%cache%aux_cache%idomain_g%at(ientity)
                    ielement_l = self%cache%aux_cache%ielement_l%at(ientity)
                    iface      = self%cache%aux_cache%iface%at(ientity)
                else
                    idomain_l  = self%cache%ref_cache%idomain_l%at(ientity)
                    ielement_l = self%cache%ref_cache%ielement_l%at(ientity)
                    iface      = self%cache%ref_cache%iface%at(ientity)
                end if

                elem_info = worker%mesh%get_element_info(idomain_l,ielement_l)
                call worker%set_element(elem_info)


                ! Compute differential interpolator for mesh sensititivites
                call worker%mesh%compute_derivatives_dx(elem_info,self%dtype)


                ! Update worker cache
                call worker_cache_handler%update(worker,self%data%eqnset,       &
                                                 self%data%bc_state_group,      &
                                                 components = self%component,   &
                                                 face = iface,                  & 
                                                 differentiate = self%dtype,    &
                                                 lift = .false.)
                ! Set up worker
                idepend = 1
                idiff   = DIAG
                

                if (self%component == 'interior faces') then
                    worker%function_info%seed = face_compute_seed(worker%mesh,          &
                                                                  elem_info%idomain_l,  &
                                                                  elem_info%ielement_l, &
                                                                  iface,idepend,idiff,worker%itime)
                    call worker%set_face(iface)
                
                else if (self%component == 'element') then
                    worker%function_info%seed = element_compute_seed(worker%mesh,          & 
                                                                     elem_info%idomain_l,  &
                                                                     elem_info%ielement_l, &
                                                                     idepend,idiff,worker%itime)
                end if 
                

                worker%function_info%idepend = idepend
                worker%function_info%idiff   = idiff
                worker%function_info%dtype   = self%dtype

                ! Compute functional on the entity (element/face)
                if (compute_auxiliary) call func%compute_auxiliary(worker,self%cache)
                if (compute_reference) call func%compute_functional(worker,self%cache)

                ! Release differential interpolators allocated memory
                call worker%mesh%release_derivatives_dx(elem_info,self%dtype)

            end do !ientity

            ! Communicate integrals through processors
            call self%cache%comm(self%geom)
            
            ! Finalize functional on auxiliary/reference geometry
            if (compute_auxiliary) call func%finalize_auxiliary(worker,self%cache)
            if (compute_reference) call func%finalize_functional(worker,self%cache)

        end associate
    
    end subroutine compute
    !***************************************************************************************************







    !>  This procedure stores the computed functional based on the type of the functional itself:
    !!      - functional derivatives (derivatives of the functional wrt dQ or dX are stored)
    !!      - functional             (only the real value of the functional is of interest and stored)
    !!      - auxiliary              (store the auxiliary value needed for the upcoming actuall functional)
    !!
    !!  - Functional derivatives are stored in the adjoint%Jq vector
    !!  - Functional real values are stored in the sdata%functional%func
    !!  - Auxiliary value is stored in the functional itself
    !! 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/25/2018
    !!
    !!  Added dX linearization capability
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/5/2018
    !!
    !!  restructured
    !!
    !!  @author matteo ugolotti
    !!  @date   3/7/2019
    !!
    !!
    !!  TODO: maybe this can be improved
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine store(self)
        class(functional_cache_handler_t),  intent(inout)       :: self
        
        logical                 :: proceed_storing, store_functional,          &
                                   store_derivatives_dQ, store_derivatives_dX,  &
                                   store_auxiliary, compatibility
        type(integral_cache_t)  :: global_integral
        type(AD_D)              :: finalized_integral

        real(rk)                :: time
        Integer(ik)             :: idomain_l, ielement_l, ientity, i_int, step


        ! Prepare the type of store we need to execute
        if ( (self%geom == 'reference') .and. (self%dtype == dQ_DIFF) ) then
            store_derivatives_dQ = .true.
            store_derivatives_dX = .false.
            store_functional     = .true.
        end if
        if ( (self%geom == 'reference') .and. (self%dtype == dX_DIFF) ) then
            store_derivatives_dQ = .false.
            store_derivatives_dX = .true.
            store_functional     = .false.
        end if
        if ( (self%geom == 'reference') .and. (self%dtype == NO_DIFF) ) then
            store_derivatives_dQ = .false.
            store_derivatives_dX = .false.
            store_functional     = .true.
        end if
        if ( self%geom == 'auxiliary' ) then
            store_derivatives_dQ = .false.
            store_derivatives_dX = .false.
            store_functional     = .false.
        end if
       

        associate (func => self%data%functional_group%fcl_entities(self%ifunc)%func )

            ! Update local integrals with overall value coming from all ranks and store
            ! dQ derivatives in the sdata%adjoint%Jq
            if (store_derivatives_dQ) then
                self%data%sdata%adjoint%Jq(self%ifunc) = func%store_deriv(self%cache)
                call self%data%sdata%adjoint%Jq(self%ifunc)%assemble()
            end if


            ! Update local integrals with overall value coming from all ranks and store
            ! dX derivatives in the Jx(ifunc) vector
            if (store_derivatives_dX) then
                self%data%sdata%adjointx%Jx(self%ifunc) = func%store_deriv(self%cache)
                call self%data%sdata%adjoint%Jq(self%ifunc)%assemble()
            end if


            ! Store real value of the functional
            if (store_functional) then
                
                ! Store value
                call self%data%sdata%functional%func(self%ifunc)%push_back(func%store_value(self%cache))
                
                ! Store time and step
                ! Istep and itime are the same for all the functional
                ! We do not want to keep pushing times and steps, we do it once for the first
                ! functional.
                step = self%data%time_manager%istep
                if (self%ifunc == 1) call self%data%sdata%functional%step%push_back(step)
                time = self%data%time_manager%t
                if (self%ifunc == 1) call self%data%sdata%functional%time%push_back(time)

            end if

        end associate


    end subroutine store
    !***************************************************************************************************



    !>  Nullify all the pointers before moving to the next geometry or functional
    !!
    !!  @author matteo ugolotti
    !!  @date   3/8/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        class(functional_cache_handler_t),  intent(inout)       :: self

        ! Deassociate pointers in functional cache handler
        if (associated(self%cache)) nullify(self%cache)
        if (associated(self%data))  nullify(self%data)

    end subroutine destructor
    !***************************************************************************************************


end module type_functional_cache_handler
