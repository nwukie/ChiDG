module type_equation_set
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: INTERIOR, CHIMERA, DIAG, &
                                              BOUNDARY_ADVECTIVE_FLUX, &
                                              BOUNDARY_DIFFUSIVE_FLUX, &
                                              VOLUME_ADVECTIVE_FLUX, VOLUME_DIFFUSIVE_FLUX, &
                                              BC_FLUX, XI_MIN, XI_MAX, ETA_MIN, ETA_MAX,    &
                                              ZETA_MIN, ZETA_MAX, BOUNDARY
    use mod_operators,                  only: operator_factory
    use mod_models,                     only: model_factory
    use mod_DNAD_tools,                 only: element_compute_seed, face_compute_seed

    use type_operator,                  only: operator_t
    use type_operator_wrapper,          only: operator_wrapper_t
    use type_model,                     only: model_t
    use type_model_wrapper,             only: model_wrapper_t
    use type_properties,                only: properties_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_mesh,                  only: mesh_t
    use type_bc,                        only: bc_t
    use type_solverdata,                only: solverdata_t
    use type_element_info,              only: element_info_t
    use type_face_info,                 only: face_info_t
    use type_function_info,             only: function_info_t
    use type_pseudo_timestep,           only: pseudo_timestep_t, default_pseudo_timestep_t
    implicit none
    private



    !>  Abstract equation-set type. Can be extended to implement a concrete equation set.
    !!      - Contains name and number of equations.
    !!      - Contains properties type with equations and definitions
    !!      - Contains arrays of flux components
    !!
    !!  When a new equation set is defined. It should extend from this abstract type. 
    !!  It must then implement the 'init' function where equations/variables can be added, 
    !!  and fluxes can be added to the definition of the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Reorganized to contain operators instead of fluxes
    !!
    !----------------------------------------------------------------------------------------
    type, public :: equation_set_t

        ! Name
        character(:),               allocatable :: name 

        ! ID
        integer(ik)                             :: eqn_ID

        ! Properties/Fields
        type(properties_t)                      :: prop

        ! Operators
        type(operator_wrapper_t),   allocatable :: boundary_advective_operator(:)
        type(operator_wrapper_t),   allocatable :: boundary_diffusive_operator(:)
        type(operator_wrapper_t),   allocatable :: volume_advective_operator(:)
        type(operator_wrapper_t),   allocatable :: volume_diffusive_operator(:) 
        type(operator_wrapper_t),   allocatable :: bc_operator(:)

        ! Models
        type(model_wrapper_t),      allocatable :: models(:)

        ! Pseudo time-step calculator
        class(pseudo_timestep_t),   allocatable :: pseudo_timestep

        ! Data for the flux and source functions. Ex how many. This gets passed to a container 
        ! in sdata that keeps track of whether these have been executed or not.
        type(equationset_function_data_t)       :: function_data

    contains

        procedure   :: set_name                             
        procedure   :: get_name                             

        procedure   :: add_operator 
        procedure   :: add_model
        procedure   :: add_pseudo_timestep                  


        procedure   :: compute_boundary_advective_operators 
        procedure   :: compute_boundary_diffusive_operators 
        procedure   :: compute_volume_advective_operators   
        procedure   :: compute_volume_diffusive_operators   
        procedure   :: compute_bc_operators                 

        procedure   :: compute_pseudo_timestep

        !procedure   :: get_boundary_ndependent_elements    
        !procedure   :: get_volume_ndependent_elements       
        procedure   :: get_face_ncompute
        procedure   :: get_element_ncompute

        procedure   :: nmodels

    end type equation_set_t
    !**************************************************************************************







    !> Interface definitions
    abstract interface
        subroutine self_interface(self)
            import equation_set_t
            class(equation_set_t), intent(inout) :: self
        end subroutine
    end interface


contains





    !> Set name of the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  name_string     Character string indicating the name of the equation set
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_name(self,ename)
        class(equation_set_t),  intent(inout)   :: self
        character(*),           intent(in)      :: ename

        self%name = ename

    end subroutine set_name
    !**************************************************************************************




    !>  Return the name of the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_name(self) result(ename)
        class(equation_set_t),   intent(in)   :: self

        character(:),   allocatable :: ename

        ename = self%name

    end function get_name
    !**************************************************************************************








    !>  Add a spatial operator to the equation set.
    !!
    !!  This accepts a string that is the name of an operator to add. This function then 
    !!  searched the operator register for an operator that matches the incoming string. 
    !!  If one is found, it is created and added to the equation set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine add_operator(self,string)
        class(equation_set_t),  intent(inout)   :: self
        character(*),           intent(in)      :: string

        class(operator_t),          allocatable :: new_operator
        class(operator_wrapper_t),  allocatable :: temp(:)
        integer(ik)                             :: ierr, iflux, operator_type, ifield, imodel


        !
        ! Create new operator
        !
        allocate(new_operator, source=operator_factory%produce(string), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Get operator type
        !
        operator_type = new_operator%get_operator_type()

        
        !
        ! Add to correct operator array
        !
        if (operator_type == BOUNDARY_ADVECTIVE_FLUX) then



            ! Allocate temporary flux array with one additional slot
            if (allocated(self%boundary_advective_operator)) then

                ! Allocate
                allocate(temp(size(self%boundary_advective_operator) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do iflux = 1,size(self%boundary_advective_operator)
                    allocate(temp(iflux)%op,source=self%boundary_advective_operator(iflux)%op, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else
                ! Allocate new slot
                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if




        else if (operator_type == BOUNDARY_DIFFUSIVE_FLUX) then

            ! Allocate temporary flux array with one additional slot
            if (allocated(self%boundary_diffusive_operator)) then

                ! Allocate
                allocate(temp(size(self%boundary_diffusive_operator) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do iflux = 1,size(self%boundary_diffusive_operator)
                    allocate(temp(iflux)%op,source=self%boundary_diffusive_operator(iflux)%op, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else
                ! Allocate new slot
                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if




        else if (operator_type == VOLUME_ADVECTIVE_FLUX) then



            ! Allocate temporary flux array with one additional slot
            if (allocated(self%volume_advective_operator)) then

                ! Allocate
                allocate(temp(size(self%volume_advective_operator) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do iflux = 1,size(self%volume_advective_operator)
                    allocate(temp(iflux)%op,source=self%volume_advective_operator(iflux)%op, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else

                ! Allocate new slot
                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if





        else if (operator_type == VOLUME_DIFFUSIVE_FLUX) then



            ! Allocate temporary flux array with one additional slot
            if (allocated(self%volume_diffusive_operator)) then

                ! Allocate
                allocate(temp(size(self%volume_diffusive_operator) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do iflux = 1,size(self%volume_diffusive_operator)
                    allocate(temp(iflux)%op,source=self%volume_diffusive_operator(iflux)%op, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else

                ! Allocate new slot
                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if








        else if (operator_type == BC_FLUX) then


            ! Allocate temporary flux array with one additional slot
            if (allocated(self%bc_operator)) then

                ! Allocate
                allocate(temp(size(self%bc_operator) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do iflux = 1,size(self%bc_operator)
                    allocate(temp(iflux)%op,source=self%bc_operator(iflux)%op, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else

                ! Allocate new slot
                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if




        else
            call chidg_signal_one(FATAL,"equation_set%add_operator: 'Operator type was not valid'", operator_type)
        end if





        !
        ! Assign new operator
        !
        allocate(temp(size(temp))%op, source=new_operator, stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        ! Copy extended temp array to equation set storage
        !
        if (operator_type == BOUNDARY_ADVECTIVE_FLUX) then
            self%boundary_advective_operator = temp
            self%function_data%nboundary_advective_flux = size(self%boundary_advective_operator)


        else if (operator_type == BOUNDARY_DIFFUSIVE_FLUX) then
            self%boundary_diffusive_operator = temp
            self%function_data%nboundary_diffusive_flux = size(self%boundary_diffusive_operator)

        else if (operator_type == VOLUME_ADVECTIVE_FLUX) then
            self%volume_advective_operator = temp
            self%function_data%nvolume_advective_flux = size(self%volume_advective_operator)

        else if (operator_type == VOLUME_DIFFUSIVE_FLUX) then
            self%volume_diffusive_operator = temp
            self%function_data%nvolume_diffusive_flux = size(self%volume_diffusive_operator)

        else if (operator_type == BC_FLUX) then
            self%bc_operator = temp

        else
            call chidg_signal_one(FATAL,"equation_set%add_operator: 'Operator type was not valid'", operator_type)
        end if




        ! Turn on primary fields from the new operator
        do ifield = 1,new_operator%nprimary_fields()
            call self%prop%add_primary_field(new_operator%get_primary_field(ifield))
        end do

        ! Turn on auxiliary fields from the new operator
        do ifield = 1,new_operator%nauxiliary_fields()
            call self%prop%add_auxiliary_field(new_operator%get_auxiliary_field(ifield))
        end do

        ! Turn on models from the new operator
        do imodel = 1,new_operator%nmodels()
            call self%add_model(new_operator%get_model(imodel))
        end do

    end subroutine add_operator
    !****************************************************************************************









    !>  Add a model to the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine add_model(self,string)
        class(equation_set_t),  intent(inout)   :: self
        character(*),           intent(in)      :: string
        
        integer(ik)                         :: imodel, ierr, ifield
        logical                             :: already_added
        class(model_t),         allocatable :: new_model
        type(model_wrapper_t),  allocatable :: temp(:)


        !
        ! Check that the model wasn't already added.
        !
        already_added = .false.
        do imodel = 1,self%nmodels()
            already_added = (trim(string) == self%models(imodel)%model%get_name())
            if (already_added) exit
        end do



        !
        ! Add the model if it wasn't added already by another operator.
        !
        if (.not. already_added) then

            allocate(new_model, source=model_factory%produce(string), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Extend storage
            !
            if (allocated(self%models)) then

                allocate(temp(size(self%models) + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current flux components to temp array
                do imodel = 1,size(self%models)
                    allocate(temp(imodel)%model,source=self%models(imodel)%model, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end do

            else

                allocate(temp(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if



            !
            ! Allocate new model to end of extended array
            !
            allocate(temp(size(temp))%model, source=new_model, stat=ierr)
            if (ierr /= 0) call AllocationError


            ! Turn on model fields from the new model
            do ifield = 1,new_model%nmodel_fields()
                call self%prop%add_model_field(new_model%get_model_field(ifield))
            end do


            !
            ! Move extended temp allocation to data type
            !
            call move_alloc(from=temp, to=self%models)

        end if


    end subroutine add_model
    !***************************************************************************************








    !>  Set a pseudo time-step calculator.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine add_pseudo_timestep(self,calculator)
        class(equation_set_t),      intent(inout)   :: self
        class(pseudo_timestep_t),   intent(in)      :: calculator

        integer(ik) :: ierr

        if (allocated(self%pseudo_timestep)) deallocate(self%pseudo_timestep)

        allocate(self%pseudo_timestep, source=calculator, stat=ierr)
        if (ierr /= 0) call AllocationError 

    end subroutine add_pseudo_timestep
    !****************************************************************************************








    !>  Loops through the attached boundary advective operator functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine compute_boundary_advective_operators(self,worker,differentiate)
        class(equation_set_t),   intent(inout)  :: self
        type(chidg_worker_t),    intent(inout)  :: worker
        logical,                 intent(in)     :: differentiate

        integer(ik),    allocatable :: compute_pattern(:)
        integer(ik)                 :: nfcn, ifcn, icompute, ncompute, ipattern, idiff
        logical                     :: interior_face, chimera_face, compute, &
                                       compute_function, differentiate_function

        associate( mesh  => worker%mesh,                    &
                   idom  => worker%element_info%idomain_l,  &
                   ielem => worker%element_info%ielement_l, &
                   iface => worker%iface,                   &
                   prop  => self%prop )



        !
        ! Only call 'Boundary Advective' operators for face type:
        !   - INTERIOR
        !   - CHIMERA
        !
        interior_face = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )
        chimera_face  = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
        compute       = (interior_face .or. chimera_face)


        if (compute) then

            !
            ! Determine pattern to compute functions. Depends on if we are differentiating 
            ! or not. These will be used to set idiff, indicating the differentiation
            ! direction.
            !
            if (differentiate) then
                ! compute function, wrt internal, external(in direction iface) states
                compute_pattern = [DIAG,iface]
            else
                ! compute function, but do not differentiate
                compute_pattern = [0]
            end if


            !
            ! Execute compute pattern
            !
            do ipattern = 1,size(compute_pattern)

                !
                ! Set differentiation indicator
                !
                idiff = compute_pattern(ipattern)

                !
                ! Determine the number of times to execute the functions.
                ! Depends on if we are differentiating or not.
                !
                ncompute = self%get_face_ncompute(mesh,worker%face_info(),idiff)


                if (allocated(self%boundary_advective_operator)) then
                    nfcn = size(self%boundary_advective_operator)
                    do ifcn = 1,nfcn

                        worker%function_info%type   = BOUNDARY_ADVECTIVE_FLUX
                        worker%function_info%ifcn   = ifcn
                        worker%function_info%idiff  = idiff

                        compute_function       = worker%solverdata%function_status%compute_function(   worker%face_info(), worker%function_info )
                        differentiate_function = worker%solverdata%function_status%linearize_function( worker%face_info(), worker%function_info )
                        

                        if ( compute_function .or. differentiate_function ) then
                            !
                            ! If we are differentiating, compute boundary flux once for each 
                            ! contributing element:
                            !   - For conforming faces ncompute == 1. 
                            !   - For Chimera faces ncompute is potentially > 1
                            !
                            do icompute = 1,ncompute
                                worker%function_info%seed    = face_compute_seed(mesh,idom,ielem,iface,icompute,idiff)
                                worker%function_info%idepend = icompute

                                call self%boundary_advective_operator(ifcn)%op%compute(worker,prop)

                            end do !icompute
                        end if


                    end do ! ifcn

                end if ! allocated boundary_advective_flux

            end do !ipattern

        end if !compute



        end associate

    end subroutine compute_boundary_advective_operators
    !***************************************************************************************









    !>  Loops through the attached boundary diffusive operator functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_boundary_diffusive_operators(self,worker,differentiate)
        class(equation_set_t),  intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        logical,                intent(in)      :: differentiate

        integer(ik),    allocatable :: compute_pattern(:)
        integer(ik)                 :: nfcn, ifcn, ipattern, icompute, ncompute, idiff, &
                                       ChiID, eqn_m, eqn_p
        logical                     :: compute_function, linearize_function,            &
                                       interior_face, chimera_face, compute, skip

        associate( mesh  => worker%mesh,                    &
                   idom  => worker%element_info%idomain_l,  &
                   ielem => worker%element_info%ielement_l, &
                   iface => worker%iface,                   &
                   prop  => self%prop)



        !
        ! Only call 'Boundary Advective' operators for face type:
        !   - INTERIOR
        !   - CHIMERA
        !
        interior_face = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )
        chimera_face  = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA  )
        compute       = (interior_face .or. chimera_face)


        !
        ! Determine if multi-fidelity interface. Check that equation sets match:
        !   - interior faces: should have the same equation set. No need to skip
        !   - chimera faces: could be different. If they are, we want to skip 
        !
        skip = .false.
        if (interior_face) then
            skip = .false.

        else if (chimera_face) then

            ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
            eqn_m = self%eqn_ID
            eqn_p = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_eqn_ID%at(1)
            skip = (eqn_m /= eqn_p)

        end if




        if (compute .and. (.not. skip)) then

            !
            ! Determine pattern to compute functions. Depends on if we are differentiating 
            ! or not. These will be used to set idiff, indicating the differentiation
            ! direction.
            !
            if (differentiate) then
                ! compute function, wrt internal, external(in direction iface) states
                compute_pattern = [DIAG,iface]
            else
                ! compute function, but do not differentiate
                compute_pattern = [0]
            end if


            !
            ! Execute compute pattern
            !
            do ipattern = 1,size(compute_pattern)

            
                !
                ! Set differentiation indicator
                !
                idiff = compute_pattern(ipattern)


                !
                ! Get number of elements we are linearizing with respect to
                !
                ncompute = self%get_face_ncompute(mesh,worker%face_info(),idiff)


                if (allocated(self%boundary_diffusive_operator)) then
                    nfcn = size(self%boundary_diffusive_operator)
                    do ifcn = 1,nfcn

                        worker%function_info%type   = BOUNDARY_DIFFUSIVE_FLUX
                        worker%function_info%ifcn   = ifcn
                        worker%function_info%idiff  = idiff

                        compute_function     = worker%solverdata%function_status%compute_function(   worker%face_info(), worker%function_info )
                        linearize_function   = worker%solverdata%function_status%linearize_function( worker%face_info(), worker%function_info )
                        

                        if ( compute_function .or. linearize_function ) then
                            !
                            ! If we are differentiating, compute boundary flux once for each 
                            ! contributing element:
                            !   - For conforming faces ncompute == 1. 
                            !   - For Chimera faces ncompute is potentially > 1
                            !
                            do icompute = 1,ncompute
                                worker%function_info%seed    = face_compute_seed(mesh,idom,ielem,iface,icompute,idiff)
                                worker%function_info%idepend = icompute

                                call self%boundary_diffusive_operator(ifcn)%op%compute(worker,prop)

                            end do
                        end if


                    end do ! ifcn

                end if ! boundary_diffusive_flux loop

            end do !ipattern

        end if ! compute



        end associate

    end subroutine compute_boundary_diffusive_operators
    !***************************************************************************************













    !>  Loops through the attached volume advective operator functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_volume_advective_operators(self,worker,differentiate)
        class(equation_set_t),      intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        logical,                    intent(in)      :: differentiate

        integer(ik),    allocatable :: compute_pattern(:)
        integer(ik)                 :: nfcn, ifcn, icompute, ncompute, idiff, ipattern
        type(function_info_t)       :: function_info
        logical                     :: compute

        associate( mesh      => worker%mesh,                    &
                   elem_info => worker%element_info,            &
                   idom      => worker%element_info%idomain_l,  &
                   ielem     => worker%element_info%ielement_l, &
                   prop      => self%prop)


        !
        ! Determine pattern to compute functions. Depends on if we are differentiating 
        ! or not. These will be used to set idiff, indicating the differentiation
        ! direction.
        !
        if (differentiate) then
            ! compute function, wrt internal state
            compute_pattern = [DIAG]
        else
            ! compute function, but do not differentiate
            compute_pattern = [0]
        end if



        !
        ! Execute compute pattern
        !
        do ipattern = 1,size(compute_pattern)

            !
            ! Set differentiation indicator
            !
            idiff = compute_pattern(ipattern)

            !
            ! Volume advective flux only depends on the interior element
            !
            icompute = 1

            if (allocated(self%volume_advective_operator)) then
                nfcn = size(self%volume_advective_operator)
                do ifcn = 1,nfcn

                    worker%function_info%type    = VOLUME_ADVECTIVE_FLUX
                    worker%function_info%ifcn    = ifcn
                    worker%function_info%idiff   = idiff
                    worker%function_info%idepend = icompute
                    worker%function_info%seed    = element_compute_seed(mesh,idom,ielem,icompute,idiff)

                    call self%volume_advective_operator(ifcn)%op%compute(worker,prop)

                end do ! ifcn
            end if

        end do ! ipattern

        end associate

    end subroutine compute_volume_advective_operators
    !***************************************************************************************










    !>  Loops through the attached volume diffusive operator functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_volume_diffusive_operators(self,worker,differentiate)
        class(equation_set_t),      intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        logical,                    intent(in)      :: differentiate

        integer(ik),    allocatable :: compute_pattern(:)
        integer(ik)                 :: nfcn, ifcn, icompute, ncompute, ipattern, idiff
        type(function_info_t)       :: function_info
        logical                     :: diff_none, diff_interior, diff_exterior
        logical                     :: linearize_me, compute_function

        associate( mesh      => worker%mesh,                    &
                   elem_info => worker%element_info,            &
                   idom      => worker%element_info%idomain_l,  &
                   ielem     => worker%element_info%ielement_l, &
                   prop      => self%prop)


        !
        ! Determine pattern to compute functions. Depends on if we are differentiating 
        ! or not. These will be used to set idiff, indicating the differentiation
        ! direction.
        !
        if (differentiate) then
            ! compute function, wrt (all exterior)/interior states
            compute_pattern = [1,2,3,4,5,6,DIAG]
        else
            ! compute function, but do not differentiate
            compute_pattern = [0]
        end if

        !
        ! Execute compute pattern
        !
        do ipattern = 1,size(compute_pattern)

        
            !
            ! Set differentiation indicator
            !
            idiff = compute_pattern(ipattern)

            diff_none = (idiff == 0)
            diff_interior = (idiff == DIAG)
            diff_exterior = ( (idiff == 1) .or. (idiff == 2) .or. &
                              (idiff == 3) .or. (idiff == 4) .or. &
                              (idiff == 5) .or. (idiff == 6) )



            !
            ! Get number of elements we are linearizing with respect to
            !
            ncompute = self%get_element_ncompute(mesh,elem_info,idiff)



            if (diff_none) then
                compute_function = .true.
            else if (diff_interior) then
                compute_function = .true.
            else if (diff_exterior) then
                compute_function = ( (mesh%domain(elem_info%idomain_l)%faces(elem_info%ielement_l,idiff)%ftype == INTERIOR) .or. &
                                     (mesh%domain(elem_info%idomain_l)%faces(elem_info%ielement_l,idiff)%ftype == CHIMERA) )
            end if




                                

            if (compute_function) then

                if (allocated(self%volume_diffusive_operator)) then
                    nfcn = size(self%volume_diffusive_operator)
                    do ifcn = 1,nfcn
        
                        !
                        ! Compute boundary flux once for each donor. 
                        !   - For interior faces ndepend == 1. 
                        !   - For Chimera faces ndepend is potentially > 1.
                        !
                        do icompute = 1,ncompute
        
                            worker%function_info%type    = VOLUME_DIFFUSIVE_FLUX
                            worker%function_info%ifcn    = ifcn
                            worker%function_info%idiff   = idiff
                            worker%function_info%idepend = icompute
                            worker%function_info%seed    = element_compute_seed(mesh,idom,ielem,icompute,idiff)
        
                            call self%volume_diffusive_operator(ifcn)%op%compute(worker,prop)
        
                        end do
        
                    end do ! ifcn
                end if

            end if

        end do !ipattern

        end associate

    end subroutine compute_volume_diffusive_operators
    !***************************************************************************************













    !>  Loops through the attached boundary advective operator functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_bc_operators(self,worker,bc,differentiate)
        class(equation_set_t),      intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(bc_t),                 intent(inout)   :: bc(:)
        logical,                    intent(in)      :: differentiate

        integer(ik),    allocatable :: compute_pattern(:)
        integer(ik)                 :: nfcn, ifcn, icompute, ncompute, bc_ID, patch_ID, &
                                       face_ID, idiff, ipattern
        logical                     :: boundary_face, compute

        associate( mesh  => worker%mesh,                    &
                   idom  => worker%element_info%idomain_l,  &
                   ielem => worker%element_info%ielement_l, &
                   iface => worker%iface,                   &
                   prop  => self%prop )

        !
        ! Only call the following routines for interior faces -- ftype == 0
        ! Furthermore, only call the routines if we are computing derivatives for the 
        ! neighbor of iface or for the current element(DIAG). This saves a lot of 
        ! unnecessary compute_boundary calls.
        !
        boundary_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == BOUNDARY)
        compute       = (boundary_face)

        if (compute) then

            !
            ! Determine pattern to compute functions. Depends on if we are differentiating 
            ! or not. These will be used to set idiff, indicating the differentiation
            ! direction.
            !
            if (differentiate) then
                ! compute function, wrt (all exterior)/interior states
                compute_pattern = [DIAG]
            else
                ! compute function, but do not differentiate
                compute_pattern = [0]
            end if

            
            do ipattern = 1,size(compute_pattern)

                !
                ! Set differentiation indicator
                !
                idiff = compute_pattern(ipattern)

                !
                ! Get number of elements we are linearizing with respect to
                !
                ncompute = self%get_face_ncompute(mesh,worker%face_info(),idiff)


                !
                ! Get index of boundary condition, patch, patch face. 
                !
                bc_ID    = mesh%domain(idom)%faces(ielem,iface)%bc_ID
                patch_ID = mesh%domain(idom)%faces(ielem,iface)%patch_ID
                face_ID  = mesh%domain(idom)%faces(ielem,iface)%face_ID


                if (allocated(self%bc_operator)) then
                    nfcn = size(self%bc_operator)
                    do ifcn = 1,nfcn

                        worker%function_info%type   = BC_FLUX
                        worker%function_info%ifcn   = ifcn
                        worker%function_info%idiff  = idiff

                        !
                        ! Compute boundary flux once for each donor. 
                        !   - For interior faces ndepend == 1. 
                        !   - For Chimera faces ndepend is potentially > 1.
                        !   - For BC faces ndepend is potentially > 1.
                        !
                        do icompute = 1,ncompute


                            !
                            ! Get coupled element to linearize against.
                            !
                            worker%function_info%seed%idomain_g  = bc(bc_ID)%bc_patch(patch_ID)%idomain_g_coupled(face_ID)%at(icompute)
                            worker%function_info%seed%idomain_l  = bc(bc_ID)%bc_patch(patch_ID)%idomain_l_coupled(face_ID)%at(icompute)
                            worker%function_info%seed%ielement_g = bc(bc_ID)%bc_patch(patch_ID)%ielement_g_coupled(face_ID)%at(icompute)
                            worker%function_info%seed%ielement_l = bc(bc_ID)%bc_patch(patch_ID)%ielement_l_coupled(face_ID)%at(icompute)
                            worker%function_info%seed%iproc      = bc(bc_ID)%bc_patch(patch_ID)%proc_coupled(face_ID)%at(icompute)
                            worker%function_info%idepend         = icompute

                            call self%bc_operator(ifcn)%op%compute(worker,prop)

                        end do


                    end do ! ifcn

                end if ! boundary_advective_flux loop

            end do ! ipattern

        end if !compute

        end associate

    end subroutine compute_bc_operators
    !**************************************************************************************











    !>  Return the number of times to compute a function on a face. Depends on if 
    !!  we are differentiating or not, indicated by the value of 'idiff'.
    !!
    !!  Often, this may be just 1. Because a flux differentiated wrt its owner element 
    !!  just has 1 dependent element; itself. A flux differentiated wrt its neighbor 
    !!  element has just 1 dependent element; it's neighbor. If however, a face is a 
    !!  Chimera face, then when differentiated with respect to its donors, ncompute will 
    !!  be equal to the number of donor elements that are providing information to 
    !!  that face.
    !!
    !!  Additionally, if we are not differentiating(idiff == 0) we at least want to 
    !!  compute the function once so ncompute = 1.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/15/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_face_ncompute(self,mesh,face_info,idiff) result(ncompute)
        class(equation_set_t),  intent(in)   :: self
        type(mesh_t),       intent(in)   :: mesh
        type(face_info_t),      intent(in)   :: face_info
        integer(ik),            intent(in)   :: idiff

        integer(ik) :: ChiID, ncompute
        logical     :: compute_wrt_none, compute_wrt_interior, compute_wrt_exterior, &
                       chimera_face, bc_face

        associate( idom  => face_info%idomain_l,  &
                   ielem => face_info%ielement_l, &
                   iface => face_info%iface )

        
        compute_wrt_none     = (idiff == 0)
        compute_wrt_interior = (idiff == DIAG)
        compute_wrt_exterior = (idiff == XI_MIN   .or. idiff == XI_MAX  .or. &
                                idiff == ETA_MIN  .or. idiff == ETA_MAX .or. &
                                idiff == ZETA_MIN .or. idiff == ZETA_MAX)

        if (compute_wrt_none) then
            ncompute = 1

        else if (compute_wrt_interior) then
            ncompute = 1

        else if (compute_wrt_exterior) then

            chimera_face  = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
            bc_face       = ( mesh%domain(idom)%faces(ielem,iface)%ftype == BOUNDARY)

            if ( chimera_face ) then
                ChiID    = mesh%domain(idom)%faces(ielem,iface)%ChiID
                ncompute = mesh%domain(idom)%chimera%recv%data(ChiID)%ndonors()

            else if ( bc_face ) then
                ncompute = mesh%domain(idom)%faces(ielem,iface)%bc_ndepend

            else
                ! Standard conforming neighbor, only one dependent element.
                ncompute = 1
            end if

        end if


        end associate

    end function get_face_ncompute
    !**************************************************************************************





    !>  Return the number of elements that a volume function depends on in the linearization 
    !!  direction 'idiff'.
    !!
    !!  Often, this may be just 1. Because a flux linearized wrt its owner element just has 
    !!  1 dependent element; itself. A flux linearized wrt its neighbor element has just 
    !!  1 dependent element; it's neighbor. If however, a face is a Chimera face, then 
    !!  when linearized with respect to its donors, ndepend will be equal to the number of
    !!  donor elements that are providing information to that face.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_element_ncompute(self,mesh,elem_info,idiff) result(ncompute)
        class(equation_set_t),  intent(in)   :: self
        type(mesh_t),       intent(in)   :: mesh
        type(element_info_t),   intent(in)   :: elem_info
        integer(ik),            intent(in)   :: idiff

        integer(ik) :: ChiID, iface, ncompute
        logical     :: compute_wrt_none, compute_wrt_interior, compute_wrt_exterior, &
                       chimera_face

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l )


        compute_wrt_none     = (idiff == 0)
        compute_wrt_interior = (idiff == DIAG)
        compute_wrt_exterior = (idiff == XI_MIN   .or. idiff == XI_MAX  .or. &
                                idiff == ETA_MIN  .or. idiff == ETA_MAX .or. &
                                idiff == ZETA_MIN .or. idiff == ZETA_MAX)

        if (compute_wrt_none) then
            ncompute = 1

        else if (compute_wrt_interior) then
            ncompute = 1

        else if (compute_wrt_exterior) then
            ! Search iface in the direction of linearization
            iface = idiff

            chimera_face  = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )

            if ( chimera_face ) then
                ! only need to compute multiple times when we need the linearization of 
                ! the chimera neighbors
                ChiID  = mesh%domain(idom)%faces(ielem,iface)%ChiID
                ncompute = mesh%domain(idom)%chimera%recv%data(ChiID)%ndonors()
            else
                ! Standard conforming neighbor, only one dependent element.
                ncompute = 1
            end if

        end if


        end associate

    end function get_element_ncompute
    !**************************************************************************************









    !>  Compute a pseudo-timestep.
    !!
    !!  If none was allocated, use the default pseudo time-step based on element volume.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute_pseudo_timestep(self,idomain,mesh,sdata,cfln,itime)
        class(equation_set_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: idomain
        type(mesh_t),       intent(inout)   :: mesh
        type(solverdata_t),     intent(inout)   :: sdata
        real(rk),               intent(in)      :: cfln(:)
        integer(ik),            intent(in)      :: itime

        type(default_pseudo_timestep_t) :: default_timestep

        if (allocated(self%pseudo_timestep)) then
            call self%pseudo_timestep%compute(idomain,mesh,self%prop,sdata,cfln,itime)
        else
            call default_timestep%compute(idomain,mesh,self%prop,sdata,cfln,itime)
        end if

    end subroutine compute_pseudo_timestep
    !***************************************************************************************









    !>  Return the number of models allocated in the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function nmodels(self) result(nmodels_)
        class(equation_set_t),  intent(in)  :: self

        integer(ik) :: nmodels_


        if (allocated(self%models)) then
            nmodels_ = size(self%models)
        else
            nmodels_ = 0
        end if

    end function nmodels
    !**************************************************************************************


    


end module type_equation_set
