module type_equationset
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: CHIMERA, DIAG, BOUNDARY_ADVECTIVE_FLUX, BOUNDARY_DIFFUSIVE_FLUX, &
                                              VOLUME_ADVECTIVE_FLUX, VOLUME_DIFFUSIVE_FLUX, XI_MIN, XI_MAX, &
                                              ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_equation,                  only: equation_t
    use type_properties,                only: properties_t
    use type_boundary_flux_wrapper,     only: boundary_flux_wrapper_t
    use type_volume_flux_wrapper,       only: volume_flux_wrapper_t
    use type_volume_flux,               only: volume_flux_t
    use type_boundary_flux,             only: boundary_flux_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_mesh,                      only: mesh_t
    use type_solverdata,                only: solverdata_t
    use type_element_info,              only: element_info_t
    use type_face_info,                 only: face_info_t
    use type_function_info,             only: function_info_t
    use mod_DNAD_tools,                 only: element_compute_seed, face_compute_seed
    implicit none
    private



    !>  Abstract equation-set type. Can be extended to implement a concrete equation set.
    !!      - Contains name and number of equations.
    !!      - Contains properties type with equations and material(ex. fluid) properties and definitions
    !!      - Contains arrays of flux components
    !!
    !!  When a new equation set is defined. It should extend from this abstract type. It must then 
    !!  implement the 'init' function where equations/variables can be added, and fluxes can be 
    !!  added to the definition of the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!  @note   Added support for diffusion terms
    !!
    !-------------------------------------------------------------------------------------------------
    type, public, abstract :: equationset_t
        character(len=:), allocatable       :: name 
        integer(ik)                         :: neqns

        ! Equation set properties
        class(properties_t), allocatable    :: prop


        ! Arrays of flux functions
        type(boundary_flux_wrapper_t),  allocatable   :: boundary_advective_flux(:)
        type(volume_flux_wrapper_t),    allocatable   :: volume_advective_flux(:)
        type(boundary_flux_wrapper_t),  allocatable   :: boundary_diffusive_flux(:)
        type(volume_flux_wrapper_t),    allocatable   :: volume_diffusive_flux(:) 


        ! Data for the flux and source functions. Ex how many.
        ! This gets passed to a container in sdata that keeps track of whether these
        ! have been executed or not.
        type(equationset_function_data_t)               :: function_data


    contains

        procedure(self_interface),  deferred  :: init   !< Initialization. This must be implemented by any actual equation set.

        procedure   :: set_name                         !< Set the name for the set of equations
        procedure   :: get_name                         !< Return the name fo the set of equations

        procedure   :: add_equation                     !< Add an equation, it's string, and index
        procedure   :: add_properties                   !< Add a properties type to the equation set

        procedure   :: add_boundary_advective_flux      !< Add a boundary advective function
        procedure   :: add_boundary_diffusive_flux      !< Add a boundary diffusive function
        procedure   :: add_volume_advective_flux        !< Add a volume advective function
        procedure   :: add_volume_diffusive_flux        !< Add a volume diffusive function

        procedure   :: compute_boundary_advective_flux  !< Compute all the boundary advective functions
        procedure   :: compute_boundary_diffusive_flux  !< Compute all the boundary diffusive functions
        procedure   :: compute_volume_advective_flux    !< Compute all the volume advective functions
        procedure   :: compute_volume_diffusive_flux    !< Compute all the volume diffusive functions

        procedure   :: get_boundary_ndependent_elements !< return the number of elements that a boundary function is depending on
        procedure   :: get_volume_ndependent_elements   !< return the number of elements that a volume function is depending on

    end type equationset_t
    !**************************************************************************************************







    !> Interface definitions
    abstract interface
        subroutine self_interface(self)
            import equationset_t
            class(equationset_t), intent(inout) :: self
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
    !---------------------------------------------------------------------------------------------------------
    subroutine set_name(self,name_string)
        class(equationset_t),   intent(inout)   :: self
        character(len=*),       intent(in)      :: name_string

        self%name = name_string

    end subroutine set_name
    !*********************************************************************************************************




    !>  Return the name of the equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------------
    function get_name(self) result(ename)
        class(equationset_t),   intent(in)   :: self

        character(len=:),   allocatable :: ename


        ename = self%name

    end function get_name
    !*********************************************************************************************************






    !> Add equationset%properties component
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  prop    properties_t class to be added
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine add_properties(self,prop)
        class(equationset_t),   intent(inout)   :: self
        class(properties_t),    intent(in)      :: prop
    
        integer(ik)     :: ierr

        if (allocated(self%prop)) then
            !
            ! If self%prop is already allocated, that is strange since only one is allowed per eqnset. Warn it is being replaced.
            !
            call chidg_signal(WARN,"equationset%add_properties: properties component was already allocated. Replacing current definition with incoming component.")

            !
            ! Deallocate current component.
            !
            deallocate(self%prop)

        end if


        !
        ! Allocate properties component
        !
        allocate(self%prop, source=prop, stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine add_properties
    !*********************************************************************************************************











    !>  Procedure to adding equations to the equation set properties
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  varstring   String defining the variable associated with the equation being added
    !!  @param[in]  varindex    The index of the equation in the given set. 
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine add_equation(self,varstring,varindex)
        class(equationset_t),   intent(inout)  :: self
        character(*),           intent(in)     :: varstring
        integer(ik),            intent(in)     :: varindex

        type(equation_t), allocatable    :: temp(:)
        integer(ik) :: ieq, ierr



        ! Check that properties storage has been allocated
        if (.not. allocated(self%prop)) call chidg_signal(FATAL,"Properties storage has not yet been allocated in the Equation Set. This must be done before adding equations since they are stored in the properties component")




        !
        ! If there are already equations allocated, reallocate and add new equation
        !
        if (allocated(self%prop%eqns)) then
            !
            ! Allocate temp eqn array with one extra slot for new eqn
            !
            allocate(temp(size(self%prop%eqns) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy current eqns to first temp slots
            !
            do ieq = 1,size(self%prop%eqns)
                temp(ieq) = self%prop%eqns(ieq)
            end do


            !
            ! Add new eqn to last slot
            !
            temp(size(temp))%name = varstring
            temp(size(temp))%ind  = varindex


            !
            ! Store temp equation array to equation properties
            !
            self%prop%eqns = temp

        !
        ! If there are no equations allocated, allocate one slot and set data
        !
        else
            !
            ! Allocate equation
            !
            allocate(self%prop%eqns(1), stat=ierr)
            if (ierr /= 0) call AllocationError

            self%prop%eqns(1)%name = varstring
            self%prop%eqns(1)%ind  = varindex

        end if



        !
        ! Resize neqns
        !
        self%neqns = size(self%prop%eqns)

    end subroutine add_equation
    !***************************************************************************************************************










    !> Add components to volume_advective_flux array
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  flux    Volume advective flux component to be added
    !!
    !-----------------------------------------------------------------------------------------------------------------------------
    subroutine add_volume_advective_flux(self,flux)
        class(equationset_t),   intent(inout)   :: self
        class(volume_flux_t),   intent(in)      :: flux
    
        class(volume_flux_wrapper_t), allocatable   :: temp(:)
        integer(ik)     :: ierr, iflux

        if (allocated(self%volume_advective_flux)) then

            !
            ! Allocate temporary flux array with one additional slot
            !
            allocate(temp(size(self%volume_advective_flux) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy current flux components to temp array
            !
            do iflux = 1,size(self%volume_advective_flux)
                allocate(temp(iflux)%flux,source=self%volume_advective_flux(iflux)%flux, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do


            !
            ! Add new flux to last slot
            !
            allocate(temp(size(temp))%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy temp array back to equationset
            !
            self%volume_advective_flux = temp

        else
            !
            ! Allocate one slot
            !
            allocate(self%volume_advective_flux(1), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Allocate flux component from source
            !
            allocate(self%volume_advective_flux(1)%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

        end if


        !
        ! Update function data
        !
        self%function_data%nvolume_advective_flux = size(self%volume_advective_flux)

    end subroutine add_volume_advective_flux
    !*****************************************************************************************************************************












    !> Add components to volume_advective_flux array
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  flux    Volume advective flux component to be added
    !!
    !-----------------------------------------------------------------------------------------------------------------------------
    subroutine add_volume_diffusive_flux(self,flux)
        class(equationset_t),   intent(inout)   :: self
        class(volume_flux_t),   intent(in)      :: flux
    
        class(volume_flux_wrapper_t), allocatable   :: temp(:)
        integer(ik)     :: ierr, iflux

        if (allocated(self%volume_diffusive_flux)) then

            !
            ! Allocate temporary flux array with one additional slot
            !
            allocate(temp(size(self%volume_diffusive_flux) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy current flux components to temp array
            !
            do iflux = 1,size(self%volume_diffusive_flux)
                allocate(temp(iflux)%flux,source=self%volume_diffusive_flux(iflux)%flux, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do


            !
            ! Add new flux to last slot
            !
            allocate(temp(size(temp))%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy temp array back to equationset
            !
            self%volume_diffusive_flux = temp

        else
            !
            ! Allocate one slot
            !
            allocate(self%volume_diffusive_flux(1), stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Allocate flux component from source
            !
            allocate(self%volume_diffusive_flux(1)%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

        end if


        !
        ! Update function data
        !
        self%function_data%nvolume_diffusive_flux = size(self%volume_diffusive_flux)

    end subroutine add_volume_diffusive_flux
    !*****************************************************************************************************************************








    !> Add components to boundary_advective_flux array
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  flux    Boundary advective flux type to be added
    !!
    !-----------------------------------------------------------------------------------------------------------------------------
    subroutine add_boundary_advective_flux(self,flux)
        class(equationset_t),   intent(inout)   :: self
        class(boundary_flux_t), intent(in)      :: flux
    
        class(boundary_flux_wrapper_t), allocatable   :: temp(:)
        integer(ik)     :: ierr, iflux

        if (allocated(self%boundary_advective_flux)) then

            !
            ! Allocate temporary flux array with one additional slot
            !
            allocate(temp(size(self%boundary_advective_flux) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Copy current flux components to temp array
            !
            do iflux = 1,size(self%boundary_advective_flux)
                allocate(temp(iflux)%flux,source=self%boundary_advective_flux(iflux)%flux, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do


            !
            ! Add new flux to last slot
            !
            allocate(temp(size(temp))%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy temp array back to equationset
            !
            self%boundary_advective_flux = temp

        else

            !
            ! Allocate new slot
            !
            allocate(self%boundary_advective_flux(1), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Allocate flux in new wrapper slot
            !
            allocate(self%boundary_advective_flux(1)%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

        end if



        !
        ! Update function data
        !
        self%function_data%nboundary_advective_flux = size(self%boundary_advective_flux)


    end subroutine add_boundary_advective_flux
    !*****************************************************************************************************************************











    !> Add components to boundary_diffusive_flux array
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!  @param[in]  flux    Boundary diffusive flux type to be added
    !!
    !-----------------------------------------------------------------------------------------------------------------------------
    subroutine add_boundary_diffusive_flux(self,flux)
        class(equationset_t),   intent(inout)   :: self
        class(boundary_flux_t), intent(in)      :: flux
    
        class(boundary_flux_wrapper_t), allocatable   :: temp(:)
        integer(ik)     :: ierr, iflux

        if (allocated(self%boundary_diffusive_flux)) then

            !
            ! Allocate temporary flux array with one additional slot
            !
            allocate(temp(size(self%boundary_diffusive_flux) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Copy current flux components to temp array
            !
            do iflux = 1,size(self%boundary_diffusive_flux)
                allocate(temp(iflux)%flux,source=self%boundary_diffusive_flux(iflux)%flux, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do


            !
            ! Add new flux to last slot
            !
            allocate(temp(size(temp))%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

            !
            ! Copy temp array back to equationset
            !
            self%boundary_diffusive_flux = temp

        else

            !
            ! Allocate new slot
            !
            allocate(self%boundary_diffusive_flux(1), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Allocate flux in new wrapper slot
            !
            allocate(self%boundary_diffusive_flux(1)%flux, source=flux, stat=ierr)
            if (ierr /= 0) call AllocationError

        end if



        !
        ! Update function data
        !
        self%function_data%nboundary_diffusive_flux = size(self%boundary_diffusive_flux)


    end subroutine add_boundary_diffusive_flux
    !*****************************************************************************************************************************








    !>  Loops through the attached boundary advective flux functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine compute_boundary_advective_flux(self,mesh,sdata,face_info,iblk)
        class(equationset_t),   intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(solverdata_t),     intent(inout)   :: sdata
        type(face_info_t),      intent(inout)   :: face_info
        integer(ik),            intent(in)      :: iblk

        integer(ik)             :: nfcn, ifcn, idepend, ndepend
        type(function_info_t)   :: function_info
        logical                 :: compute_function, linearize_function

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface, prop => self%prop )

        !
        ! Get number of elements we are linearizing with respect to
        !
        ndepend = self%get_boundary_ndependent_elements(mesh,face_info,iblk)



        if (allocated(self%boundary_advective_flux)) then
            nfcn = size(self%boundary_advective_flux)
            do ifcn = 1,nfcn

                function_info%type   = BOUNDARY_ADVECTIVE_FLUX
                function_info%ifcn   = ifcn
                function_info%iblk   = iblk

                compute_function     = sdata%function_status%compute_function(   face_info, function_info )
                linearize_function   = sdata%function_status%linearize_function( face_info, function_info )
                

                if ( compute_function .or. linearize_function ) then
                    !
                    ! Compute boundary flux once for each donor. 
                    !   - For interior faces ndepend == 1. 
                    !   - For Chimera faces ndepend is potentially > 1.
                    !
                    do idepend = 1,ndepend
                        function_info%seed    = face_compute_seed(mesh,idom,ielem,iface,idepend,iblk)
                        function_info%idepend = idepend

                        call self%boundary_advective_flux(ifcn)%flux%compute(mesh,sdata,prop,face_info,function_info)
                    end do
                end if


            end do ! ifcn

        end if ! boundary_advective_flux loop

        end associate

    end subroutine compute_boundary_advective_flux
    !*****************************************************************************************************************









    !>  Loops through the attached boundary diffusive flux functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine compute_boundary_diffusive_flux(self,mesh,sdata,face_info,iblk)
        class(equationset_t),   intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(solverdata_t),     intent(inout)   :: sdata
        type(face_info_t),      intent(inout)   :: face_info
        integer(ik),            intent(in)      :: iblk

        integer(ik)             :: nfcn, ifcn, idepend, ndepend
        type(function_info_t)   :: function_info
        logical                 :: compute_function, linearize_function

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface, prop => self%prop)


        !
        ! Get number of elements we are linearizing with respect to
        !
        ndepend = self%get_boundary_ndependent_elements(mesh,face_info,iblk)


        if (allocated(self%boundary_diffusive_flux)) then
            nfcn = size(self%boundary_diffusive_flux)
            do ifcn = 1,nfcn

                function_info%type   = BOUNDARY_DIFFUSIVE_FLUX
                function_info%ifcn   = ifcn
                function_info%iblk   = iblk

                compute_function     = sdata%function_status%compute_function(   face_info, function_info )
                linearize_function   = sdata%function_status%linearize_function( face_info, function_info )
                

                if ( compute_function .or. linearize_function ) then
                    !
                    ! Compute boundary flux once for each donor. 
                    !   - For interior faces ndepend == 1. 
                    !   - For Chimera faces ndepend is potentially > 1.
                    !
                    do idepend = 1,ndepend
                        function_info%seed    = face_compute_seed(mesh,idom,ielem,iface,idepend,iblk)
                        function_info%idepend = idepend

                        call self%boundary_diffusive_flux(ifcn)%flux%compute(mesh,sdata,prop,face_info,function_info)
                    end do
                end if


            end do ! ifcn

        end if ! boundary_diffusive_flux loop

        end associate

    end subroutine compute_boundary_diffusive_flux
    !*****************************************************************************************************************













    !>  Loops through the attached volume advective flux functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine compute_volume_advective_flux(self,mesh,sdata,elem_info,iblk)
        class(equationset_t),       intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(solverdata_t),         intent(inout)   :: sdata
        type(element_info_t),       intent(inout)   :: elem_info
        integer(ik),                intent(in)      :: iblk

        integer(ik)             :: nfcn, ifcn, idepend
        type(function_info_t)   :: function_info
        logical                 :: compute_flux

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l, prop => self%prop)


        !
        ! Volume advective flux only depends on the interior element
        !
        idepend = 1


        compute_flux = ( (iblk == DIAG) .and. allocated(self%volume_advective_flux) )
        if (compute_flux) then
            nfcn = size(self%volume_advective_flux)
            do ifcn = 1,nfcn

                function_info%type    = VOLUME_ADVECTIVE_FLUX
                function_info%ifcn    = ifcn
                function_info%iblk    = iblk
                function_info%idepend = idepend
                function_info%seed    = element_compute_seed(mesh,idom,ielem,idepend,iblk)

                call self%volume_advective_flux(ifcn)%flux%compute(mesh,sdata,prop,elem_info,function_info)

            end do ! ifcn
        end if

        end associate


    end subroutine compute_volume_advective_flux
    !******************************************************************************************************************










    !>  Loops through the attached volume diffusive flux functions and computes them 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine compute_volume_diffusive_flux(self,mesh,sdata,elem_info,iblk)
        class(equationset_t),       intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(solverdata_t),         intent(inout)   :: sdata
        type(element_info_t),       intent(inout)   :: elem_info
        integer(ik),                intent(in)      :: iblk

        integer(ik)             :: nfcn, ifcn, idepend, ndepend
        type(function_info_t)   :: function_info
        logical                 :: compute_function, linearize_function

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l, prop => self%prop)


        !
        ! Get number of elements we are linearizing with respect to
        !
        ndepend = self%get_volume_ndependent_elements(mesh,elem_info,iblk)


        if (allocated(self%volume_diffusive_flux)) then
            nfcn = size(self%volume_diffusive_flux)
            do ifcn = 1,nfcn

                !
                ! Compute boundary flux once for each donor. 
                !   - For interior faces ndepend == 1. 
                !   - For Chimera faces ndepend is potentially > 1.
                !
                do idepend = 1,ndepend

                    function_info%type    = VOLUME_DIFFUSIVE_FLUX
                    function_info%ifcn    = ifcn
                    function_info%iblk    = iblk
                    function_info%idepend = idepend
                    function_info%seed    = element_compute_seed(mesh,idom,ielem,idepend,iblk)

                    call self%volume_diffusive_flux(ifcn)%flux%compute(mesh,sdata,prop, elem_info, function_info)

                end do

            end do ! ifcn
        end if

        end associate


    end subroutine compute_volume_diffusive_flux
    !******************************************************************************************************************








    !>  Return the number of elements that a boundary function depends on in the linearization direction 'iblk'.
    !!
    !!  Often, this may be just 1. Because a flux linearized wrt its owner element just has 1 dependent element; itself.
    !!  A flux linearized wrt its neighbor element has just 1 dependent element; it's neighbor. If however, a face
    !!  is a Chimera face, then when linearized with respect to its donors, ndepend will be equal to the number of
    !!  donor elements that are providing information to that face.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/15/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------
    function get_boundary_ndependent_elements(self,mesh,face_info,iblk) result(ndepend)
        class(equationset_t),   intent(in)   :: self
        type(mesh_t),           intent(in)   :: mesh(:)
        type(face_info_t),      intent(in)   :: face_info
        integer(ik),            intent(in)   :: iblk

        integer(ik) :: ChiID, ndepend
        logical     :: depend_me, depend_neighbor, chimera_face

        associate( idom => face_info%idomain_l, ielem => face_info%ielement_l, iface => face_info%iface )

        depend_me       = (iblk == DIAG)
        depend_neighbor = (iblk == XI_MIN   .or. iblk == XI_MAX  .or. &
                           iblk == ETA_MIN  .or. iblk == ETA_MAX .or. &
                           iblk == ZETA_MIN .or. iblk == ZETA_MAX)


        if (depend_me) then

            ndepend = 1



        else if (depend_neighbor) then

            chimera_face  = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )

            if ( chimera_face ) then
                ! only need to compute multiple times when we need the linearization of the chimera neighbors
                ChiID  = mesh(idom)%faces(ielem,iface)%ChiID
                ndepend = mesh(idom)%chimera%recv%data(ChiID)%ndonors()
            else
                ! Standard conforming neighbor, only one dependent element.
                ndepend = 1
            end if

        end if


        end associate

    end function get_boundary_ndependent_elements
    !*****************************************************************************************************************





    !>  Return the number of elements that a volume function depends on in the linearization direction 'iblk'.
    !!
    !!  Often, this may be just 1. Because a flux linearized wrt its owner element just has 1 dependent element; itself.
    !!  A flux linearized wrt its neighbor element has just 1 dependent element; it's neighbor. If however, a face
    !!  is a Chimera face, then when linearized with respect to its donors, ndepend will be equal to the number of
    !!  donor elements that are providing information to that face.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------------
    function get_volume_ndependent_elements(self,mesh,elem_info,iblk) result(ndepend)
        class(equationset_t),   intent(in)   :: self
        type(mesh_t),           intent(in)   :: mesh(:)
        type(element_info_t),   intent(in)   :: elem_info
        integer(ik),            intent(in)   :: iblk

        integer(ik) :: ChiID, ndepend, iface
        logical     :: depend_me, depend_neighbor, chimera_face

        associate( idom => elem_info%idomain_l, ielem => elem_info%ielement_l )


        depend_me       = (iblk == DIAG)
        depend_neighbor = (iblk == XI_MIN   .or. iblk == XI_MAX  .or. &
                           iblk == ETA_MIN  .or. iblk == ETA_MAX .or. &
                           iblk == ZETA_MIN .or. iblk == ZETA_MAX)


        if (depend_me) then

            ndepend = 1



        else if (depend_neighbor) then
            ! Search iface in the direction of linearization
            iface = iblk

            chimera_face  = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )

            if ( chimera_face ) then
                ! only need to compute multiple times when we need the linearization of the chimera neighbors
                ChiID  = mesh(idom)%faces(ielem,iface)%ChiID
                ndepend = mesh(idom)%chimera%recv%data(ChiID)%ndonors()
            else
                ! Standard conforming neighbor, only one dependent element.
                ndepend = 1
            end if

        end if


        end associate

    end function get_volume_ndependent_elements
    !*****************************************************************************************************************


end module type_equationset
