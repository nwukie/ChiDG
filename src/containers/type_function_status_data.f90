module type_function_status_data
#include <messenger.h>
    use mod_kinds,                      only: ik
    use mod_constants,                  only: NFUNCTION_TYPES, NFACES, NBLK
    use type_domain,                    only: domain_t
    use type_equationset_function_data, only: equationset_function_data_t
    implicit none



    
    !> Contains arrays of logicals that indicate the status of a given function. This class
    !! is used to manage the contributions and linearizations of functions.
    !!
    !! So, in this way, we can compute a flux on a face and add the contribution to both elements
    !! by registering that the function contribution was added. This container holds the registration
    !! information as arrays of logicals.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: function_status_data_t

        logical,    allocatable     ::  function_computed(:,:,:,:)                 !< fcn_type, ielem, iface, ifcn
        logical,    allocatable     ::  function_equation_computed(:,:,:,:,:)      !< fcn_type, ielem, iface, ifcn, ieqn
        logical,    allocatable     ::  function_linearized(:,:,:,:,:)             !< fcn_type, ielem, iface, ifcn, iblk
        logical,    allocatable     ::  function_equation_linearized(:,:,:,:,:,:)  !< fcn_type, ielem, iface, ifcn, iblk, ieqn

    contains

        procedure   :: init         !< Initialize storage for function status data
        procedure   :: clear        !< Reset function status data

    end type function_status_data_t
    !****************************************************************************************









contains







    !> Initialize storage arrays.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self, domain, function_data)
        class(function_status_data_t),      intent(inout)   :: self
        type(domain_t),                     intent(in)      :: domain
        type(equationset_function_data_t),  intent(in)      :: function_data

        integer(ik) :: nelem, nfcn, neqn, ierr

        nelem = domain%nelem
        neqn  = domain%neqns


        !
        ! Allocate boundary advective flux status storage
        !
        nfcn = max(function_data%nboundary_advective_flux, &
                   function_data%nboundary_diffusive_flux, &
                   function_data%nvolume_advective_flux,   &
                   function_data%nvolume_diffusive_flux,   &
                   function_data%nvolume_advective_source, &
                   function_data%nvolume_diffusive_source)




        ! Allocate
        allocate(self%function_computed(            NFUNCTION_TYPES,nelem,NFACES,nfcn), stat=ierr)
        if (ierr /= 0) call AllocationError
            
        allocate(self%function_equation_computed(   NFUNCTION_TYPES,nelem,NFACES,nfcn,neqn), stat=ierr)
        if (ierr /= 0) call AllocationError

        allocate(self%function_linearized(          NFUNCTION_TYPES,nelem,NFACES,nfcn,NBLK), stat=ierr)
        if (ierr /= 0) call AllocationError

        allocate(self%function_equation_linearized( NFUNCTION_TYPES,nelem,NFACES,nfcn,NBLK,neqn), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initialize
        self%function_computed = .false.
        self%function_equation_computed = .false.
        self%function_linearized = .false.
        self%function_equation_linearized = .false.




    end subroutine init
    !****************************************************************************************













    !> Reset all tracking of function status data. Should be executed before any evaluation of the 
    !! right-hand side.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine clear(self)
        class(function_status_data_t),  intent(inout)   :: self


        self%function_computed = .false.
        self%function_equation_computed = .false.
        self%function_linearized = .false.
        self%function_equation_linearized = .false.


    end subroutine clear
    !***************************************************************************************




















end module type_function_status_data
