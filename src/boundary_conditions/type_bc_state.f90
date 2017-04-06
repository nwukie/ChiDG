module type_bc_state
#include <messenger.h>
    use mod_kinds,                  only: rk, ik

    use type_bcproperty_set,        only: bcproperty_set_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use type_mesh,              only: mesh_t
    use type_bc_patch,              only: bc_patch_t
    use mpi_f08,                    only: mpi_comm
    implicit none




    !>  Abstract base-type for computing a boundary condition state
    !!
    !!      - contains a procedure for computing the bc state: compute_bc_state
    !!      - compute_bc_state is deferred and so must be implemented by any new bc_state_t
    !!      - bc_state_t also contains properties that can hold parameters and functions 
    !!        that have been set for the boundary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/9/2016
    !!  @note   Changed boundary condition to compute a bc state
    !!
    !--------------------------------------------------------------------------------------------
    type, public, abstract :: bc_state_t

        character(:),   allocatable :: name
        character(:),   allocatable :: family

        ! Boundary condition options
        type(bcproperty_set_t)      :: bcproperties

    contains

        procedure(bc_state_init),       deferred :: init
        procedure(bc_state_compute),    deferred :: compute_bc_state

        procedure   :: init_bc_specialized
        procedure   :: init_bc_coupling

        procedure   :: set_name
        procedure   :: get_name
        procedure   :: set_family
        procedure   :: get_family

        procedure   :: set_fcn               ! Set a particular function definition for a specified bcfunction_t
        procedure   :: set_fcn_option        ! Set function-specific options for a specified bcfunction_t

        procedure   :: get_nproperties       ! Return the number of properties associated with the boundary condition.
        procedure   :: get_property_name     ! Return the name of a property given a property index.
        procedure   :: get_noptions          ! Return the number of available options for a given property, specified by a property index.
        procedure   :: get_option_key        ! Return the key for an option, given a property index and subsequent option index.
        procedure   :: get_option_value      ! Return the value of a given key, inside of a specified property.

    end type bc_state_t
    !*********************************************************************************************




    abstract interface
        subroutine bc_state_init(self)
            import bc_state_t

            class(bc_state_t),  intent(inout)   :: self
        end subroutine
    end interface



    abstract interface
        subroutine bc_state_compute(self,worker,prop)
            import bc_state_t
            import chidg_worker_t
            import properties_t

            class(bc_state_t),      intent(inout)   :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
        end subroutine
    end interface


contains



    !>  Default specialized initialization procedure. This is called from the base bc%init procedure
    !!  and can be overwritten by derived types to implement specialized initiailization details.
    !!
    !!  By default, this routine does nothing. However, a particular bc_state_t could reimplement
    !!  this routine to perform some specialized initialization calculations during initialization.
    !!
    !!  For example, a point pressure outlet boundary condition may want to find a particular 
    !!  quadrature node to set pressure at. init_bc_specialized could be defined for that
    !!  bc_state_t implementation to search the quadrature nodes over all the bc_patch faces
    !!  to find the correct node to set the pressure at.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/21/2017
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_specialized(self,mesh,bc_patch,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),   intent(in)      :: mesh
        type(bc_patch_t),   intent(in)      :: bc_patch(:)
        type(mpi_comm),     intent(in)      :: bc_COMM



    end subroutine init_bc_specialized
    !**********************************************************************************************





    !>  Default boundary coupling initialization routine. 
    !!
    !!  Default initializes coupling for a given element to just itself and no coupling with 
    !!  other elements on the boundary. For a boundary condition that is coupled across the face
    !!  this routine can be overwritten to set the coupling information specific to the boundary 
    !!  condition.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!  @date   2/27/2017   updated for multiple patches
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,bc_patch)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),   intent(in)      :: mesh
        type(bc_patch_t),   intent(inout)   :: bc_patch(:)

        integer(ik) :: ipatch, iface_bc, idomain_g, idomain_l, ielement_g, ielement_l



        !
        ! For each patch, loop through faces and set default element coupling.
        ! Default is that each face is coupled only with its owner element.
        ! So, strictly local coupling.
        !
        do ipatch = 1,size(bc_patch)
            do iface_bc = 1,bc_patch(ipatch)%nfaces()


                !
                ! Get block-element index of current iface_bc
                !
                idomain_g  = bc_patch(ipatch)%idomain_g(iface_bc)
                idomain_l  = bc_patch(ipatch)%idomain_l(iface_bc)
                ielement_g = bc_patch(ipatch)%ielement_g(iface_bc)
                ielement_l = bc_patch(ipatch)%ielement_l(iface_bc)

                
                !
                ! Add the element index as the only dependency.
                !
                call bc_patch(ipatch)%add_coupled_element(iface_bc, idomain_g,  &
                                                                    idomain_l,  &
                                                                    ielement_g, &
                                                                    ielement_l, &
                                                                    IRANK)


            end do ! iface_bc
        end do !ipatch


    end subroutine init_bc_coupling
    !**********************************************************************************************


















    !>  Set a function for a specified property.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcprop  String specifying a bcproperty_t to edit.
    !!  @param[in]  fcn     String specifying the concrete function_t to set.
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set_fcn(self,bcprop,fcn)
        class(bc_state_t),      intent(inout)   :: self
        character(*),           intent(in)      :: bcprop
        character(*),           intent(in)      :: fcn


        call self%bcproperties%set_fcn(bcprop,fcn)


    end subroutine set_fcn
    !*********************************************************************************************









    !>  Set a function option for a specified property.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcprop  String specifying a bcproperty_t to edit.
    !!  @param[in]  option  String specifying a particular option within bcproperty_f%fcn to edit
    !!  @param[in]  val     Real value to be set for the option.
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine set_fcn_option(self,bcprop,option,val)
        class(bc_state_t),      intent(inout)   :: self
        character(*),           intent(in)      :: bcprop
        character(*),           intent(in)      :: option
        real(rk),               intent(in)      :: val

        call self%bcproperties%set_fcn_option(bcprop,option,val)

    end subroutine set_fcn_option
    !************************************************************************************************








    !>  Return number of properties available in the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    function get_nproperties(self) result(nprop)
        class(bc_state_t),    intent(in)  :: self

        integer(ik) :: nprop

        nprop = self%bcproperties%get_nproperties()

    end function get_nproperties
    !***************************************************************************************************








    !>  Return a property name string, given the index of the property in the boundary condition.
    !!
    !!  This probably works best by first calling get_nproperties, and then iterating through the 
    !!  number of available properties to get their names.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop   Integer specifying the index of the property to be queried.
    !!  @result     pname   String of the property name associated with the index.
    !!
    !---------------------------------------------------------------------------------------------------
    function get_property_name(self,iprop) result(pname)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop

        character(len=:),   allocatable :: pname

        pname = self%bcproperties%get_property_name(iprop) 

    end function get_property_name
    !***************************************************************************************************












    !>  Return an option key, given a property index and option index. 
    !!
    !!  One probably calls get_noptions(iprop)
    !!  first, to get the number of available options for the function currently set for the property 'iprop'.
    !!  Then one can loop over the number of available options and return their availble names dynamically.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!  @param[in]  iprop       Integer index of a property to modify.
    !!  @param[in]  ioption     Integer index of an option inside bcproperty%fcn
    !!  @result     key         String(key) corresponding to the option index (ioption)
    !!
    !----------------------------------------------------------------------------------------------------
    function get_option_key(self,iprop,ioption) result(key)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop
        integer(ik),        intent(in)  :: ioption

        character(len=:),   allocatable :: key

        key = self%bcproperties%bcprop(iprop)%get_option_key(ioption)

    end function get_option_key
    !****************************************************************************************************










    !>  Return an option value, given a property index and option key.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop       Integer index of a property to modify.
    !!  @param[in]  key         String(key) specifying the option to be queried.
    !!  @result     val         Returned value of the selected key.
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_option_value(self,iprop,key) result(val)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop
        character(*),       intent(in)  :: key

        real(rk)        :: val

        val = self%bcproperties%bcprop(iprop)%get_option_value(key)

    end function get_option_value
    !****************************************************************************************************









    !>  Return the number of available options, given a property index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop       Integer index of a property to query.
    !!  @result     noption     Returned number of options available for the property. Dependends on 
    !!                          the function that is set for the property.
    !!
    !---------------------------------------------------------------------------------------------------
    function get_noptions(self,iprop) result(noptions)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop

        integer(ik)     :: noptions

        noptions = self%bcproperties%bcprop(iprop)%get_noptions()

    end function get_noptions
    !***************************************************************************************************





    !>  Set the boundary condition name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_name(self,bcname)
        class(bc_state_t),  intent(inout)   :: self
        character(*),       intent(in)      :: bcname

        self%name = trim(bcname)

    end subroutine set_name
    !***************************************************************************************************








    !>  Return the boundary condition name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function get_name(self) result(bcname)
        class(bc_state_t),    intent(in)  :: self

        character(:),   allocatable :: bcname

        bcname = self%name

    end function get_name
    !***************************************************************************************************







    !>  Set the boundary condition family.
    !!
    !!  Allowable families:
    !!      - Inlet
    !!      - Oulet
    !!      - Wall
    !!      - Symmetry
    !!      - Periodic
    !!      - Farfield
    !!      - Scalar
    !!      - Extrapolation
    !!      - Empty
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/21/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_family(self,bc_family)
        class(bc_state_t),  intent(inout)   :: self
        character(*),       intent(in)      :: bc_family

        character(:),   allocatable :: user_msg


        !
        ! Check incoming bc_state family
        !
        if ( (trim(bc_family) == 'Inlet')           .or. &
             (trim(bc_family) == 'Outlet')          .or. &
             (trim(bc_family) == 'Wall')            .or. &
             (trim(bc_family) == 'Symmetry')        .or. &
             (trim(bc_family) == 'Periodic')        .or. &
             (trim(bc_family) == 'Farfield')        .or. &
             (trim(bc_family) == 'Scalar')          .or. &
             (trim(bc_family) == 'Extrapolation')   .or. &
             (trim(bc_family) == 'Empty') ) then


            self%family = trim(bc_family)

        else

             user_msg = "bc_state%set_family: An invalid Family was trying to be set for the &
                         bc_state. Valid Families are: 'Inlet', 'Outlet', 'Wall', Symmetry', &
                         'Periodic', 'Farfield', 'Scalar', 'Extrapolation'."
             call chidg_signal_one(FATAL,user_msg,trim(bc_family))

        end if



    end subroutine set_family
    !***************************************************************************************************








    !>  Return the boundary condition family.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/21/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function get_family(self) result(bc_family)
        class(bc_state_t),    intent(in)  :: self

        character(:),   allocatable :: bc_family, user_msg


        user_msg = "bc_state%get_family: It looks like the Family component for the bc_state was &
                    not set. Make sure self%set_family('my_family') is being called in the bc_state &
                    initialization procedure."
        if (.not. allocated(self%family)) call chidg_signal(FATAL,user_msg)

        bc_family = self%family

    end function get_family
    !***************************************************************************************************





end module type_bc_state
