module type_bc_state
#include <messenger.h>
    use mod_kinds,                  only: rk, ik

    use type_bcproperty_set,        only: bcproperty_set_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use type_mesh,                  only: mesh_t
    use type_bc_patch,              only: bc_patch_t
    implicit none




    !> Abstract base-type for boundary conditions
    !!  - contains a list of associated element indices
    !!  - contains a list of face indices
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

        character(len=:), allocatable   :: name

        ! Boundary condition options
        type(bcproperty_set_t)          :: bcproperties

    contains

        procedure(bc_state_init),       deferred :: init
        procedure(bc_state_compute),    deferred :: compute_bc_state

        procedure   :: init_bc_coupling

        procedure   :: set_name
        procedure   :: get_name

        procedure   :: set_fcn               !< Set a particular function definition for a specified bcfunction_t
        procedure   :: set_fcn_option        !< Set function-specific options for a specified bcfunction_t


        procedure   :: get_nproperties       !< Return the number of properties associated with the boundary condition.
        procedure   :: get_property_name     !< Return the name of a property given a property index.
        procedure   :: get_noptions          !< Return the number of available options for a given property, specified by a property index.
        procedure   :: get_option_key        !< Return the key for an option, given a property index and subsequent option index.
        procedure   :: get_option_value      !< Return the value of a given key, inside of a specified property.

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
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,bc_patch)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(in)      :: mesh
        type(bc_patch_t),   intent(inout)   :: bc_patch

        integer(ik) :: iface_bc, ielem



        !
        ! Loop through elements and set default coupling information
        !
        do iface_bc = 1,bc_patch%nfaces()


            !
            ! Get block-element index of current iface_bc
            !
            ielem = bc_patch%ielement_l(iface_bc)

            
            !
            ! Add the element index as the only dependency.
            !
            call bc_patch%coupled_elements(iface_bc)%push_back(ielem)


        end do ! iface_bc


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
        class(bc_state_t),            intent(inout)   :: self
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
        class(bc_state_t),            intent(inout)   :: self
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
        class(bc_state_t),    intent(in)  :: self
        integer(ik),    intent(in)  :: iprop

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
        class(bc_state_t),    intent(inout)   :: self
        integer(ik),    intent(in)      :: iprop
        integer(ik),    intent(in)      :: ioption

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
        class(bc_state_t),    intent(inout)  :: self
        integer(ik),    intent(in)  :: iprop
        character(*),   intent(in)  :: key

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
        class(bc_state_t),    intent(inout)  :: self
        integer(ik),    intent(in)  :: iprop

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
        class(bc_state_t),    intent(inout)   :: self
        character(*),   intent(in)      :: bcname

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

        character(len=:), allocatable :: bcname

        bcname = self%name

    end function get_name
    !***************************************************************************************************







end module type_bc_state
