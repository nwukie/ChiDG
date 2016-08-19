module type_bc
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, BOUNDARY, CHIMERA, ORPHAN, BC_BLK, ZERO, ONE, TWO, RKTOL
    use mod_chidg_mpi,              only: IRANK

    use type_mesh,                  only: mesh_t
    use type_point,                 only: point_t
    use type_ivector,               only: ivector_t
    use type_solverdata,            only: solverdata_t
    use type_properties,            only: properties_t
    use type_face_info,             only: face_info_t
    use type_function_info,         only: function_info_t
    use type_bcproperty_set,        only: bcproperty_set_t
    use type_function,              only: function_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    implicit none
    private




    !> Abstract base-type for boundary conditions
    !!  - contains a list of associated element indices
    !!  - contains a list of face indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !--------------------------------------------------------------------------------------------
    type, public, abstract :: bc_t

        character(len=:),   allocatable :: name
        logical,    public              :: isInitialized = .false.  !< Logical switch for indicating the boundary condition initializaiton status

        !
        ! Boundary condition geometry
        !
        integer(ik),        allocatable :: dom(:)                   !< Indices of domains
        integer(ik),        allocatable :: elems(:)                 !< Indices of elements associated with boundary condition. Block-local indices. Local to partition.
        integer(ik),        allocatable :: faces(:)                 !< Indices of the boundary face for elements elems(ielems)
        type(ivector_t),    allocatable :: coupled_elems(:)         !< For each element on the boundary, a vector of element block-indices coupled with the current element.

        !        integer(ik),        allocatable :: elems_g(:)               !< Indices of elements associated with boundary condition. Block-global indices. 

        !
        ! Boundary condition options
        !
        !type(bcparameter_set_t)    :: bcparameters
        type(bcproperty_set_t)      :: bcproperties

    contains



        procedure                               :: init                     !< Boundary condition initialization
        procedure                               :: init_spec                !< Call specialized initialization routine
        procedure                               :: init_boundary_coupling   !< Initialize book-keeping for coupling interaction between elements.
        procedure(compute_interface), deferred  :: compute                  !< Implements boundary condition function
        procedure                               :: apply                    !< Apply bc function over bc elements



        procedure   :: set_name                                 !< Set the boundary condition name
        procedure   :: get_name                                 !< Return the boundary condition name

        
        procedure   :: add_options                              !< Specialized by each bc_t implementation. Adds options available

        procedure   :: set_fcn                                  !< Set a particular function definition for a specified bcfunction_t
        procedure   :: set_fcn_option                           !< Set function-specific options for a specified bcfunction_t


        procedure   :: get_nproperties                          !< Return the number of properties associated with the boundary condition.
        procedure   :: get_property_name                        !< Return the name of a property given a property index.
        procedure   :: get_noptions                             !< Return the number of available options for a given property, specified by a property index.
        procedure   :: get_option_key                           !< Return the key for an option, given a property index and subsequent option index.
        procedure   :: get_option_value                         !< Return the value of a given key, inside of a specified property.
        procedure   :: get_ncoupled_elems                       !< Return the number of elements coupled with a specified boundary element.

    end type bc_t
    !*********************************************************************************************



    abstract interface
        subroutine compute_interface(self,mesh,sdata,prop,face,fcn)
            use mod_kinds,  only: ik
            import bc_t
            import mesh_t
            import solverdata_t
            import properties_t
            import face_info_t
            import function_info_t

            class(bc_t),            intent(inout)   :: self
            type(mesh_t),           intent(in)      :: mesh(:)
            type(solverdata_t),     intent(inout)   :: sdata
            class(properties_t),    intent(inout)   :: prop
            type(face_info_t),      intent(in)      :: face
            type(function_info_t),  intent(in)      :: fcn
        end subroutine
    end interface



contains

    !> Initialize boundary condition routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!  @param[in]  mesh    mesh_t object containing elements and faces
    !!  @param[in]  iface   block face index to which the boundary condition is being applied
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,bconnectivity)
        class(bc_t),                    intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        type(boundary_connectivity_t),  intent(in)      :: bconnectivity

        type(point_t)                   :: pnt, point_one, point_two, point_three
        character(len=:),   allocatable :: bcname        
        real(rk)                        :: time, x, y, z
        integer(ik)                     :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc, & 
                                           xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end, & 
                                           ixi, ieta, izeta, ierr, ielem, ielem_test, nface_nodes, iface, inode, i, nfaces_bc, iface_bc

        logical,    allocatable         :: node_matched(:), xi_face, eta_face, zeta_face
        
        integer(ik), allocatable    :: element_nodes(:)
        integer(ik)                 :: face_node


        !
        ! Get number of elements/faces associated with boundary condition.
        !
        nelem_bc = bconnectivity%get_nfaces()




        !
        ! Detect number of bcfaces included in the mesh. Could only be a partial match because of parallel partitioning. Sets nfaces_bc
        !
        nfaces_bc = 0
        do ielem_bc = 1,nelem_bc

            ! Allocate array to register if bc node is included in element node list
            nface_nodes = size(bconnectivity%data(ielem_bc)%data)
            if ( allocated(node_matched) ) deallocate(node_matched)
            allocate(node_matched(nface_nodes), stat=ierr)
            if (ierr /= 0) call AllocationError


            ! Search for element with common nodes
            do ielem = 1,mesh%nelem

                ! Loop through bc nodes and see if current element contains them
                node_matched = .false.
                do inode = 1,nface_nodes
                    element_nodes = mesh%elems(ielem)%connectivity%get_element_nodes()
                    face_node     = bconnectivity%data(ielem_bc)%get_face_node(inode)
                    if ( any(element_nodes == face_node) ) then
                        node_matched(inode) = .true.
                    end if
                end do

                ! If all match, increment number of boundary faces to allocate
                if ( all(node_matched) ) then
                    nfaces_bc = nfaces_bc + 1 
                end if

            end do ! ielem

        end do ! ielem_bc








        !
        ! Allocate storage for element and face indices
        !
        allocate(self%elems(nfaces_bc), self%faces(nfaces_bc), self%coupled_elems(nfaces_bc), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Loop through each face in bc connectivity and call initialization for each face in local mesh
        !
        nelem_bc = bconnectivity%get_nfaces()
        iface_bc = 0
        do ielem_bc = 1,nelem_bc

            nface_nodes = size(bconnectivity%data(ielem_bc)%data)
            if ( allocated(node_matched) ) deallocate(node_matched)
            allocate(node_matched(nface_nodes), stat=ierr)
            if (ierr /= 0) call AllocationError

            
            !
            ! Search for local element with common nodes
            !
            do ielem = 1,mesh%nelem

                ! Loop through bc nodes and see if current element has them all
                node_matched = .false.
                do inode = 1,nface_nodes
                    element_nodes = mesh%elems(ielem)%connectivity%get_element_nodes()
                    face_node     = bconnectivity%data(ielem_bc)%get_face_node(inode)
                    if ( any(element_nodes == face_node) ) then
                        node_matched(inode) = .true.
                    end if
                end do

        

                ! If all match, set element/face indices in boundary condition list.
                if ( all(node_matched) ) then
                    iface_bc = iface_bc + 1

                    !
                    ! Set element index
                    !
                    self%elems(iface_bc)   = ielem



                    !
                    ! Detect element face index associated with the boundary condition
                    !
                    ! Get xi,eta,zeta for three points, defining a face
                    !
                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(1))%c3_
                    point_one = mesh%elems(ielem)%computational_point(x,y,z)


                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(2))%c3_
                    point_two = mesh%elems(ielem)%computational_point(x,y,z)


                    x = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c1_
                    y = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c2_
                    z = mesh%nodes(bconnectivity%data(ielem_bc)%data(nface_nodes))%c3_
                    point_three = mesh%elems(ielem)%computational_point(x,y,z)


                    xi_face   = ( (abs(point_one%c1_ - point_two%c1_) < 1000.*RKTOL)  .and. (abs(point_one%c1_ - point_three%c1_) < 1000.*RKTOL) )
                    eta_face  = ( (abs(point_one%c2_ - point_two%c2_) < 1000.*RKTOL)  .and. (abs(point_one%c2_ - point_three%c2_) < 1000.*RKTOL) )
                    zeta_face = ( (abs(point_one%c3_ - point_two%c3_) < 1000.*RKTOL)  .and. (abs(point_one%c3_ - point_three%c3_) < 1000.*RKTOL) )


                    ! Determine iface by value of xi,eta,zeta
                    if ( xi_face ) then
                        if (point_one%c1_ < ZERO) then
                            iface = XI_MIN
                        else
                            iface = XI_MAX
                        end if
                    else if ( eta_face ) then
                        if (point_one%c2_ < ZERO) then
                            iface = ETA_MIN
                        else
                            iface = ETA_MAX
                        end if
                    else if ( zeta_face ) then
                        if (point_one%c3_ < ZERO) then
                            iface = ZETA_MIN
                        else
                            iface = ZETA_MAX
                        end if
                    else
                        call chidg_signal(FATAL,"bc%init: could not determine element face associated with the boundary")
                    end if


                    ! Set face index
                    self%faces(iface_bc) = iface






                    ! Set face type - 'ftype'
                    bcname = self%get_name()
                    if ( trim(bcname) == 'periodic' ) then
                        !
                        ! Set to ORPHAN face so it will be recognized as chimera in the detection process.
                        !
                        mesh%faces(ielem,iface)%ftype = ORPHAN

                        !
                        ! Set periodic offset from boundary condition to the face. To be used in detection of gq_donor.
                        !
                        if ( self%bcproperties%compute('type', time, pnt) == ONE ) then
                            mesh%faces(ielem,iface)%periodic_type = 'cartesian'
                        else if ( self%bcproperties%compute('type', time, pnt) == TWO ) then
                            mesh%faces(ielem,iface)%periodic_type = 'cylindrical'
                        end if

                        mesh%faces(ielem,iface)%chimera_offset_x     = self%bcproperties%compute('offset_x', time, pnt) ! time, pnt and do nothing here, but interface for function requires them.
                        mesh%faces(ielem,iface)%chimera_offset_y     = self%bcproperties%compute('offset_y', time, pnt)
                        mesh%faces(ielem,iface)%chimera_offset_z     = self%bcproperties%compute('offset_z', time, pnt)
                        mesh%faces(ielem,iface)%chimera_offset_theta = self%bcproperties%compute('offset_theta', time, pnt)

                    else
                        !
                        ! Set face to boundary condition face
                        !
                        mesh%faces(ielem,iface)%ftype = BOUNDARY
                    end if





                    ! End search
                    exit

                end if

            end do ! ielem





        end do !ibc_face





        !
        ! Call user-specialized boundary condition initialization
        !
        call self%init_spec(mesh)


        !
        ! Call user-specialized boundary coupling initialization
        !
        call self%init_boundary_coupling(mesh)



        self%isInitialized = .true. ! Set initialization confirmation

    end subroutine init
    !**********************************************************************************************
    








    !> Default specialized initialization procedure. This is called from the base bc%init procedure
    !! and can be overwritten by derived types to implement specialized initiailization details.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[inout]   mesh        mesh_t object containing elements and faces
    !!  @param[in]      iface       block face index to which the boundary condition is being applied
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_spec(self,mesh)
        class(bc_t),            intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh




    end subroutine init_spec
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
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_boundary_coupling(self,mesh)
        class(bc_t),    intent(inout)   :: self
        type(mesh_t),   intent(in)      :: mesh

        integer(ik) :: ielem_bc, ielem



        !
        ! Loop through elements and set default coupling information
        !
        do ielem_bc = 1,size(self%elems)


            !
            ! Get block-element index of current ielem_bc
            !
            ielem = self%elems(ielem_bc)

            
            !
            ! Add the element index as the only dependency.
            !
            call self%coupled_elems(ielem_bc)%push_back(ielem)


        end do ! ielem_bc


    end subroutine init_boundary_coupling
    !**********************************************************************************************









    !>  Function for returning the number of elements coupled with a specified boundary element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
    function get_ncoupled_elems(self,ielem) result(ncoupled_elems)
        class(bc_t),    intent(inout)   :: self
        integer(ik),    intent(in)      :: ielem

        integer(ik) :: ncoupled_elems


        ncoupled_elems = self%coupled_elems(ielem)%size()


    end function get_ncoupled_elems
    !**********************************************************************************************
















    !>  Apply boundary condition to the mesh and solution
    !!      - Loops through the associated elements(faces) and calls the specialized bc_t%compute
    !!        procedure for computing the rhs and linearization.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[in]      mesh    mesh_t defining elements and faces
    !!  @param[inout]   sdata   solverdata_t containing solution, rhs, and linearization(lin) data
    !!  @param[in]      iblk    Block of the linearization for the current element that is being computed (XI_MIN, XI_MAX, eta.)
    !!  @param[inout]   prop    properties_t object containing equationset properties and material_t objects
    !!
    !---------------------------------------------------------------------------------------------
    subroutine apply(self,mesh,sdata,prop,idomain_l)
        class(bc_t),            intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        class(solverdata_t),    intent(inout)   :: sdata
        class(properties_t),    intent(inout)   :: prop
        integer(ik),            intent(in)      :: idomain_l

        integer(ik) :: ielem_bc, ielement_l, ielement_c, iface, idonor, iflux, icoupled_elem, ncoupled_elems

        type(face_info_t)       :: face
        type(function_info_t)   :: fcn

        !
        ! Loop through associated boundary condition elements and call compute routine for the boundary flux calculation
        !
        do ielem_bc = 1,size(self%elems)
            ielement_l  = self%elems(ielem_bc)   ! Get index of the element being operated on
            iface       = self%faces(ielem_bc)   ! Get face index of element 'ielem' that is being operated on


            face%idomain_g  = mesh(idomain_l)%elems(ielement_l)%idomain_g
            face%idomain_l  = mesh(idomain_l)%elems(ielement_l)%idomain_l
            face%ielement_g = mesh(idomain_l)%elems(ielement_l)%ielement_g
            face%ielement_l = mesh(idomain_l)%elems(ielement_l)%ielement_l
            face%iface      = iface

            fcn%ifcn     = 0       ! Boundary conditions are not tracked.
            fcn%idepend  = 0       ! Chimera interface not applicable on boundary condition.
            fcn%idiff    = BC_BLK  ! Indicates to storage routine in LHS to store in BC section.

            
            !
            ! For current element, get number of coupled elements.
            !
            ncoupled_elems = self%get_ncoupled_elems(ielem_bc)


            !
            ! Compute current element function enough times to linearize all the coupled elements.
            ! If no coupling accross the face, the ncoupled_elems=1 for just the local interior element.
            !
            do icoupled_elem = 1,ncoupled_elems

                !
                ! Get coupled element to linearize against.
                !
                ielement_c = self%coupled_elems(ielem_bc)%at(icoupled_elem)
                !face%seed%idomain_g  = mesh(idomain_l)%elems(ielement_c)%idomain_g
                !face%seed%idomain_l  = mesh(idomain_l)%elems(ielement_c)%idomain_l
                !face%seed%ielement_g = mesh(idomain_l)%elems(ielement_c)%ielement_g
                !face%seed%ielement_l = mesh(idomain_l)%elems(ielement_c)%ielement_l
                !face%seed%iproc      = IRANK
                fcn%seed%idomain_g  = mesh(idomain_l)%elems(ielement_c)%idomain_g
                fcn%seed%idomain_l  = mesh(idomain_l)%elems(ielement_c)%idomain_l
                fcn%seed%ielement_g = mesh(idomain_l)%elems(ielement_c)%ielement_g
                fcn%seed%ielement_l = mesh(idomain_l)%elems(ielement_c)%ielement_l
                fcn%seed%iproc      = IRANK

                !
                ! For the current boundary element(face), call specialized compute procedure.
                !
                call self%compute(mesh,sdata,prop,face,fcn)

            end do !ielem_c


        end do !ielem_bc


    end subroutine apply
    !********************************************************************************************




















    !> Default options initialization procedure. This is called at the creation of a boundary condition
    !! in create_bc to set the options of a concrete bc_t. This function can be overwritten by a concrete
    !! bc_t to set case-specific options; parameters and functions.
    !!
    !!      - add entries to self%bcfunctions
    !!      - add entries to self%bcparameters
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(bc_t),            intent(inout)   :: self




    end subroutine add_options
    !********************************************************************************************










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
        class(bc_t),            intent(inout)   :: self
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
        class(bc_t),            intent(inout)   :: self
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
        class(bc_t),    intent(in)  :: self

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
        class(bc_t),    intent(in)  :: self
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
        class(bc_t),    intent(inout)   :: self
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
        class(bc_t),    intent(inout)  :: self
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
        class(bc_t),    intent(inout)  :: self
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
        class(bc_t),    intent(inout)   :: self
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
        class(bc_t),    intent(in)  :: self

        character(len=:), allocatable :: bcname

        bcname = self%name

    end function get_name
    !***************************************************************************************************







end module type_bc
