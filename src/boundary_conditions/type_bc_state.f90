module type_bc_state
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: CARTESIAN, CYLINDRICAL, NO_ID
    use mod_chidg_mpi,          only: NRANK

    use type_bcproperty_set,    only: bcproperty_set_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_mesh,              only: mesh_t
    use type_bc_patch,          only: bc_patch_t
    use mpi_f08,                only: mpi_comm, mpi_integer, mpi_real8, mpi_integer4
    implicit none




    !>  Abstract base-type for computing a boundary condition state
    !!
    !!      - contains a procedure for computing the bc state: compute_bc_state
    !!      - compute_bc_state is deferred and so must be implemented by any new bc_state_t
    !!      - bc_state_t also contains properties that can hold parameters and functions 
    !!        that have been set for the boundary.
    !!
    !!  family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/9/2016
    !!  @note   Changed boundary condition to compute a bc state
    !!
    !-------------------------------------------------------------------------------------------
    type, public, abstract :: bc_state_t

        character(:),   allocatable :: name
        character(:),   allocatable :: family

        ! Boundary condition options
        type(bcproperty_set_t)      :: bcproperties

    contains

        procedure(bc_state_init),       deferred :: init
        procedure(bc_state_compute),    deferred :: compute_bc_state

        procedure   :: init_bc_precomm
        procedure   :: init_bc_postcomm
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

        ! Some predefined coupling strategies
        procedure   :: init_bc_coupling_global  ! All bc elements coupled with all other bc elements
        procedure   :: init_bc_coupling_local   ! Each bc element is coupled only with itself through the boundary


    end type bc_state_t
    !*******************************************************************************************




    abstract interface
        subroutine bc_state_init(self)
            import bc_state_t

            class(bc_state_t),  intent(inout)   :: self
        end subroutine
    end interface



    abstract interface
        subroutine bc_state_compute(self,worker,prop,bc_COMM)
            import bc_state_t
            import chidg_worker_t
            import properties_t
            import mpi_comm

            class(bc_state_t),      intent(inout)   :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
            type(mpi_comm),         intent(in)      :: bc_COMM
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
    subroutine init_bc_precomm(self,mesh,group_ID,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(inout)   :: mesh
        integer(ik),        intent(in)      :: group_ID
        type(mpi_comm),     intent(in)      :: bc_COMM



    end subroutine init_bc_precomm
    !**********************************************************************************************



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
    subroutine init_bc_postcomm(self,mesh,group_ID,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(inout)   :: mesh
        integer(ik),        intent(in)      :: group_ID
        type(mpi_comm),     intent(in)      :: bc_COMM



    end subroutine init_bc_postcomm
    !************************************************************************************






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
    !------------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(inout)   :: mesh
        integer(ik),        intent(in)      :: group_ID
        type(mpi_comm),     intent(in)      :: bc_COMM

        ! Default, initialize only local coupling.
        call self%init_bc_coupling_local(mesh,group_ID,bc_COMM)

    end subroutine init_bc_coupling
    !************************************************************************************




    !>  Set a function for a specified property.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcprop  String specifying a bcproperty_t to edit.
    !!  @param[in]  fcn     String specifying the concrete function_t to set.
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_fcn(self,bcprop,fcn)
        class(bc_state_t),      intent(inout)   :: self
        character(*),           intent(in)      :: bcprop
        character(*),           intent(in)      :: fcn

        call self%bcproperties%set_fcn(bcprop,fcn)

    end subroutine set_fcn
    !***********************************************************************************





    !>  Set a function option for a specified property.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcprop  String specifying a bcproperty_t to edit.
    !!  @param[in]  option  String specifying a particular option within bcproperty_f%fcn to edit
    !!  @param[in]  val     Real value to be set for the option.
    !!
    !------------------------------------------------------------------------------------
    subroutine set_fcn_option(self,bcprop,option,val)
        class(bc_state_t),      intent(inout)   :: self
        character(*),           intent(in)      :: bcprop
        character(*),           intent(in)      :: option
        real(rk),               intent(in)      :: val

        call self%bcproperties%set_fcn_option(bcprop,option,val)

    end subroutine set_fcn_option
    !************************************************************************************





    !>  Return number of properties available in the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_nproperties(self) result(nprop)
        class(bc_state_t),    intent(in)  :: self

        integer(ik) :: nprop

        nprop = self%bcproperties%get_nproperties()

    end function get_nproperties
    !********************************************************************************








    !>  Return a property name string, given the index of the property in the boundary 
    !!  condition.
    !!
    !!  This probably works best by first calling get_nproperties, and then iterating 
    !!  through the number of available properties to get their names.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop   Integer specifying the index of the property to be queried.
    !!  @result     pname   String of the property name associated with the index.
    !!
    !--------------------------------------------------------------------------------
    function get_property_name(self,iprop) result(pname)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop

        character(len=:),   allocatable :: pname

        pname = self%bcproperties%get_property_name(iprop) 

    end function get_property_name
    !********************************************************************************




    !>  Return an option key, given a property index and option index. 
    !!
    !!  One probably calls get_noptions(iprop) first, to get the number of available 
    !!  options for the function currently set for the property 'iprop'. Then one can 
    !!  loop over the number of available options and return their availble names 
    !!  dynamically.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop       Integer index of a property to modify.
    !!  @param[in]  ioption     Integer index of an option inside bcproperty%fcn
    !!  @result     key         String(key) corresponding to the option index (ioption)
    !!
    !--------------------------------------------------------------------------------
    function get_option_key(self,iprop,ioption) result(key)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop
        integer(ik),        intent(in)  :: ioption

        character(len=:),   allocatable :: key

        key = self%bcproperties%bcprop(iprop)%get_option_key(ioption)

    end function get_option_key
    !********************************************************************************




    !>  Return an option value, given a property index and option key.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop       Integer index of a property to modify.
    !!  @param[in]  key         String(key) specifying the option to be queried.
    !!  @result     val         Returned value of the selected key.
    !!
    !---------------------------------------------------------------------------------
    function get_option_value(self,iprop,key) result(val)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop
        character(*),       intent(in)  :: key

        real(rk)        :: val

        val = self%bcproperties%bcprop(iprop)%get_option_value(key)

    end function get_option_value
    !*********************************************************************************




    !>  Return the number of available options, given a property index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop       Integer index of a property to query.
    !!  @result     noption     Returned number of options available for the property. 
    !!                          Dependends on the function that is set for the property.
    !!
    !---------------------------------------------------------------------------------
    function get_noptions(self,iprop) result(noptions)
        class(bc_state_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iprop

        integer(ik)     :: noptions

        noptions = self%bcproperties%bcprop(iprop)%get_noptions()

    end function get_noptions
    !*********************************************************************************





    !>  Set the boundary condition name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !---------------------------------------------------------------------------------
    subroutine set_name(self,bcname)
        class(bc_state_t),  intent(inout)   :: self
        character(*),       intent(in)      :: bcname

        self%name = trim(bcname)

    end subroutine set_name
    !*********************************************************************************




    !>  Return the boundary condition name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    function get_name(self) result(bcname)
        class(bc_state_t),    intent(in)  :: self

        character(:),   allocatable :: bcname

        bcname = self%name

    end function get_name
    !**********************************************************************************



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
    !----------------------------------------------------------------------------------
    subroutine set_family(self,bc_family)
        class(bc_state_t),  intent(inout)   :: self
        character(*),       intent(in)      :: bc_family

        character(:),   allocatable :: user_msg

        ! Check incoming bc_state family
        if ( (trim(bc_family) == 'Inlet')           .or. &
             (trim(bc_family) == 'Outlet')          .or. &
             (trim(bc_family) == 'Wall')            .or. &
             (trim(bc_family) == 'Symmetry')        .or. &
             (trim(bc_family) == 'Periodic')        .or. &
             (trim(bc_family) == 'Farfield')        .or. &
             (trim(bc_family) == 'Scalar')          .or. &
             (trim(bc_family) == 'Mesh Motion')          .or. &
             (trim(bc_family) == 'Extrapolation')   .or. &
             (trim(bc_family) == 'Empty') ) then

            self%family = trim(bc_family)

        else

             user_msg = "bc_state%set_family: An invalid Family was trying to be set for the &
                         bc_state. Valid Families are: 'Inlet', 'Outlet', 'Wall', Symmetry', &
                         'Periodic', 'Farfield', 'Scalar', 'Mesh Motion', 'Extrapolation'."
             call chidg_signal_one(FATAL,user_msg,trim(bc_family))

        end if

    end subroutine set_family
    !*********************************************************************************





    !>  Return the boundary condition family.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/21/2016
    !!
    !---------------------------------------------------------------------------------
    function get_family(self) result(bc_family)
        class(bc_state_t),    intent(in)  :: self

        character(:),   allocatable :: bc_family, user_msg

        user_msg = "bc_state%get_family: It looks like the Family component for the &
                    bc_state was not set. Make sure self%set_family('my_family') is &
                    being called in the bc_state initialization procedure."
        if (.not. allocated(self%family)) call chidg_signal(FATAL,user_msg)

        bc_family = self%family

    end function get_family
    !**********************************************************************************





    !>  Initialize boundary group coupling.
    !!
    !!  Each element is coupled with every other element that belongs to the boundary
    !!  condition. This coupling occurs because each face uses an
    !!  average pressure that is computed over the group. The average pressure
    !!  calculation couples every element on the group. This coupling is initialized
    !!  here.
    !!
    !!  Coupling initialization:
    !!      1: each process loops through its local faces, initializes coupling
    !!         of all local faces with all other local faces.
    !!
    !!      2: loop through ranks in bc_COMM
    !!          a: iproc broadcasts information about its coupling to bc_COMM
    !!          b: all other procs receive from iproc and initialize parallel coupling
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling_global(self,mesh,group_ID,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(inout)   :: mesh
        integer(ik),        intent(in)      :: group_ID
        type(mpi_comm),     intent(in)      :: bc_COMM

        integer(ik) :: patch_ID, face_ID, elem_ID, patch_ID_coupled, face_ID_coupled,   &
                       idomain_g, idomain_l, ielement_g, ielement_l, iface,             &
                       bc_IRANK, bc_NRANK, ierr, iproc, nbc_elements,                   &
                       ielem, nfields, nterms_s, dof_start, dof_local_start, ngq, ibc

        integer(ik) :: idomain_g_coupled, idomain_l_coupled, ielement_g_coupled, ielement_l_coupled, &
                       iface_coupled, proc_coupled, send_size_a, send_size_b, send_size_c, send_size_d

        integer(ik) :: etype, nnodes, nterms_c, ntime, pelem_ID, interpolation_level,    &
                       coordinate_system, element_location(5), element_data(9), spacedim, inode

        real(rk),       allocatable :: interp_coords_def(:,:)
        real(rk),       allocatable :: areas(:)
        real(rk)                    :: total_area


        character(:),   allocatable :: coord_system
        real(rk),       allocatable :: nodes(:,:), nodes_def(:,:), nodes_vel(:,:), nodes_disp(:,:)
        integer(ik),    allocatable :: connectivity(:)



        !
        ! For each face, initialize coupling with all faces on the current processor.
        !
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            mesh%bc_patch_group(group_ID)%patch(patch_ID)%spatial_coupling  = 'Global'
            mesh%bc_patch_group(group_ID)%patch(patch_ID)%temporal_coupling = 'Global'
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                !
                ! Loop through, initialize coupling with all other patches/faces
                !
                do patch_ID_coupled = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID_coupled = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%nfaces()


                        !
                        ! Get block-element index of current face_ID_coupled
                        !
                        idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_g()
                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_l()
                        ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_g(face_ID_coupled)
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_l(face_ID_coupled)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%iface(     face_ID_coupled)


                        nfields           = mesh%domain(idomain_l)%elems(ielement_l)%nfields
                        nterms_s          = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s
                        dof_start         = mesh%domain(idomain_l)%elems(ielement_l)%dof_start
                        dof_local_start   = mesh%domain(idomain_l)%elems(ielement_l)%dof_local_start
                        total_area        = mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area
                        areas             = mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas
                        interp_coords_def = mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def



                        !
                        ! For the face (patch_ID,face_ID) add the element on (patch_ID_coupled,face_ID_coupled)
                        !
                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g,  &
                                                                                                        idomain_l,  &
                                                                                                        ielement_g, &
                                                                                                        ielement_l, &
                                                                                                        iface,      &
                                                                                                        IRANK)


                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g,       &
                                                                                                             ielement_g,      &
                                                                                                             nfields,         &
                                                                                                             ntime,           &
                                                                                                             nterms_s,        &
                                                                                                             dof_start,       &
                                                                                                             dof_local_start, &
                                                                                                             total_area,      &
                                                                                                             areas,           &
                                                                                                             interp_coords_def)

                    end do ! face_ID_couple
                end do ! patch_ID_couple

            end do ! face_ID
        end do ! patch_ID


        !
        ! Get bc_NRANK, bc_IRANK from bc_COMM
        !
        call MPI_Comm_Size(bc_COMM, bc_NRANK, ierr)
        call MPI_Comm_Rank(bc_COMM, bc_IRANK, ierr)


        !
        ! Initialize coupling with faces on other processors
        !
        do iproc = 0,bc_NRANK-1


            ! Send local elements out
            if (iproc == bc_IRANK) then


                nbc_elements = mesh%bc_patch_group(group_ID)%nfaces()
                call MPI_Bcast(IRANK,        1, MPI_INTEGER, iproc, bc_COMM, ierr)
                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)


                do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                        ! Broadcast element for coupling
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g(),         1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l(),         1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID), 1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID), 1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID),      1, MPI_INTEGER, iproc, bc_comm, ierr)

                        ! Broadcast auxiliary data
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%nfields,          1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,         1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%dof_start,        1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area, 1, MPI_INTEGER, iproc, bc_comm, ierr)

                        ngq = size(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def,1)
                        call MPI_Bcast(ngq,                                                                          1, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas,          ngq, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,1),      ngq, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,2),      ngq, MPI_INTEGER, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,3),      ngq, MPI_INTEGER, iproc, bc_comm, ierr)

                        ! Send information to construct mesh%parallel_element entry
                        send_size_a = size(mesh%domain(idomain_l)%elems(ielement_l)%connectivity)
                        send_size_b = size(mesh%domain(idomain_l)%elems(ielement_l)%node_coords)
                        send_size_c = size(mesh%domain(idomain_l)%elems(ielement_l)%node_coords_def)
                        send_size_d = size(mesh%domain(idomain_l)%elems(ielement_l)%node_coords_vel)
                
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%element_location,             5, mpi_integer4, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%element_data,                 9, mpi_integer4, iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%node_coords,        send_size_b, mpi_real8,    iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%node_coords_def,    send_size_c, mpi_real8,    iproc, bc_comm, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%elems(ielement_l)%node_coords_vel,    send_size_d, mpi_real8,    iproc, bc_comm, ierr)

                    end do ! face_ID
                end do ! patch_ID
            



            !
            ! All other processors recieve
            !
            else

                call MPI_Bcast(proc_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)

                ! For the face (patch_ID,face_ID) add each element from the sending proc
                do ielem = 1,nbc_elements

                    ! Receive coupled element
                    call MPI_BCast(idomain_g_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(idomain_l_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(ielement_g_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(ielement_l_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(iface_coupled,      1, MPI_INTEGER, iproc, bc_COMM, ierr)

                    ! Receive auxiliary data
                    call MPI_BCast(nfields,   1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(nterms_s,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(dof_start, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(total_area,1, MPI_INTEGER, iproc, bc_COMM, ierr)

                    call MPI_BCast(ngq, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    if (allocated(areas) ) deallocate(areas, interp_coords_def)
                    allocate(areas(ngq), interp_coords_def(ngq,3), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    call MPI_BCast(areas,                  ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,1), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,2), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,3), ngq, MPI_REAL8, iproc, bc_COMM, ierr)

                    ! Receive information to construct mesh%parallel_element entry
                    ! element_location = [idomain_g, idomain_l, ielement_g, ielement_l, iproc]
                    call mpi_bcast(element_location, 5, mpi_integer4, iproc, bc_comm, ierr)
                    idomain_g  = element_location(1)
                    ielement_g = element_location(3)

                    ! element_data = [element_type, spacedim, coordinate_system, nfields, nterms_s, nterms_c, ntime, interpolation_level]
                    call mpi_bcast(element_data, 9, mpi_integer4, iproc, bc_comm, ierr)
                    etype               = element_data(1)
                    spacedim            = element_data(2)
                    coordinate_system   = element_data(3)
                    nfields             = element_data(4)
                    nterms_s            = element_data(5)
                    nterms_c            = element_data(6)
                    ntime               = element_data(7)
                    interpolation_level = element_data(8)
                    dof_start           = element_data(9)
                    nnodes = (etype+1)*(etype+1)*(etype+1)

                    
                    ! Allocate buffers and receive: nodes, displacements, and velocities. 
                    ! These quantities are located at the element support nodes, not interpolation
                    ! nodes.
                    if (allocated(nodes)) deallocate(nodes, nodes_def, nodes_vel, connectivity)
                    allocate(nodes(       nnodes,3), &
                             nodes_def(   nnodes,3), &
                             nodes_vel(   nnodes,3), &
                             connectivity(nnodes  ), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    call mpi_bcast(nodes,     nnodes*3, mpi_real8, iproc, bc_comm, ierr)
                    call mpi_bcast(nodes_def, nnodes*3, mpi_real8, iproc, bc_comm, ierr)
                    call mpi_bcast(nodes_vel, nnodes*3, mpi_real8, iproc, bc_comm, ierr)

                    ! Compute node displacements
                    nodes_disp = nodes_def - nodes

                    ! Build local connectivity
                    !   : we construct the parallel element using just a local ordering
                    !   : so connectivity starts at 1 and goes to the number of nodes in 
                    !   : the element, nnodes.
                    !   :
                    !   :   connectivity = [1, 2, 3, 4, 5, 6, 7, 8 ...]
                    !   :
                    !   : We assume here that the displacements and velocities are ordered
                    !   : appropriately.
                    do inode = 1,nnodes
                        connectivity(inode) = inode
                    end do

                    ! Check for existing parallel element. If one does not
                    ! exist, get an identifier for a new parallel element.
                    pelem_ID = mesh%find_parallel_element(idomain_g,ielement_g)
                    if (pelem_ID == NO_ID) pelem_ID = mesh%new_parallel_element()

                    ! Initialize element geometry
                    select case(coordinate_system)
                        case(CARTESIAN)
                            coord_system = 'Cartesian'
                        case(CYLINDRICAL)
                            coord_system = 'Cylindrical'
                        case default
                            call chidg_signal(FATAL,"element%comm_recv: invalid coordinate system.")
                    end select

                    ! Construct/initialize/reinitialize parallel element
                    if (.not. mesh%parallel_element(pelem_ID)%geom_initialized) then
                        call mesh%parallel_element(pelem_ID)%init_geom(nodes,connectivity,etype,element_location,trim(coord_system))
                    end if

                    call mesh%parallel_element(pelem_ID)%init_sol('Quadrature',interpolation_level,nterms_s,nfields,ntime,dof_start,dof_local_start=NO_ID)
                    call mesh%parallel_element(pelem_ID)%set_displacements_velocities(nodes_disp,nodes_vel)
                    call mesh%parallel_element(pelem_ID)%update_interpolations_ale()




                    !
                    ! Each face on the current proc adds the off-processor element to their list 
                    ! of coupled elems
                    !
                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                        do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g_coupled,     &
                                                                                                            idomain_l_coupled,     &
                                                                                                            ielement_g_coupled,    &
                                                                                                            ielement_l_coupled,    &
                                                                                                            iface_coupled,         &
                                                                                                            proc_coupled)

                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g_coupled,     &
                                                                                                                 ielement_g_coupled,    &
                                                                                                                 nfields,               &
                                                                                                                 ntime,                 &
                                                                                                                 nterms_s,              &
                                                                                                                 dof_start,             &
                                                                                                                 NO_ID,                 &
                                                                                                                 total_area,            &
                                                                                                                 areas,                 &
                                                                                                                 interp_coords_def)

                        end do ! face_ID
                    end do ! patch_ID

                end do !ielem

            end if

            call MPI_Barrier(bc_COMM,ierr)
        end do

    end subroutine init_bc_coupling_global
    !******************************************************************************************



    !>  Default boundary coupling initialization routine. 
    !!
    !!  Default initializes coupling for a given element to just itself and no coupling with 
    !!  other elements on the boundary. For a boundary condition that is coupled across the face
    !!  this routine can be overwritten to set the coupling information specific to the boundary 
    !!  condition.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_coupling_local(self,mesh,group_ID,bc_COMM)
        class(bc_state_t),  intent(inout)   :: self
        type(mesh_t),       intent(inout)   :: mesh
        integer(ik),        intent(in)      :: group_ID
        type(mpi_comm),     intent(in)      :: bc_COMM

        integer(ik) :: patch_ID, face_ID, idomain_g, idomain_l,         &
                       ielement_g, ielement_l, iface, nfields, nterms_s,  &
                       dof_start



        !
        ! For each patch, loop through faces and set default element coupling.
        ! Default is that each face is coupled only with its owner element.
        ! So, strictly local coupling.
        !
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                !
                ! Get block-element index of current iface_bc
                !
                idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
                
                !
                ! Add the element index as the only dependency.
                !
                call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g,  &
                                                                                                idomain_l,  &
                                                                                                ielement_g, &
                                                                                                ielement_l, &
                                                                                                iface,      &
                                                                                                IRANK)


                call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID,    &
                                                                                            idomain_g,  &
                                                                                            ielement_g, &
                                                                                            mesh%domain(idomain_l)%elems(ielement_l)%nfields,                   &
                                                                                            mesh%domain(idomain_l)%elems(ielement_l)%ntime,                     &
                                                                                            mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,                  &
                                                                                            mesh%domain(idomain_l)%elems(ielement_l)%dof_start,                 &
                                                                                            mesh%domain(idomain_l)%elems(ielement_l)%dof_local_start,           &
                                                                                            mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area,          &
                                                                                            mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas,  &
                                                                                            mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def)

            end do ! face_ID
        end do ! patch_ID

    end subroutine init_bc_coupling_local
    !**********************************************************************************************

















end module type_bc_state
