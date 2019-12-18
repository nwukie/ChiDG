module mod_test_functional_utilities
#include <messenger.h>
    
    use hdf5
    use h5lt

    use mod_hdf_utilities,         only: open_file_hdf, create_functional_group_hdf,                &
                                         create_functional_hdf, set_functional_reference_geom_hdf,  &
                                         set_functional_auxiliary_geom_hdf, set_functional_LS_hdf,  &
                                         close_file_hdf, remove_functional_hdf,                     &
                                         open_functional_group_hdf, close_functional_group_hdf
    
    
    implicit none

contains


    !>  Add a functional to a test mesh file.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/12/2017
    !!
    !!  @param[in]  filename            Name of the mesh to which we need to add the functional
    !!  @param[in]  functional_name     The name of the registered functional to be added
    !!  @param[in]  reference_geom      Name of the reference geometry
    !!  @param[in]  auxiliary_geom      Name of the auxiliary geometry
    !!  @param[in]  linear_solver       Name of the linear solver
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine meshfile_add_functional(filename,functional_name,reference_geom,auxiliary_geom,linear_solver)
        character(*),                           intent(in)  :: filename
        character(*),                           intent(in)  :: functional_name
        character(*),                           intent(in)  :: reference_geom
        character(*),    optional,              intent(in)  :: auxiliary_geom
        character(*),    optional,              intent(in)  :: linear_solver

        integer(HID_T)     :: fid

        ! Open hdf meshfile
        fid = open_file_hdf(filename)

        ! Create the functional group if necessary
        call create_functional_group_hdf(fid)

        ! Create functional
        call create_functional_hdf(fid,functional_name)

        ! Define reference geometry
        call set_functional_reference_geom_hdf(fid,functional_name,reference_geom)
    
        ! Define auxiliary geometry
        if (present(auxiliary_geom)) call set_functional_auxiliary_geom_hdf(fid,functional_name,auxiliary_geom)

        ! Define linear solver
        if (present(linear_solver)) call set_functional_LS_hdf(fid,functional_name,linear_solver)

        ! Close hdf file
        call close_file_hdf(fid)

    end subroutine meshfile_add_functional
    !---------------------------------------------------------------------------------------------





    !>  Edit a functional to a test mesh file.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/12/2017
    !!
    !!  @param[in]  filename            Name of the mesh to which we need to add the functional
    !!  @param[in]  functional_name     The name of the registered functional to be added
    !!  @param[in]  attribute           Attribute to be changed 
    !!  @param[in]  new_entry           New entry
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine meshfile_edit_functional(filename,functional_name,attribute,new_entry)
        character(*),                           intent(in)  :: filename
        character(*),                           intent(in)  :: functional_name
        character(*),                           intent(in)  :: attribute
        character(*),                           intent(in)  :: new_entry

        integer(HID_T)     :: fid
        character(:), allocatable   :: user_msg

        ! Open hdf meshfile
        fid = open_file_hdf(filename)

        select case (trim(attribute))
            case("reference geometry")
                
                ! Define reference geometry
                call set_functional_reference_geom_hdf(fid,functional_name,new_entry)
            
            case("auxiliary geometry")
                
                ! Define auxiliary geometry
                call set_functional_auxiliary_geom_hdf(fid,functional_name,new_entry)
            
            case("linear solver")
    
                ! Define linear solver
                call set_functional_LS_hdf(fid,functional_name,new_entry)

            case default
                user_msg = "edit_functional2meshfile: There was no valid case that matched the incoming string"
                call chidg_signal(FATAL,user_msg)
        
        end select
        
        ! Close hdf file
        call close_file_hdf(fid)


    end subroutine meshfile_edit_functional
    !---------------------------------------------------------------------------------------------






    !>  Remove a functional from a test mesh file.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/12/2017
    !!
    !!  @param[in]  filename            Name of the mesh to which we need to remove the functional
    !!  @param[in]  functional_name     The name of the registered functional to be removed
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine meshfile_remove_functional(filename,functional_name)
        character(*),                           intent(in)  :: filename
        character(*),                           intent(in)  :: functional_name

        integer(HID_T)     :: fid, fcl_id
        character(:), allocatable   :: user_msg

        ! Open hdf meshfile
        fid = open_file_hdf(filename)
        fcl_id = open_functional_group_hdf(fid)

        ! Remove functional
        call remove_functional_hdf(fcl_id,trim(functional_name))

        ! Close functional group
        call close_functional_group_hdf(fcl_id)
        ! Close hdf file
        call close_file_hdf(fid)

    end subroutine meshfile_remove_functional
    !---------------------------------------------------------------------------------------------


end module mod_test_functional_utilities
