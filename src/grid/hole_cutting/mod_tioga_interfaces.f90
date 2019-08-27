module mod_tioga_interfaces
    use mpi_f08,        only: mpi_comm
    use iso_c_binding,  only: c_int
    implicit none



interface
    subroutine tioga_init_f90_(int_comm) bind(C)
        use iso_c_binding,  only: c_int
        implicit none
        integer(c_int), intent(in)  :: int_comm
    end subroutine
end interface

interface 
    subroutine tioga_registergrid_data_mb(bid, btag, nnodes, xyz, iblank, nwall_nodes, noverset_nodes, wall_nodes, overset_nodes, nelement_types, nvertices_per_element, nelements, connectivity) bind(C,name="tioga_registergrid_data_mb_")
        use iso_c_binding,  only: c_int, c_double
        implicit none
        integer(c_int), intent(in)  :: bid
        integer(c_int), intent(in)  :: btag
        integer(c_int), intent(in)  :: nnodes
        real(c_double), intent(in)  :: xyz(:)       ! [3*nnodes]
        integer(c_int), intent(in)  :: iblank(:)    ! [nnodes?]
        integer(c_int), intent(in)  :: nwall_nodes
        integer(c_int), intent(in)  :: noverset_nodes
        integer(c_int), intent(in)  :: wall_nodes(:)
        integer(c_int), intent(in)  :: overset_nodes(:)
        integer(c_int), intent(in)  :: nelement_types
        integer(c_int), intent(in)  :: nvertices_per_element
        integer(c_int), intent(in)  :: nelements
        integer(c_int), intent(in)  :: connectivity(:,:)    ! [conn_size, nelements]
    end subroutine
end interface

interface
    subroutine tioga_preprocess_grids() bind(C,name="tioga_preprocess_grids_")
        implicit none
    end subroutine
end interface

interface
    subroutine tioga_performconnectivity() bind(C,name="tioga_performconnectivity_")
        implicit none
    end subroutine
end interface

contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/27/2019
    !!
    !----------------------------------------------------------------------------
    subroutine tioga_init_f90(comm)
        type(mpi_comm), intent(in)  :: comm

        call tioga_init_f90_(int(comm%mpi_val,c_int))

    end subroutine tioga_init_f90
    !****************************************************************************









end module mod_tioga_interfaces
