!>  Explicit interfaces for calling TIOGA procedures.
!!
!!  These are already exposed by TIOGA and can be called 
!!  without these interfaces, but that does not admit any
!!  compile-time checking of errors. Calling through explicit
!!  these interfaces is a safer approach.
!!
!!
!!  @author Nathan A. Wukie (AFRL)
!!  @date   8/27/2019
!!
!!
!!  Subset of TIOGA procedures interfaced here:
!!  -------------------------------------------
!!  tioga_init_f90
!!  tioga_registergrid_data_mb
!!  tioga_preprocess_grids
!!  tioga_performconnectivity
!!  tioga_setnfringe
!!  tioga_reduce_fringes
!!
!!  More procedures exist in the TIOGA file tiogaInterface.C
!!
!----------------------------------------------------------------
module mod_tioga_interfaces
    use mpi_f08,        only: mpi_comm
    use iso_c_binding,  only: c_int
    implicit none


interface
    subroutine tioga_init_f90_(int_comm) bind(C,name="tioga_init_f90_")
        use iso_c_binding,  only: c_int
        implicit none
        integer(c_int), intent(in)  :: int_comm
    end subroutine tioga_init_f90_
end interface

interface 
    subroutine tioga_registergrid_data_mb(bid, btag, nnodes, xyz, iblank, nwall_nodes, noverset_nodes, wall_nodes, overset_nodes, nelement_types, nvertices_per_element, nelements, connectivity) bind(C,name="tioga_registergrid_data_mb_")
        use iso_c_binding,  only: c_int, c_double, c_ptr
        implicit none
        integer(c_int), intent(in), target  :: bid
        integer(c_int), intent(in), target  :: btag
        integer(c_int), intent(in), target  :: nnodes
        type(c_ptr),    intent(in), value   :: xyz       ! [3*nnodes]
        type(c_ptr),    intent(in), value   :: iblank    ! [nnodes?]
        integer(c_int), intent(in), target  :: nwall_nodes
        integer(c_int), intent(in), target  :: noverset_nodes
        type(c_ptr),    intent(in), value   :: wall_nodes
        type(c_ptr),    intent(in), value   :: overset_nodes
        integer(c_int), intent(in), target  :: nelement_types
        integer(c_int), intent(in), target  :: nvertices_per_element
        integer(c_int), intent(in), target  :: nelements
        type(c_ptr),    intent(in), value   :: connectivity 
    end subroutine tioga_registergrid_data_mb
end interface

interface
    subroutine tioga_preprocess_grids() bind(C,name="tioga_preprocess_grids_")
        implicit none
    end subroutine tioga_preprocess_grids
end interface

interface
    subroutine tioga_performconnectivity() bind(C,name="tioga_performconnectivity_")
        implicit none
    end subroutine tioga_performconnectivity
end interface

interface
    subroutine tioga_setnfringe(nfringe) bind(C,name="tioga_setnfringe_")
        use iso_c_binding,  only: c_int
        implicit none
        integer(c_int), intent(in), target :: nfringe
    end subroutine tioga_setnfringe
end interface

interface
    subroutine tioga_reduce_fringes() bind(C,name="tioga_reduce_fringes_")
        implicit none
    end subroutine tioga_reduce_fringes
end interface

contains

    !>  Pass mpi_comm datatype instead of integer. Manage getting integer 
    !!  identifier and passing to TIOGA.
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
