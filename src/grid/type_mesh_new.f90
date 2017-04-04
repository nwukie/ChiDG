module type_mesh_new
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_chidg_mpi,              only: IRANK, NRANK, GLOBAL_MASTER
    use mpi_f08

    use type_domain,                only: domain_t
    use type_ivector,               only: ivector_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_element_connectivity,  only: element_connectivity_t
    implicit none
    private


    !> Data type for mesh information
    !!      - contains array of elements, array of faces for each element
    !!      - calls initialization procedure for elements and faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/27/2016
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: mesh_new_t

        type(domain_t), allocatable :: dom(:)

    contains


    end type mesh_new_t
    !*****************************************************************************************





contains





end module type_mesh_new
