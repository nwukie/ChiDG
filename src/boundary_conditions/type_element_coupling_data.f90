module type_element_coupling_data
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NO_PROC, NO_ID, ZERO
    use type_point,     only: point_t
    implicit none


    !>  Data container for information regarding element coupling.
    !!
    !!  This would be used on a boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2017
    !!
    !---------------------------------------------------------------------
    type, public :: element_coupling_data_t

        ! Element coupling
        integer(ik) :: idomain_g    = NO_ID
        integer(ik) :: idomain_l    = NO_ID
        integer(ik) :: ielement_g   = NO_ID
        integer(ik) :: ielement_l   = NO_ID
        integer(ik) :: iface        = NO_ID
        integer(ik) :: proc         = NO_PROC


        ! Element parallel recv access
        integer(ik) :: recv_comm    = NO_ID
        integer(ik) :: recv_domain  = NO_ID
        integer(ik) :: recv_element = NO_ID
        integer(ik) :: recv_dof     = NO_ID

        
        ! Element data
        integer(ik)                 :: coordinate_system
        integer(ik)                 :: nfields     = 0
        integer(ik)                 :: ntime       = 0
        integer(ik)                 :: nterms_s    = 0
        integer(ik)                 :: nnodes_r    = 0
        integer(ik)                 :: dof_start   = 0
        integer(ik)                 :: dof_local_start   = 0
        real(rk)                    :: total_area  = ZERO
        real(rk),       allocatable :: areas(:)
        type(point_t),  allocatable :: quad_pts(:)

    contains

        procedure   :: set_coupling
        procedure   :: set_recv
        procedure   :: set_data

    end type element_coupling_data_t
    !*********************************************************************


contains



    !>  Set element coupling data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2017
    !!
    !!
    !---------------------------------------------------------------------
    subroutine set_coupling(self,idomain_g,idomain_l,ielement_g,ielement_l,iface,proc)
        class(element_coupling_data_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: idomain_g
        integer(ik),                    intent(in)      :: idomain_l
        integer(ik),                    intent(in)      :: ielement_g
        integer(ik),                    intent(in)      :: ielement_l
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(in)      :: proc


        self%idomain_g  = idomain_g
        self%idomain_l  = idomain_l
        self%ielement_g = ielement_g
        self%ielement_l = ielement_l
        self%iface      = iface
        self%proc       = proc

    end subroutine set_coupling
    !**********************************************************************




    !>  Set receive access indices for accessing parallel data 
    !!  in chidg_vector.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2017
    !!
    !!
    !----------------------------------------------------------------------
    subroutine set_recv(self,recv_comm,recv_domain,recv_element,recv_dof)
        class(element_coupling_data_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: recv_comm
        integer(ik),                    intent(in)      :: recv_domain
        integer(ik),                    intent(in)      :: recv_element
        integer(ik),                    intent(in)      :: recv_dof

        self%recv_comm    = recv_comm
        self%recv_domain  = recv_domain
        self%recv_element = recv_element
        self%recv_dof     = recv_dof

    end subroutine set_recv
    !**********************************************************************






    !>  Set auxiliary data for coupled element.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2017
    !!
    !----------------------------------------------------------------------
    subroutine set_data(self,nfields,ntime,nterms_s,nnodes_r,coordinate_system,dof_start,dof_local_start,total_area,areas,quad_pts)
        class(element_coupling_data_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: nfields
        integer(ik),                    intent(in)      :: ntime
        integer(ik),                    intent(in)      :: nterms_s
        integer(ik),                    intent(in)      :: nnodes_r
        integer(ik),                    intent(in)      :: coordinate_system
        integer(ik),                    intent(in)      :: dof_start
        integer(ik),                    intent(in)      :: dof_local_start
        real(rk),                       intent(in)      :: total_area
        real(rk),                       intent(in)      :: areas(:)
        type(point_t),                  intent(in)      :: quad_pts(:)

        self%nfields           = nfields
        self%ntime             = ntime
        self%nterms_s          = nterms_s
        self%nnodes_r          = nnodes_r
        self%coordinate_system = coordinate_system
        self%dof_start         = dof_start
        self%dof_local_start   = dof_local_start
        self%total_area        = total_area
        self%areas             = areas
        self%quad_pts          = quad_pts

    end subroutine set_data
    !**********************************************************************








end module type_element_coupling_data
