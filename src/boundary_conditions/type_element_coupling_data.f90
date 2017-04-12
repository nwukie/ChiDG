module type_element_coupling_data
#include <messenger.h>
    use mod_kinds,  only: rk, ik
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

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l
        integer(ik) :: proc

        integer(ik) :: recv_comm
        integer(ik) :: recv_domain
        integer(ik) :: recv_element

    contains

        procedure   :: set_coupling
        procedure   :: set_recv

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
    subroutine set_coupling(self,idomain_g,idomain_l,ielement_g,ielement_l,proc)
        class(element_coupling_data_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: idomain_g
        integer(ik),                    intent(in)      :: idomain_l
        integer(ik),                    intent(in)      :: ielement_g
        integer(ik),                    intent(in)      :: ielement_l
        integer(ik),                    intent(in)      :: proc


        self%idomain_g  = idomain_g
        self%idomain_l  = idomain_l
        self%ielement_g = ielement_g
        self%ielement_l = ielement_l
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
    subroutine set_recv(self,recv_comm,recv_domain,recv_element)
        class(element_coupling_data_t), intent(inout)   :: self
        integer(ik),                    intent(in)      :: recv_comm
        integer(ik),                    intent(in)      :: recv_domain
        integer(ik),                    intent(in)      :: recv_element

        self%recv_comm    = recv_comm
        self%recv_domain  = recv_domain
        self%recv_element = recv_element

    end subroutine set_recv
    !**********************************************************************






end module type_element_coupling_data
