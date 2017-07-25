module type_chimera_send
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    use type_element,   only: element_t
    use mod_chidg_mpi,  only: IRANK
    implicit none



    !>  Container for chimera donor elements that need to be sent off-processor.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !-----------------------------------------------------------------------------------
    type, public :: chimera_send_t

        ! Local donor element
        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l

        ! Processors element is being sent to
        type(ivector_t) :: send_procs

    contains

        procedure   :: nsend_procs

    end type chimera_send_t
    !***********************************************************************************

    interface chimera_send
        module procedure chimera_send
    end interface





contains


    !>  Constructor for chimera_send_t
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !----------------------------------------------------------------------------------
    function chimera_send(idomain_g, idomain_l, ielement_g, ielement_l) result(instance)
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_g
        integer(ik),    intent(in)  :: ielement_l

        type(chimera_send_t)    :: instance

        instance%idomain_g  = idomain_g
        instance%idomain_l  = idomain_l
        instance%ielement_g = ielement_g
        instance%ielement_l = ielement_l

    end function chimera_send
    !***********************************************************************************




    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------
    function nsend_procs(self) result(nsend_)
        class(chimera_send_t), intent(in)  :: self

        integer(ik) :: nsend_

        nsend_ = self%send_procs%size()

    end function nsend_procs
    !***********************************************************************************




end module type_chimera_send
