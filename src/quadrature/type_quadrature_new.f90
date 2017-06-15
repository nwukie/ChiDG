module type_quadrature_new
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_legendre_gauss, only: lg_nodes, lg_weights
    implicit none






    !>  An object for returning quadrature nodes and weights.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/14/2017
    !!
    !--------------------------------------------------------------------
    type, public :: quadrature_new_t

    contains

        procedure   :: nodes
        procedure   :: weights

    end type quadrature_new_t
    !********************************************************************



contains




    !>  Return a 1D quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/14/2017
    !!
    !!  @param[in]  rule    String selecting the quadrature rule.
    !!  @param[in]  res     Number of nodes requested in the quadrature rule.
    !!
    !--------------------------------------------------------------------
    function nodes(self,rule,res) result(nodes_)
        class(quadrature_new_t),    intent(in)  :: self
        character(*),               intent(in)  :: rule
        integer(ik),                intent(in)  :: res

        real(rk)    :: nodes_(res)


        select case(trim(rule))
            case('Legendre-Gauss')
                nodes_ = lg_nodes(res)

            case default
                call chidg_signal_one(FATAL,"quadrature_new%nodes: We didn't find a quadrature rule that matched the selection.",trim(rule))

        end select

    end function nodes
    !*********************************************************************







    !>  Return a 1D quadrature weight set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/14/2017
    !!
    !!  @param[in]  rule    String selecting the quadrature rule.
    !!  @param[in]  res     Number of weights requested in the quadrature rule.
    !!
    !--------------------------------------------------------------------
    function weights(self,rule,res) result(weights_)
        class(quadrature_new_t),    intent(in)  :: self
        character(*),               intent(in)  :: rule
        integer(ik),                intent(in)  :: res

        real(rk)    :: weights_(res)


        select case(trim(rule))
            case('Legendre-Gauss')
                weights_ = lg_weights(res)

            case default
                call chidg_signal_one(FATAL,"quadrature_new%weights: We didn't find a quadrature rule that matched the selection.",trim(rule))

        end select

    end function weights
    !*********************************************************************












end module type_quadrature_new
