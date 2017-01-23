module type_chidg_discretization
    use type_chidg_data,    only: chidg_data_t
    implicit none






    !>  Handles computing the spatio-temporal contributions to the right-hand side and
    !!  left-hand side data containers (rhs,lhs).
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2017
    !!
    !---------------------------------------------------------------------------------------
    type, public :: chidg_discretization_t



    contains

        procedure   :: update

    end type chidg_discretization_t
    !***************************************************************************************






contains




    !>  Update contributions from spatial and temporal residual terms.
    !!
    !!      For both temporal and spatial derivatives:
    !!          - Contribute to rhs
    !!          - Optionally, contribute to lhs for implicit problems
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine update(self,data)
        class(chidg_discretization_t),  intent(in)      :: self
        type(chidg_data_t),             intent(inout)   :: data









    end subroutine update
    !***************************************************************************************



end module type_chidg_discretization
