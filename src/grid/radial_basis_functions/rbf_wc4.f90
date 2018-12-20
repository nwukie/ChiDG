module rbf_wc4
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_radial_basis_function,  only: radial_basis_function_t
    use mod_rbf_tools
    implicit none








    !>  Wendland's c4 compactly supported RBF 
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(radial_basis_function_t), public :: wc4_rbf

        
    contains

        procedure   :: init
        procedure   :: compute
        procedure   :: compute_grad

    end type wc4_rbf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(wc4_rbf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("wc4_rbf")




    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute(self,eval_node,support_node,support_radius) result(val)
        class(wc4_rbf),     intent(inout)   :: self
        real(rk),               intent(in)      :: eval_node(3)
        real(rk),               intent(in)      :: support_node(3)
        real(rk),               intent(in)      :: support_radius(3)

        real(rk)                                :: val

        real(rk) :: radius_norm
        
        !2D case ?
        !radius_norm = sqrt((val(1)-support_val(1))**TWO+(val(2)-support_val(2))**TWO)


        !3D Case - how to handle dimensionality more generally?
        radius_norm = compute_radius_norm(eval_node, support_node, support_radius)

        if (radius_norm < ONE) then

            val = (ONE-radius_norm)**6.0_rk*(35.0_rk*radius_norm**TWO+18.0_rk*radius_norm+3.0_rk)
        else
            val = ZERO
        end if


    end function compute
    !**********************************************************************************

    !>
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_grad(self,eval_node,support_node,support_radius) result(val)
        class(wc4_rbf),     intent(inout)   :: self
        real(rk),               intent(in)      :: eval_node(3)
        real(rk),               intent(in)      :: support_node(3)
        real(rk),               intent(in)      :: support_radius(3)

        real(rk)                                :: val(3)

        real(rk) :: radius_norm
        real(rk) :: radius_norm_grad(3)
        
        !2D case ?
        !radius_norm = sqrt((val(1)-support_val(1))**TWO+(val(2)-support_val(2))**TWO)


        !3D Case - how to handle dimensionality more generally?
        radius_norm         = compute_radius_norm(eval_node, support_node, support_radius)
        radius_norm_grad    = compute_radius_norm_grad(eval_node, support_node, support_radius)

        if (radius_norm < ONE) then

            val = (6.0_rk*(ONE-radius_norm)**5.0_rk*(35.0_rk*radius_norm**TWO+18.0_rk*radius_norm+3.0_rk) + &
                    (ONE-radius_norm)**6.0_rk*(70.0_rk*radius_norm+18.0_rk))*radius_norm_grad

        else
            val = ZERO
        end if


    end function compute_grad
    !**********************************************************************************


end module rbf_wc4
