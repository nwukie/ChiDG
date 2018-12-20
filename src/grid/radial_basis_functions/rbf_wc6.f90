module rbf_wc6
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_radial_basis_function,  only: radial_basis_function_t
    use mod_rbf_tools
    implicit none








    !>  Wendland's c6 compactly supported RBF 
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(radial_basis_function_t), public :: wc6_rbf

        
    contains

        procedure   :: init
        procedure   :: compute
        procedure   :: compute_grad

    end type wc6_rbf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date  8/4/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(wc6_rbf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("wc6_rbf")




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
        class(wc6_rbf),     intent(inout)   :: self
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

            val = (ONE-radius_norm)**8.0_rk*(32.0_rk*radius_norm**THREE+25.0_rk*radius_norm**TWO+8.0_rk*radius_norm+1.0_rk)
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
        class(wc6_rbf),     intent(inout)   :: self
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

            val = (8.0_rk*(ONE-radius_norm)**7.0_rk*(32.0_rk*radius_norm**THREE+25.0_rk*radius_norm**TWO+8.0_rk*radius_norm+1.0_rk) +&
                    (ONE-radius_norm)**8.0_rk*(96.0_rk*radius_norm**TWO+50.0_rk*radius_norm+8.0_rk)) *&
                    radius_norm_grad;

        else
            val = ZERO
        end if


    end function compute_grad
    !**********************************************************************************


end module rbf_wc6
