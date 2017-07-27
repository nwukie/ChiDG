module rbf_c0
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_radial_basis_function,  only: radial_basis_function_t
    implicit none
    private








    !>  Static mesh (no mesh motion).
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(radial_basis_function_t), public :: c0_rbf
        private

        
    contains

        procedure   :: init
        procedure   :: compute

    end type c0_rbf
    !********************************************************************************



contains




    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(c0_rbf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("c0_rbf")




    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date  3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute(self,node,support_node,support_radius) result(val)
        class(static_pmmf),     intent(inout)   :: self
        real(rk),               intent(in)      :: node(3)
        real(rk),               intent(in)      :: support_node(3)
        real(rk),               intent(in)      :: support_radius

        real(rk)                                :: val(3)

        real(rk) :: radius_norm
        
        !2D case ?
        !radius_norm = sqrt((val(1)-support_val(1))**TWO+(val(2)-support_val(2))**TWO)


        !3D Case - how to handle dimensionality more generally?
        radius_norm = sqrt((val(1)-support_val(1))**TWO+(val(2)-support_val(2))**TWO+(val(3)-support_val(3))**TWO)

        val = (ONE-radius_norm)**TWO


    end function compute
    !**********************************************************************************

end module rbf_c0
