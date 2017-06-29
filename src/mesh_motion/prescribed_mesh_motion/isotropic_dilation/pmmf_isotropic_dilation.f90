module pmmf_isotropic_dilation
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_prescribed_mesh_motion_function,  only: prescribed_mesh_motion_function_t
    implicit none
    private








    !>  isotropic_dilation mesh motion. 
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(prescribed_mesh_motion_function_t), public :: isotropic_dilation_pmmf
        private

        
    contains

        procedure   :: init
        procedure   :: compute_pos
        procedure   :: compute_vel

    end type isotropic_dilation_pmmf
    !********************************************************************************



contains




    !>  X = (isotropic_dilation_factor*t+1)*x
    !!  v = isotropic_dilation_factor*x
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(isotropic_dilation_pmmf),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("isotropic_dilation")


        !
        ! Set function options to default settings
        !
        call self%add_option('isotropic_dilation_factor',ZERO)
        
    end subroutine init
    !*************************************************************************



    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_pos(self,time,node) result(val)
        class(isotropic_dilation_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)

        real(rk)                                    :: factor

        factor = self%get_option_value('isotropic_dilation_factor')

        val = (factor*time+ONE)*node 
       
        
    end function compute_pos
    !**********************************************************************************






    !>
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017 
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(isotropic_dilation_pmmf),     intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)

        integer(ik)                                 :: ivar
        real(rk)                                    :: val(3)
        real(rk)                                    :: factor

        factor = self%get_option_value('isotropic_dilation_factor')

        val = factor*node 

        
    end function compute_vel
    !**********************************************************************************


end module pmmf_isotropic_dilation
