module eqn_wall_distance
#include <messenger.h>
    use mod_constants,         only: ZERO
    use type_equation_set,     only: equation_set_t
    use type_equation_builder, only: equation_builder_t
    use type_scalar,           only: scalar_t
    use DNAD_D
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: wall_distance

    contains

        procedure   :: init
        procedure   :: build

    end type wall_distance
    !******************************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(scalar_t), public :: wall_distance_model



    contains

        procedure   :: compute_mu
        procedure   :: compute_cx
        procedure   :: compute_cy
        procedure   :: compute_cz

    end type wall_distance_model
    !******************************************************************************************************




contains



    !-------------------------------------------------------------------------------
    !                      Scalar equation model for cx,cy,cz,mu
    !-------------------------------------------------------------------------------


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_mu(self,u) result(val)
        class(wall_distance_model), intent(in)  :: self
        type(AD_D),                 intent(in)  :: u

        type(AD_D)  :: val
        integer(ik) :: igq

        val = u*100.01_rk
        val = 1.0_rk

        !val%xp_ad_ = ZERO

    end function compute_mu
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cx(self,u,dudx,dudy,dudz) result(val)
        class(wall_distance_model), intent(in)  :: self
        type(AD_D),                 intent(in)  :: u
        type(AD_D),                 intent(in), optional    :: dudx
        type(AD_D),                 intent(in), optional    :: dudy
        type(AD_D),                 intent(in), optional    :: dudz

        type(AD_D)  :: val, mag

        ! Compute magnitude of gradient
        mag = sqrt(dudx*dudx + dudy*dudy + dudz*dudz)

        ! Set to constant
        val = dudx/mag

    end function compute_cx
    !*******************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cy(self,u,dudx,dudy,dudz) result(val)
        class(wall_distance_model), intent(in)  :: self
        type(AD_D),                 intent(in)  :: u
        type(AD_D),                 intent(in), optional    :: dudx
        type(AD_D),                 intent(in), optional    :: dudy
        type(AD_D),                 intent(in), optional    :: dudz

        type(AD_D) :: val, mag

        ! Compute magnitude of gradient
        mag = sqrt(dudx*dudx + dudy*dudy + dudz*dudz)

        ! Set to constant
        val = dudy/mag

    end function compute_cy
    !*******************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    impure elemental function compute_cz(self,u,dudx,dudy,dudz) result(val)
        class(wall_distance_model), intent(in)  :: self
        type(AD_D),                 intent(in)  :: u
        type(AD_D),                 intent(in), optional    :: dudx
        type(AD_D),                 intent(in), optional    :: dudy
        type(AD_D),                 intent(in), optional    :: dudz

        type(AD_D)  :: val, mag

        ! Compute magnitude of gradient
        mag = sqrt(dudx*dudx + dudy*dudy + dudz*dudz)

        ! Set to constant
        val = dudz/mag

    end function compute_cz
    !*******************************************************************************































    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_distance),   intent(inout)  :: self

        call self%set_name('Wall Distance')

    end subroutine init
    !*********************************************************************************************





    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    function build(self,blueprint) result(wall_distance_eqn)
        class(wall_distance),    intent(in)  :: self
        character(len=*),        intent(in)  :: blueprint

        type(equation_set_t)         :: wall_distance_eqn
        type(wall_distance_model)    :: wall_distance_model_instance
        

        !
        ! Set equationset name.
        !
        call wall_distance_eqn%set_name("Wall Distance")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call wall_distance_eqn%add_operator('Scalar Advection Boundary Average Operator')
                call wall_distance_eqn%add_operator('Scalar Advection Volume Operator')
                call wall_distance_eqn%add_operator('Scalar Advection LaxFriedrichs Operator')
                call wall_distance_eqn%add_operator('Scalar Advection BC Operator')
                call wall_distance_eqn%add_operator("Scalar Diffusion Boundary Average Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion Volume Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion BC Operator")
                call wall_distance_eqn%add_operator("Wall Distance Volume Source")


                call wall_distance_eqn%prop%add_scalar(wall_distance_model_instance)

            case default
                call chidg_signal_one(FATAL, "build_wall_distance: I didn't recognize the construction &
                                              parameter that was passed to build the equation set.", blueprint)

        end select


    end function build
    !*********************************************************************************************************





end module eqn_wall_distance
