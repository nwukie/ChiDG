module eqn_wall_distance
#include <messenger.h>
    use mod_constants,         only: ZERO, ONE, TWO, THREE
    use type_equation_set,     only: equation_set_t
    use type_equation_builder, only: equation_builder_t
    use type_scalar,           only: scalar_t
    use DNAD_D
    implicit none



    !>  A builder for an approximate Wall Distance equation set based on a p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(equation_builder_t), public :: wall_distance

    contains

        procedure   :: init
        procedure   :: build

    end type wall_distance
    !******************************************************************************************







    !>  A new diffusion coefficient model for the scalar equations. This facilitates
    !!  the implementation of a p-Poisson equation, which is like a normal
    !!  Poisson equation, but with a nonlinear diffusion coefficient.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(scalar_t), public :: wall_distance_model

    contains

        procedure   :: compute_mu

    end type wall_distance_model
    !******************************************************************************************



    !!
    !! parameter 'p' in the p-Poisson equation
    !!
    !!   p=2 :: linear Poisson equation
    !!   p>2 :: nonlinear p-Poisson equation
    !!
    !! Procedures:
    !!   set_p_poisson_parameter
    !!   get_p_poisson_parameter
    !!
    real(rk)    :: p = 2._rk



contains



    !-------------------------------------------------------------------------------
    !                      Scalar equation model mu
    !-------------------------------------------------------------------------------


    !>  Diffusion coefficient function for the p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    impure elemental function compute_mu(self,u,dudx,dudy,dudz) result(val)
        class(wall_distance_model), intent(in)  :: self
        type(AD_D),                 intent(in)  :: u
        type(AD_D),                 intent(in), optional    :: dudx
        type(AD_D),                 intent(in), optional    :: dudy
        type(AD_D),                 intent(in), optional    :: dudz

        type(AD_D)  :: val, mag2



        ! Compute magnitude of gradient
        mag2 = dudx*dudx + dudy*dudy + dudz*dudz

        
        if (abs(p-2._rk) > 1.e-14_rk) then
            val = mag2**((p-TWO)/TWO)
        else
            val = mag2
            val = ONE
        end if



    end function compute_mu
    !*****************************************************************************************







    !>  Set the 'p' parameter in the p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_p_poisson_parameter(p_in)
        real(rk),   intent(in)  :: p_in

        p = p_in

    end subroutine set_p_poisson_parameter
    !*****************************************************************************************



    !>  Get the 'p' parameter in the p-Poisson equation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_p_poisson_parameter() result(p_out)
        real(rk)    :: p_out

        p_out = p

    end function get_p_poisson_parameter
    !*****************************************************************************************








    !>  Wall Distance builder initialization.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_distance),   intent(inout)  :: self

        call self%set_name('Wall Distance')

    end subroutine init
    !*****************************************************************************************





    !>  Build a p-Poisson equation.
    !!
    !!  This substitutes the standard Laplace operator in the Poisson equation 
    !!  with a p-Laplace operator, creating a p-Poisson equation.
    !!
    !!  A nice property of the p-Poisson equation is that as 'p' goes to infinity,
    !!  the solution approaches the distance field, satisfying the Eikonal equation.
    !!  This allows one to obtain approximate distance fields, selecting the fidelity
    !!  of the approximation that is required.
    !!
    !!  The practical implementation of this is to create a new nonlinear diffusion 
    !!  coefficient model that will then be used in the regular scalar diffusion equation.
    !!  A source term is also implemented so that the equation is an approximation of the
    !!  distance field.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/16/2016
    !!
    !-----------------------------------------------------------------------------------------
    function build(self,blueprint) result(wall_distance_eqn)
        class(wall_distance),   intent(in)  :: self
        character(len=*),       intent(in)  :: blueprint

        character(:),   allocatable :: user_msg
        type(equation_set_t)        :: wall_distance_eqn
        type(wall_distance_model)   :: wall_distance_model_instance
        

        !
        ! Set equationset name.
        !
        call wall_distance_eqn%set_name("Wall Distance")


        !
        ! Add spatial operators
        !
        select case (trim(blueprint))

            case('default')
                call wall_distance_eqn%add_operator("Scalar Diffusion Boundary Average Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion Volume Operator")
                call wall_distance_eqn%add_operator("Scalar Diffusion BC Operator")
                call wall_distance_eqn%add_operator("Wall Distance Volume Source")


                call wall_distance_eqn%prop%add_scalar(wall_distance_model_instance)

            case default
                user_msg = "build_wall_distance: I didn't recognize the construction &
                            parameter that was passed to build the equation set."
                call chidg_signal_one(FATAL,user_msg, blueprint)

        end select


    end function build
    !*****************************************************************************************





end module eqn_wall_distance
