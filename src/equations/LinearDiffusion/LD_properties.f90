module LD_properties
    use mod_kinds,  only: rk
    use type_properties,    only: properties_t




    !>  Properties for Linear Diffusion equation.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(properties_t), public :: LD_properties_t

        real(rk)    :: c(3)

    end type LD_properties_t
    !**************************************************************************


end module LD_properties
