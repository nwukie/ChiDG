module SD_properties
    use mod_kinds,  only: rk
    use type_properties,    only: properties_t




    !>  Properties for Linear Diffusion equation.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/16/2016
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(properties_t), public :: SD_properties_t

        real(rk)    :: mu(3)

    end type SD_properties_t
    !**************************************************************************


end module SD_properties
