module EULER_properties
    use mod_kinds,          only: rk
    use type_properties,    only: properties_t
    implicit none




    !> Properties data type for the Euler equations
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !-----------------------------------------------------------------------
    type, extends(properties_t), public :: EULER_properties_t

        real(rk)    :: R = 287.15_rk


    end type EULER_properties_t
    !***********************************************************************








end module EULER_properties
