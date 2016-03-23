module PRIMLINEULER_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t
    implicit none




    type, extends(properties_t), public :: PRIMLINEULER_properties_t

        real(rk)    :: R = 287.15_rk


    end type PRIMLINEULER_properties_t








end module PRIMLINEULER_properties
