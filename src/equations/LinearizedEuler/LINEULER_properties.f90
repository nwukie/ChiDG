module LINEULER_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t
    implicit none




    type, extends(properties_t), public :: LINEULER_properties_t

        real(rk)    :: R = 287.15_rk


    end type LINEULER_properties_t








end module LINEULER_properties
