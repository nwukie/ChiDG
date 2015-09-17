module EULER_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t




    type, extends(properties_t), public :: EULER_properties_t

        real(rk)    :: R = 287.15_rk


    end type EULER_properties_t








end module
