module PRIMLINEULERAXI_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t
    implicit none




    type, extends(properties_t), public :: PRIMLINEULERAXI_properties_t

        real(rk)    :: R = 287.15_rk


    end type PRIMLINEULERAXI_properties_t








end module PRIMLINEULERAXI_properties
