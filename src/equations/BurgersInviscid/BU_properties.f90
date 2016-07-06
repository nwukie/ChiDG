module BU_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t




    type, extends(properties_t), public :: BU_properties_t

        real(rk)    :: c(3)

    end type BU_properties_t


end module BU_properties
