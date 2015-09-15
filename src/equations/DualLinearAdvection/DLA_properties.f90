module DLA_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t




    type, extends(properties_t), public :: DLA_properties_t

        real(rk)    :: c(3)

    end type DLA_properties_t


end module DLA_properties
