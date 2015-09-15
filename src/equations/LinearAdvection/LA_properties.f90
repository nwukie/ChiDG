module LA_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t




    type, extends(properties_t), public :: LA_properties_t

        real(rk)    :: c(3)

    end type LA_properties_t


end module
