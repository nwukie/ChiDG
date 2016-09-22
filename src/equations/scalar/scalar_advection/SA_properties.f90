module SA_properties
    use mod_kinds,  only: rk
    
    use type_properties,    only: properties_t




    type, extends(properties_t), public :: SA_properties_t

        real(rk)    :: c(3)

    end type SA_properties_t


end module
