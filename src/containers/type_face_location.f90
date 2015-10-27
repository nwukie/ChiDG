module type_face_location
    use mod_kinds,  only: ik




    type, public :: face_location_t

        integer(ik) :: idomain
        integer(ik) :: ielement
        integer(ik) :: iface

    end type face_location_t





end module type_face_location
