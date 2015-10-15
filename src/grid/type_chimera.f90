module type_chimera
    use mod_kinds,              only: rk, ik
    use type_chimera_receiver,  only: chimera_receiver_t
    use type_chimera_donor,     only: chimera_donor_t




    type, public :: chimera_t
    
        type(chimera_receiver_t),   allocatable :: recv(:)
        type(chimera_donor_t),      allocatable :: send(:)


    contains


    end type chimera_t









end module type_chimera
