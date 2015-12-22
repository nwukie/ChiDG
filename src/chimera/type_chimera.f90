module type_chimera
    use mod_kinds,              only: rk, ik
    use type_chimera_receiver,  only: chimera_receiver_t
    use type_chimera_donor,     only: chimera_donor_t



    !> Main interface and container for Chimera data and operations.
    !! Holds chimera send/receive sets which are used to facilitate inter-domain communication
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------------
    type, public :: chimera_t
    
        type(chimera_receiver_t)    :: recv
        type(chimera_donor_t)       :: send

    contains


    end type chimera_t
    !*************************************************************************************************


contains









end module type_chimera
