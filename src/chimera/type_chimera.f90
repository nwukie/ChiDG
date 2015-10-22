module type_chimera
    use mod_kinds,              only: rk, ik
    use type_chimera_receiver,  only: chimera_receiver_t
    use type_chimera_donor,     only: chimera_donor_t



    !> Main interface and container for Chimera data and operations.
    !! Holds chimera send/receive sets which are used to facilitate inter-domain communication
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chimera_t
    
        type(chimera_receiver_t),   allocatable :: recv(:)
        type(chimera_donor_t),      allocatable :: send(:)


    contains

        !> Set up Chimera data
        procedure   :: init

    end type chimera_t


contains







    subroutine init(self)
        class(chimera_t),   intent(inout)   :: self




    end subroutine




end module type_chimera
