! 
! rbf_info_t is a container for the information needed to
! initialize an rbf_interpolation_t.
! Information from the HDF5 gridfile is read and stored
! in an rbf_info_t instance, which is later used to 
! initialize an rbf_interpolation_t instance.
!

module type_rbf_info
    use mod_kinds,                  only: ik, rk
    implicit none

    type :: rbf_info_t
        character(len=1024)         :: rbf_name_default = "wc6"
        logical                     :: is_initialized   = .false.

        character(:), allocatable   :: rbf_name             ! If 'empty', use default
        character(:), allocatable   :: rbf_name_explicit    ! If 'empty', use default

        real(rk)                    :: radius(3)            ! If ZERO, use default

        real(rk)                    :: base_fraction    = 0.05_rk


    contains

        procedure                   :: init

    end type rbf_info_t

contains

    subroutine init(self, rbf_name, radius, rbf_name_explicit, base_fraction)
        class(rbf_info_t),                   intent(inout)       :: self
        character(*),                       intent(in)          :: rbf_name
        real(rk),                           intent(in)          :: radius(3)
        character(*),                       intent(in)          :: rbf_name_explicit
        real(rk),                           intent(in)          :: base_fraction


        self%rbf_name               = trim(rbf_name)
        self%radius                 = radius
        self%rbf_name_explicit      = trim(rbf_name_explicit)
        self%base_fraction          = base_fraction


        self%is_initialized         = .true.

    end subroutine init


end module type_rbf_info
