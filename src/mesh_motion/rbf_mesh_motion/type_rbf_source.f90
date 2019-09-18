module type_rbf_source
    use type_rbf_mm_driver,                 only: rbf_mm_driver_t
    implicit none

    !>
    !! Description: rbf_source_t contains a patch_name, indicating a set of boundary nodes,
    !!              and a RBF MM driver object, used to update displacements and velocities
    !!              of the boundary nodes.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    type, public    :: rbf_source_t

        character(:), allocatable         :: patch_name
        class(rbf_mm_driver_t), allocatable     :: driver

    contains

        procedure :: init
        procedure :: clear

    end type rbf_source_t

contains

    subroutine init(self, patch_name_in, driver_in)
        class(rbf_source_t),        intent(inout)   :: self
        character(*),           intent(in)      :: patch_name_in
        class(rbf_mm_driver_t),     intent(in)      :: driver_in

        call self%clear()

        self%patch_name = patch_name_in
        allocate(self%driver, source = driver_in)

    end subroutine init

    subroutine clear(self)
        class(rbf_source_t),        intent(inout)   :: self

        if (allocated(self%patch_name)) deallocate(self%patch_name)
        if (allocated(self%driver)) deallocate(self%driver)

    end subroutine clear



end module type_rbf_source
