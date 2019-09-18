module mod_rbf_mm_driver
#include <messenger.h>
    use mod_kinds,                              only: rk, ik
    use type_rbf_mm_driver,   only: rbf_mm_driver_t
    use type_rbf_mm_driver_vector,                        only: rbf_mm_driver_vector_t

    !
    ! Import rbf_mm_drivers
    !
    use pmm_rbf_mm_driver,                      only: rbf_mm_driver_pmm
    implicit none



    !
    ! Global vector of registered rbf_mm_drivers
    !
    type(rbf_mm_driver_vector_t)          :: registered_rbf_mm_drivers
    logical                     :: initialized = .false.

contains


    !>  Register rbf_mm_drivers in a module vector.
    !!
    !!  This allows the available rbf_mm_drivers to be queried in the same way that they 
    !!  are registered for allocation. Adapted from mod_functions
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_rbf_mm_drivers()
        integer :: nrbf_mm_drivers, irbf_mm_driver

        !
        ! Instantiate rbf_mm_drivers
        !
        type(rbf_mm_driver_pmm)                               :: pmm 
        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_rbf_mm_drivers%push_back(pmm)
       
            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            nrbf_mm_drivers = registered_rbf_mm_drivers%size()
            do irbf_mm_driver = 1,nrbf_mm_drivers
                call registered_rbf_mm_drivers%data(irbf_mm_driver)%driver%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_rbf_mm_drivers
    !********************************************************************************************











    !> Factory method for allocating concrete rbf_mm_drivers
    !!
    !!      - Allocate a concrete rbf_mm_driver_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   rbf_mm_driver      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_rbf_mm_driver(rbf_mm_driver,string)
        class(rbf_mm_driver_t),  allocatable,    intent(inout)   :: rbf_mm_driver
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, rbf_mm_driverindex


        if ( allocated(rbf_mm_driver) ) then
            deallocate(rbf_mm_driver)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        rbf_mm_driverindex = registered_rbf_mm_drivers%index_by_name(trim(string))



        !
        ! Check rbf_mm_driver was found in 'registered_rbf_mm_drivers'
        !
        if (rbf_mm_driverindex == 0) call chidg_signal_one(FATAL,"create_rbf_mm_driver: rbf_mm_driver not recognized", trim(string))



        !
        ! Allocate conrete rbf_mm_driver_t instance
        !
        allocate(rbf_mm_driver, source=registered_rbf_mm_drivers%data(rbf_mm_driverindex)%driver, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_rbf_mm_driver: error allocating rbf_mm_driver from global vector.")



        !
        ! Check boundary condition was allocated
        !
        if ( .not. allocated(rbf_mm_driver) ) call chidg_signal(FATAL,"create_rbf_mm_driver: error allocating concrete rbf_mm_driver.")



    end subroutine create_rbf_mm_driver
    !*****************************************************************************************













    !>  Print a list of the registered rbf_mm_drivers. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'rbf_mm_driver_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_rbf_mm_drivers()
        integer                         :: nrbf_mm_drivers, irbf_mm_driver
        character(len=:),   allocatable :: rbf_mm_driver_name

        nrbf_mm_drivers = registered_rbf_mm_drivers%size()


        do irbf_mm_driver = 1,nrbf_mm_drivers

            rbf_mm_driver_name = registered_rbf_mm_drivers%data(irbf_mm_driver)%driver%get_name()
            call write_line(trim(rbf_mm_driver_name))

        end do ! irbf_mm_driver


    end subroutine list_rbf_mm_drivers
    !******************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function check_rbf_mm_driver_registered(state_string) result(state_found)
        character(len=*),   intent(in)  :: state_string

        integer(ik) :: state_index
        logical     :: state_found

        ! Find boundary condition string in 'registered_bcs' vector
        state_index = registered_rbf_mm_drivers%index_by_name(trim(state_string))

        ! Set status of state_found
        state_found = (state_index /= 0)

    end function check_rbf_mm_driver_registered
    !*******************************************************************************************************






end module mod_rbf_mm_driver
