module mod_radial_basis_function
#include <messenger.h>
    use mod_kinds,                              only: rk, ik
    use type_radial_basis_function,   only: radial_basis_function_t
    use type_rbfvector,                        only: rbfvector_t

    !
    ! Import radial_basis_functions
    !
    use rbf_tps,                            only: tps_rbf
    use rbf_wc0,                            only: wc0_rbf
    use rbf_wc2,                            only: wc2_rbf
    use rbf_wc4,                            only: wc4_rbf
    use rbf_wc6,                            only: wc6_rbf
    implicit none



    !
    ! Global vector of registered radial_basis_functions
    !
    type(rbfvector_t)          :: registered_rbfs
    logical                     :: initialized = .false.

contains


    !>  Register radial_basis_functions in a module vector.
    !!
    !!  This allows the available radial_basis_functions to be queried in the same way that they 
    !!  are registered for allocation. Adapted from mod_functions
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_radial_basis_functions()
        integer :: nrbfs, irbf

        !
        ! Instantiate radial_basis_functions
        !
        type(tps_rbf)                               :: tps 
        type(wc0_rbf)                               :: wc0 
        type(wc2_rbf)                               :: wc2 
        type(wc4_rbf)                               :: wc4 
        type(wc6_rbf)                               :: wc6 

        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_rbfs%push_back(tps)
            call registered_rbfs%push_back(wc0)
            call registered_rbfs%push_back(wc2)
            call registered_rbfs%push_back(wc4)
            call registered_rbfs%push_back(wc6)
      
            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            nrbfs = registered_rbfs%size()
            do irbf = 1,nrbfs
                call registered_rbfs%data(irbf)%rbf%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_radial_basis_functions
    !********************************************************************************************











    !> Factory method for allocating concrete radial_basis_functions
    !!
    !!      - Allocate a concrete radial_basis_function_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   rbf      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_radial_basis_function(rbf,string)
        class(radial_basis_function_t),  allocatable,    intent(inout)   :: rbf
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, rbfindex


        if ( allocated(rbf) ) then
            deallocate(rbf)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        rbfindex = registered_rbfs%index_by_name(trim(string))



        !
        ! Check radial_basis_function was found in 'registered_rbfs'
        !
        if (rbfindex == 0) call chidg_signal_one(FATAL,"create_radial_basis_function: radial_basis_function not recognized", trim(string))



        !
        ! Allocate conrete radial_basis_function_t instance
        !
        allocate(rbf, source=registered_rbfs%data(rbfindex)%rbf, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_radial_basis_function: error allocating radial_basis_function from global vector.")



        !
        ! Check boundary condition was allocated
        !
        !if ( .not. allocated(rbf) ) call chidg_signal(FATAL,"create_radial_basis_function: error allocating concrete radial_basis_function.")



    end subroutine create_radial_basis_function
    !*****************************************************************************************













    !>  Print a list of the registered radial_basis_functions. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'radial_basis_function_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_radial_basis_functions()
        integer                         :: nrbfs, irbf
        character(len=:),   allocatable :: rbf_name

        nrbfs = registered_rbfs%size()


        do irbf = 1,nrbfs

            rbf_name = registered_rbfs%data(irbf)%rbf%get_name()
            call write_line(trim(rbf_name))

        end do ! irbf


    end subroutine list_radial_basis_functions
    !******************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function check_rbf_registered(state_string) result(state_found)
        character(len=*),   intent(in)  :: state_string

        integer(ik) :: state_index
        logical     :: state_found

        ! Find boundary condition string in 'registered_bcs' vector
        state_index = registered_rbfs%index_by_name(trim(state_string))

        ! Set status of state_found
        state_found = (state_index /= 0)

    end function check_rbf_registered
    !*******************************************************************************************************





end module mod_radial_basis_function
