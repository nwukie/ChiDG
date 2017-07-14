module mod_function
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_function,  only: function_t
    use type_fcnvector, only: fcnvector_t

    !
    ! Import functions
    !
    use fcn_xsquared,                       only: xsquared_f
    use fcn_ysquared,                       only: ysquared_f
    use fcn_zsquared,                       only: zsquared_f
    use fcn_xyz,                            only: xyz_f
    use fcn_radius,                         only: radius_f
    use fcn_gaussian,                       only: gaussian_f
    use fcn_constant,                       only: constant_f
    use fcn_sine,                           only: sine_f
    use fcn_polynomial,                     only: polynomial_f
    use fcn_scalar_adv_diff_bl_solution,    only: scalar_adv_diff_bl_solution_f
!    use fcn_cylindricalduct_eigenfunction,  only: cylindricalduct_eigenfunction_f
    use fcn_monopole,                       only: monopole_f
    implicit none



    !
    ! Global vector of registered functions
    !
    type(fcnvector_t)   :: registered_fcns
    logical             :: initialized = .false.

contains


    !>  Register functions in a module vector.
    !!
    !!  This allows the available functions to be queried in the same way that they 
    !!  are registered for allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_functions()
        integer :: nfcns, ifcn

        !
        ! Instantiate functions
        !
        type(xsquared_f)                        :: xsquared
        type(ysquared_f)                        :: ysquared
        type(zsquared_f)                        :: zsquared
        type(xyz_f)                             :: xyz
        type(radius_f)                          :: radius
        type(gaussian_f)                        :: gaussian
        type(constant_f)                        :: constant
        type(sine_f)                            :: sin_function
        type(polynomial_f)                      :: polynomial
        type(scalar_adv_diff_bl_solution_f)     :: scalar_adv_diff_bl_solution

        !
        ! Acoustics
        !
!        type(cylindricalduct_eigenfunction_f)   :: cylindricalduct_eigenfunction
        type(monopole_f)                        :: monopole

        

        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_fcns%push_back(xsquared)
            call registered_fcns%push_back(ysquared)
            call registered_fcns%push_back(zsquared)
            call registered_fcns%push_back(xyz)
            call registered_fcns%push_back(radius)
            call registered_fcns%push_back(gaussian)
            call registered_fcns%push_back(constant)
            call registered_fcns%push_back(sin_function)
            call registered_fcns%push_back(polynomial)
!            call registered_fcns%push_back(cylindricalduct_eigenfunction)
            call registered_fcns%push_back(monopole)
            call registered_fcns%push_back(scalar_adv_diff_bl_solution)


            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            nfcns = registered_fcns%size()
            do ifcn = 1,nfcns
                call registered_fcns%data(ifcn)%fcn%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_functions
    !********************************************************************************************











    !> Factory method for allocating concrete functions
    !!
    !!      - Allocate a concrete function_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   fcn      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_function(fcn,string)
        class(function_t),  allocatable,    intent(inout)   :: fcn
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, fcnindex


        if ( allocated(fcn) ) then
            deallocate(fcn)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        fcnindex = registered_fcns%index_by_name(trim(string))



        !
        ! Check function was found in 'registered_fcns'
        !
        if (fcnindex == 0) call chidg_signal_one(FATAL,"create_function: function not recognized", trim(string))



        !
        ! Allocate conrete function_t instance
        !
        allocate(fcn, source=registered_fcns%data(fcnindex)%fcn, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_function: error allocating function from global vector.")



        !
        ! Check boundary condition was allocated
        !
        if ( .not. allocated(fcn) ) call chidg_signal(FATAL,"create_function: error allocating concrete function.")



    end subroutine create_function
    !*****************************************************************************************













    !>  Print a list of the registered functions. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'function_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_functions()
        integer                         :: nfcns, ifcn
        character(len=:),   allocatable :: fcn_name

        nfcns = registered_fcns%size()


        do ifcn = 1,nfcns

            fcn_name = registered_fcns%data(ifcn)%fcn%get_name()
            call write_line(trim(fcn_name))

        end do ! ifcn


    end subroutine list_functions
    !******************************************************************************************









end module mod_function
