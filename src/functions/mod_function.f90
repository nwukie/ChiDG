module mod_function
    use mod_kinds,      only: rk, ik
    use atype_function, only: function_t

    ! IMPORT FUNCTIONS
    use fcn_xsquared,           only: xsquared_f
    use fcn_ysquared,           only: ysquared_f
    use fcn_zsquared,           only: zsquared_f
    use fcn_xyz,                only: xyz_f
    use fcn_gaussian,           only: gaussian_f
    use fcn_constant,           only: constant_f
    use fcn_isentropic_vortex,  only: isentropic_vortex_f
    implicit none

    ! LIST OF FUNCTION OBJECTS
    type(xsquared_f)            :: xsquared
    type(ysquared_f)            :: ysquared
    type(zsquared_f)            :: zsquared
    type(xyz_f)                 :: xyz
    type(gaussian_f)            :: gaussian
    type(constant_f)            :: constant
    type(isentropic_vortex_f)   :: isentropic_vortex


contains


    !> Factory method for allocating concrete functions
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   fcn     Incoming function base class to be allocated
    !!  @param[in]      str     Function string for selecting concrete type
    !------------------------------------------------------------------
    subroutine create_function(fcn,str)
        class(function_t), allocatable, intent(inout)   :: fcn
        character(*),                   intent(in)      :: str


        select case (trim(str))
            case ('xsquared')
                allocate(fcn,source=xsquared)

            case ('ysquared')
                allocate(fcn,source=ysquared)

            case ('zsquared')
                allocate(fcn,source=zsquared)

            case ('xyz')
                allocate(fcn,source=xyz)

            case ('gaussian')
                allocate(fcn,source=gaussian)

            case ('constant')
                allocate(fcn,source=constant)

            case ('isentropic_vortex','isentropic vortex')
                allocate(fcn,source=isentropic_vortex)

            case default
                stop "Error: assign_function -- function string not recognized"

        end select




    end subroutine




end module mod_function
