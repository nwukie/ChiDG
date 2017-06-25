!>  A submodule of type_chidg providing the 'auxiliary_driver' subroutine
!!  to drive chidg auxiliary_environment instances.
!!
!!  If you would like to provide another auxiliary driver routine: 
!!  --------------------------------------------------------------
!!      1: create a few file that is a submodule of type_chidg_driver.
!!         type_chidg_driver is itself a submodule of type_chidg so the
!!         first line should look something like:
!!  
!!         submodule (type_chidg:type_chidg_driver) my_driver
!!
!!      2: Make sure your routine in the submodule is declared as a
!!         'module subroutine'
!!
!!      3: In this file, declare your driver interface under the 
!!         "AUXILIARY DRIVER SUBMODULE INTERFACES" text.
!!
!!      4: Register your new driver under the 'auxiliary_driver' routine
!!         in this file.
!!
!!      NOTE: you can use wall_distance_driver as a template
!!
!!
!!  @author Nathan A. Wukie
!!  @date   6/25/2017
!!
!--------------------------------------------------------------------------------
submodule (type_chidg) type_chidg_driver
    implicit none


    !-----------------------------------------------
    !    AUXILIARY DRIVER SUBMODULE INTERFACES
    !-----------------------------------------------
    
    ! Provided by driver_wall_distance.f90 submodule
    interface
        module subroutine wall_distance_driver(chidg,wall_distance,file_name)
            type(chidg_t),  intent(inout)   :: chidg
            type(chidg_t),  intent(inout)   :: wall_distance
            character(*),   intent(in)      :: file_name
        end subroutine wall_distance_driver
    end interface



contains


    !>  An interface for calling auxiliary driver routines.
    !!
    !!  Called in chidg%prerun()
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2017
    !!
    !----------------------------------------------------------------------------------
    module subroutine auxiliary_driver(chidg,chidg_aux,case,file_name)
        type(chidg_t),  intent(inout)   :: chidg
        type(chidg_t),  intent(inout)   :: chidg_aux
        character(*),   intent(in)      :: case
        character(*),   intent(in)      :: file_name

        select case(trim(case))
            case('Wall Distance')
                call wall_distance_driver(chidg        =chidg,      &
                                          wall_distance=chidg_aux,  &
                                          file_name    =file_name)

        end select

    end subroutine auxiliary_driver
    !**********************************************************************************













end submodule type_chidg_driver
