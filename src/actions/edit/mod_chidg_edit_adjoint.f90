!>
!!
!!  @author Matteo Ugolotti
!!  @data   05/08/2017
!!
!!
!!
!---------------------------------------------------------------------------------
module mod_chidg_edit_adjoint
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use hdf5
    use h5lt
    use mod_hdf_utilities,              only: check_link_exists_hdf,create_functional_hdf,               &
                                              create_functional_group_hdf, open_functional_group_hdf,    &
                                              close_functional_group_hdf, get_adjoint_status_hdf,        &
                                              remove_functional_hdf, remove_functional_group_hdf,        &
                                              set_functional_auxiliary_geom_hdf, set_functional_LS_hdf,  &
                                              set_functional_reference_geom_hdf, get_nfunctionals_hdf,   &  
                                              check_functional_exists_hdf, get_functionals_names_hdf,    &
                                              get_functional_auxiliary_geom_hdf, set_adjoint_status_hdf, &
                                              get_functional_reference_geom_hdf, get_functional_LS_hdf

    use mod_functional,                 only: list_functionals, check_functional_existence
    use mod_chidg_edit_printoverview,   only: print_overview, chidg_clear_screen, print_line_separator
    use type_svector,   only: svector_t
    use mod_string,     only: string_t
    implicit none

    !-----------------------------------------------------------------------------------------
    !!
    !!  chidg_edit_adjoint
    !!  chidg_edit_adjoint_selector
    !!  chidg_toggle_adjoint
    !!  chidg_edit_adjoint_add
    !!  chidg_edit_adjoint_edit
    !!  chidg_edit_adjoint_remove
    !!  print_functional_overview
    !!
    !-----------------------------------------------------------------------------------------



contains

    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/08/2017
    !!
    !-----------------------------------------------------------------------------
    subroutine chidg_edit_adjoint(fid)
    integer(HID_T), intent(in)  :: fid

    logical         :: run_fcl_edit
    integer(ik)     :: int_input, ierr

    run_fcl_edit = .true.

    do while (run_fcl_edit)
        
        call chidg_clear_screen()
        call print_overview(fid)
        call print_functional_overview(fid)

        call write_line(' ')
        call write_line('Select command:')
        call write_line("1: toggle adjoint","2: add functional","3: edit functional","4: remove functional", "0: exit", columns=.true., column_width=25, color='blue' )


        read(*,'(I8)', iostat=ierr) int_input
        if ( (ierr/=0) .or. (abs(int_input)>3) ) print*,"Invalid input: expecting 0, 1, 2, 3 or 4."
        
        ! Select case and run specific procedure
        select case (int_input)
            case (0)
                run_fcl_edit = .false.
            case (1)
                call chidg_edit_adjoint_selector(fid,"adjoint")
            case (2)
                call chidg_edit_adjoint_selector(fid,"add")
            case (3)
                call chidg_edit_adjoint_selector(fid,"edit")
            case (4)
                call chidg_edit_adjoint_selector(fid,"remove")
            case default

        end select

    end do

    end subroutine chidg_edit_adjoint
    !*****************************************************************************










    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/08/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine chidg_edit_adjoint_selector(fid,operation)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: operation

        character(1024)     :: functional_str
        logical             :: run_add,run_edit,run_remove,list_functional, func_exists

        ! Print overview
        call chidg_clear_screen()
        call print_overview(fid)
        call print_functional_overview(fid)


        ! Add / Edit / Remove cases
        select case (operation)
            case ("adjoint")
                call chidg_toggle_adjoint(fid)

            case ("add")
                
                list_functional = .false.
                run_add = .true.
                do while (run_add)
                    
                    ! In case of help, clean and print the screen again with the list 
                    if (list_functional) then
                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_functional_overview(fid)
                        call list_functionals()
                    end if


                    call write_line(' ')
                    call write_line('Enter new functional name (? to list):',color='blue')
                    read(*,"(A1024)") functional_str
                    
                    list_functional = (trim(functional_str) == "?")     ! Check for help
                    run_add = (trim(functional_str)/="")                ! Check for exit
                    
                    if ( run_add .and. (.not. list_functional)) then
                        call chidg_edit_functional_add(fid,trim(functional_str),run_add,list_functional)
                    end if

                end do


            case ("edit")
                
                run_edit = .true.
                do while (run_edit)

                    call write_line(' ')
                    call write_line('Enter name of the functional to edit:',color='blue')
                    
                    read(*,"(A1024)") functional_str
                    
                    run_edit = (trim(functional_str)/="")   ! Check for exit
                    ! Check that the user typed an active functional
                    func_exists = check_functional_exists_hdf(fid,functional_str)
                                            
                    if ( run_edit .and. func_exists) then
                        call chidg_edit_functional_edit(fid,trim(functional_str))
                        run_edit=.false.
                    end if

                end do


            case ("remove")
                
                run_remove = .true.
                do while (run_remove)

                    call write_line(' ')
                    call write_line('Enter the name of the function to remove:',color='blue')
                    
                    read(*,"(A1024)") functional_str
                    
                    ! Check for user exit
                    run_remove = (trim(functional_str)/="")
                    ! Check that the user typed an active functional
                    func_exists = check_functional_exists_hdf(fid,functional_str)
                                            
                    if ( run_remove .and. func_exists) then
                        call chidg_edit_functional_remove(fid,trim(functional_str))
                        run_remove=.false.
                    end if

                end do

            case default

        end select

    end subroutine chidg_edit_adjoint_selector
    !****************************************************************************************** 
    






    !> Subroutine to toggle adjoint mode (ON/OFF)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   07/07/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine chidg_toggle_adjoint(fid)
        integer(HID_T), intent(in)  :: fid

        character(1024)         :: enter
        character, allocatable  :: adjoint_mode
        logical                 :: adjoint_status
        integer(ik)             :: nfunc

        ! Print overview
        call chidg_clear_screen()
        call print_overview(fid)
        call print_functional_overview(fid)
    
        ! Get current status
        adjoint_status = get_adjoint_status_hdf(fid)

        ! check if any functional has been added
        nfunc = get_nfunctionals_hdf(fid)

        if (nfunc == 0) then
            call write_line('Warning message:',color='blue')
            call write_line('Toggle ignored: No functionals have been added. Add functionals before toggling adjoint. Press enter to continue.',color='red')
            read(*,"(A1024)") enter
        else
            ! Toggle adjoint status
            if (adjoint_status) then
                call set_adjoint_status_hdf(fid,"OFF")
            else
                call set_adjoint_status_hdf(fid,"ON")
            end if
        end if
    
    end subroutine chidg_toggle_adjoint
    !****************************************************************************************








    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/08/2017
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine chidg_edit_functional_add(fid,fcl_name,enter_again,list_fcls)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: fcl_name
        logical,            intent(out) :: enter_again
        logical,            intent(out) :: list_fcls

        logical             :: fcl_group_exists,fcl_exists, fcl_registered
        integer(HID_T)      :: fcl_id
        integer(ik)         :: ierr
        character(len=1024) :: enter

        ! Check if the functional entered is registered
        fcl_registered = check_functional_existence(fcl_name)


        if (fcl_registered) then

            ! Create the functional group if it does not exist.
            call create_functional_group_hdf(fid)
            
            ! Check if the functional exists, if not create it
            call create_functional_hdf(fid,trim(fcl_name))
            
            ! Go back to the functional menu
            enter_again = .false.
            list_fcls   = .false.

        else

            call write_line("WARNING: The functional entered is not registered!",color='blue')
            call write_line("Press enter to proceed and select a registered functional in the list displayed.",color='blue')

            read(*,'(A1024)') enter

            ! Go back to add menu and list the functional available
            enter_again = .true.
            list_fcls   = .true.

        end if

    end subroutine chidg_edit_functional_add
    !*******************************************************************************************




    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/08/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine chidg_edit_functional_remove(fid,fcl_name)
        integer(HID_T),     intent(in):: fid
        character(*),       intent(in):: fcl_name

        logical         :: functional_on,fcl_exists
        integer(HID_T)  :: fcl_id
        integer(ik)     :: ierr
        
        ! Check if there are any functionals
        functional_on = (get_nfunctionals_hdf(fid) /= 0)

        if (functional_on) then
            
            ! Open functional group
            fcl_id = open_functional_group_hdf(fid)
            
            ! Check if the functional exists, if not create it
            call remove_functional_hdf(fcl_id,trim(fcl_name))
            
            ! Close functional group
            call close_functional_group_hdf(fcl_id)

            ! Check if the number of functionals became 0, if so delete the functional group
            ! If the number onf function is 0, functional mode is automatically turned OFF.
            functional_on = ( get_nfunctionals_hdf(fid) == 0)
            if (functional_on) call remove_functional_group_hdf(fid)

        end if
        
    end subroutine chidg_edit_functional_remove
    !*******************************************************************************************



    
    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/08/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine chidg_edit_functional_edit(fid,fcl_name)
        integer(HID_T),     intent(in)  :: fid
        character(*),       intent(in)  :: fcl_name
        
        logical             :: run_edit
        character(len=1024) :: edit_option
        integer(ik)         :: input_value
        character(len=1024) :: input_string
        
        run_edit = .true.
        do while (run_edit)
       
            call chidg_clear_screen()
            call print_overview(fid)
            call print_functional_overview(fid,active_func=fcl_name)
       
            call write_line(' ')
            call write_line('1: Edit reference geometry','2: Edit auxiliary geometry','3: Edit linear solver',columns=.true.,column_width=30,color='blue')
            
            read(*,"(A1024)") edit_option
            
            run_edit = (trim(edit_option) /= "")
                  
            if (run_edit) then
                select case (edit_option)
                    case ("1")
                        
                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_functional_overview(fid,active_func=fcl_name)
               
                        call write_line(' ')
                        call write_line('Enter the new reference geometry:',color='blue')
                        read(*,"(A1024)") input_string

                        if (trim(input_string) /= '') call set_functional_reference_geom_hdf(fid,fcl_name,trim(input_string))
                        

                    case ("2")

                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_functional_overview(fid,active_func=fcl_name)
               
                        call write_line(' ')
                        call write_line('Enter the new auxiliary geometry:',color='blue')
                        read(*,"(A1024)") input_string

                        if (trim(input_string) /= '') call set_functional_auxiliary_geom_hdf(fid,fcl_name,trim(input_string))

                    case ("3")

                        call chidg_clear_screen()
                        call print_overview(fid)
                        call print_functional_overview(fid,active_func=fcl_name)
               
                        call write_line(' ')
                        call write_line('Enter the new linear solver:',color='blue')
                        read(*,"(A1024)") input_string

                        if (trim(input_string) /= '') call set_functional_LS_hdf(fid,fcl_name,trim(input_string))

                    case default


                end select

            end if

        end do


    end subroutine chidg_edit_functional_edit
    !*******************************************************************************************







    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/8/2017
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine print_functional_overview(fid,active_func)
        integer(HID_T), intent(in)            :: fid
        character(*),   intent(in), optional  :: active_func

        logical                             :: adj_status_on
        integer(ik)                         :: num_functionals,ifunc
        character(len=1024)                 :: func_name,ref_geom,aux_geom,func_LS
        type(svector_t)                     :: func_names
        
        call write_line(' ')
        call write_line(' ')
        
        ! Check "Adjoint Mode" attribute status, if the attribute does not exists it creates it
        adj_status_on = get_adjoint_status_hdf(fid)
        call add_to_line(":Functionals  --", ltrim=.false.)
        if (adj_status_on) then
            call add_to_line("Adjoint(on)", ltrim=.false., color='green')
        else
            call add_to_line("Adjoint(off)", ltrim=.false., color='red')
        end if
        call send_line()

        call print_line_separator()
        call write_line('Functional', 'Reference geometry', 'Auxiliary geometry', 'Linear lolver', delimiter='  :  ', columns=.True., column_width=20)
        call print_line_separator()
        
        ! Check if any functional is registered
        num_functionals = get_nfunctionals_hdf(fid)

        if (num_functionals == 0) then
            call write_line('-','-','-','-', delimiter='  :  ',columns=.True., column_width=20)
        
        else
            
            func_names = get_functionals_names_hdf(fid)
            do ifunc = 1,num_functionals
                
                func_name   = func_names%data_(ifunc)%get()
                aux_geom    = get_functional_auxiliary_geom_hdf(fid,func_name)
                ref_geom    = get_functional_reference_geom_hdf(fid,func_name)
                func_LS     = get_functional_LS_hdf(fid,func_name)
                
                if ( present(active_func) .and. (active_func==func_name) ) then
                    call write_line(func_name,ref_geom,aux_geom,func_LS,delimiter='  :  ',columns=.True., column_width=20, color='pink')
                else
                    call write_line(func_name,ref_geom,aux_geom,func_LS,delimiter='  :  ',columns=.True., column_width=20)
                end if

            end do !ifunc
        end if
        
        call write_line(" ")

    end subroutine print_functional_overview
    !*******************************************************************************************************




end module mod_chidg_edit_adjoint
