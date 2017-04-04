module type_chidg_manager
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_string,         only: string_t
    use mod_chidg_mpi,      only: IRANK, NRANK, GLOBAL_MASTER, ChiDG_COMM
    use type_chidg,         only: chidg_t
    use type_svector,       only: svector_t
    use mod_wall_distance,  only: wall_distance_driver
    use mpi_f08,            only: MPI_Reduce, MPI_LOGICAL, MPI_LOR
    implicit none




    !>  Handle chidg_t instances for pre/post processing.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !!
    !-------------------------------------------------------------------------
    type, public :: chidg_manager_t


    contains

        procedure   :: process

    end type chidg_manager_t
    !*************************************************************************







contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/23/2016
    !!
    !!
    !-------------------------------------------------------------------------
    subroutine process(self,chidg)
        class(chidg_manager_t), intent(in)      :: self
        type(chidg_t),          intent(inout)   :: chidg

        character(:),   allocatable :: user_msg
        integer(ik)                 :: ifield, ierr, iproc
        type(svector_t)             :: auxiliary_fields_local
        type(string_t)              :: field_name

        logical     :: has_wall_distance, all_have_wall_distance
        
        !
        ! Assemble auxiliary fields across processors
        !
        auxiliary_fields_local = chidg%data%get_auxiliary_field_names()
        ! Send
        has_wall_distance = .false.
        do iproc = 0,NRANK-1
            do ifield = 1,auxiliary_fields_local%size()
                field_name = auxiliary_fields_local%at(ifield)
                has_wall_distance = (field_name%get() == 'Wall Distance : p-Poisson')
                if (has_wall_distance) exit
            end do

        end do !iproc


        call MPI_Reduce(has_wall_distance, all_have_wall_distance, 1, MPI_LOGICAL, MPI_LOR, GLOBAL_MASTER, ChiDG_COMM, ierr)


        if (all_have_wall_distance) then
            call wall_distance_driver(chidg,'wall_distance.h5')
        end if




!        !
!        ! Initialize auxiliary fields for ChiDG
!        !
!        auxiliary_fields_local = chidg%data%get_auxiliary_field_names()
!        do ifield = 1,auxiliary_fields_local%size()
!
!            field_name = auxiliary_fields_local%at(ifield)
!
!
!            select case ( field_name%get() )
!
!                case ('Wall Distance : p-Poisson')
!                    call wall_distance_driver(chidg,'wall_distance.h5')
!
!
!                case default
!                    user_msg = "chidg_manager%process: The current ChiDG simulation environment &
!                                seems to require an auxiliary field that we don't know how to provide. &
!                                These fields are set in the operator_t's that use them. So, make sure &
!                                the field is spelled correctly and also make sure that there is a rule &
!                                in the chidg_manager_t for obtaining/computing the field."
!                    call chidg_signal_one(FATAL,user_msg,trim(field_name%get()) )
!
!            end select
!
!
!
!        end do





    end subroutine process
    !*************************************************************************




end module type_chidg_manager
