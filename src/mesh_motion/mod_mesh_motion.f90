module mod_mesh_motion
#include <messenger.h>
    use mod_kinds,                              only: rk, ik
    use type_mesh_motion,   only: mesh_motion_t
    use type_mmvector,                        only: mmvector_t

    !
    ! Import mesh_motions
    !
    use type_prescribed_mesh_motion,            only: prescribed_mesh_motion_t
    use type_rbf_mesh_motion,                   only: rbf_mesh_motion_t
    implicit none



    !
    ! Global vector of registered mesh_motions
    !
    type(mmvector_t)          :: registered_mms
    logical                     :: initialized = .false.

contains


    !>  Register mesh_motions in a module vector.
    !!
    !!  This allows the available mesh_motions to be queried in the same way that they 
    !!  are registered for allocation. Adapted from mod_functions
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_mesh_motions()
        integer :: nmms, imm

        !
        ! Instantiate mesh_motions
        !
        type(prescribed_mesh_motion_t)                                  :: prescribed_mesh_motion
        type(rbf_mesh_motion_t)                                         :: rbf_mesh_motion

        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_mms%push_back(prescribed_mesh_motion)
            call registered_mms%push_back(rbf_mesh_motion)
      
            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            nmms = registered_mms%size()
            do imm = 1,nmms
                call registered_mms%data(imm)%mm%init_name()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_mesh_motions
    !********************************************************************************************











    !> Factory method for allocating concrete mesh_motions
    !!
    !!      - Allocate a concrete mesh_motion_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   mm      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_mesh_motion(mm,string)
        class(mesh_motion_t),  allocatable,    intent(inout)   :: mm
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, mmindex


        if ( allocated(mm) ) then
            deallocate(mm)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        mmindex = registered_mms%index_by_name(trim(string))



        !
        ! Check mesh_motion was found in 'registered_mms'
        !
        if (mmindex == 0) call chidg_signal_one(FATAL,"create_mesh_motion: mesh_motion not recognized", trim(string))



        !
        ! Allocate conrete mesh_motion_t instance
        !
        allocate(mm, source=registered_mms%data(mmindex)%mm, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_mesh_motion: error allocating mesh_motion from global vector.")



        !
        ! Check boundary condition was allocated
        !
        if ( .not. allocated(mm) ) call chidg_signal(FATAL,"create_mesh_motion: error allocating concrete mesh_motion.")



    end subroutine create_mesh_motion
    !*****************************************************************************************













    !>  Print a list of the registered mesh_motions. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'mesh_motion_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_mesh_motions()
        integer                         :: nmms, imm
        character(len=:),   allocatable :: mm_name

        nmms = registered_mms%size()


        do imm = 1,nmms

            mm_name = registered_mms%data(imm)%mm%get_family_name()
            call write_line(trim(mm_name))

        end do ! imm


    end subroutine list_mesh_motions
    !******************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function check_mm_registered(state_string) result(state_found)
        character(len=*),   intent(in)  :: state_string

        integer(ik) :: state_index
        logical     :: state_found

        ! Find boundary condition string in 'registered_bcs' vector
        state_index = registered_mms%index_by_name(trim(state_string))

        ! Set status of state_found
        state_found = (state_index /= 0)

    end function check_mm_registered
    !*******************************************************************************************************






end module mod_mesh_motion
