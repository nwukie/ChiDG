module mod_prescribed_mesh_motion
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_prescribed_mesh_motion,  only: prescribed_mesh_motion_t
    use type_pmmvector, only: pmmvector_t

    !
    ! Import prescribed_mesh_motions
    !
    use pmm_static,                       only: static_pmm
    implicit none



    !
    ! Global vector of registered prescribed_mesh_motions
    !
    type(pmmvector_t)   :: registered_pmms
    logical             :: initialized = .false.

contains


    !>  Register prescribed_mesh_motions in a module vector.
    !!
    !!  This allows the available prescribed_mesh_motions to be queried in the same way that they 
    !!  are registered for allocation. Adapted from mod_functions
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_prescribed_mesh_motions()
        integer :: npmms, ipmm

        !
        ! Instantiate prescribed_mesh_motions
        !
        type(static_pmm)                        :: static
        

        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_pmms%push_back(static)
       
            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            npmms = registered_pmms%size()
            do ipmm = 1,npmms
                call registered_pmms%data(ipmm)%pmm%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_prescribed_mesh_motions
    !********************************************************************************************











    !> Factory method for allocating concrete prescribed_mesh_motions
    !!
    !!      - Allocate a concrete prescribed_mesh_motion_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   pmm      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_prescribed_mesh_motion(pmm,string)
        class(prescribed_mesh_motion_t),  allocatable,    intent(inout)   :: pmm
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, pmmindex


        if ( allocated(pmm) ) then
            deallocate(pmm)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        pmmindex = registered_pmms%index_by_name(trim(string))



        !
        ! Check prescribed_mesh_motion was found in 'registered_pmms'
        !
        if (pmmindex == 0) call chidg_signal_one(FATAL,"create_prescribed_mesh_motion: prescribed_mesh_motion not recognized", trim(string))



        !
        ! Allocate conrete prescribed_mesh_motion_t instance
        !
        allocate(pmm, source=registered_pmms%data(pmmindex)%pmm, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_prescribed_mesh_motion: error allocating prescribed_mesh_motion from global vector.")



        !
        ! Check boundary condition was allocated
        !
        if ( .not. allocated(pmm) ) call chidg_signal(FATAL,"create_prescribed_mesh_motion: error allocating concrete prescribed_mesh_motion.")



    end subroutine create_prescribed_mesh_motion
    !*****************************************************************************************













    !>  Print a list of the registered prescribed_mesh_motions. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'prescribed_mesh_motion_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_prescribed_mesh_motions()
        integer                         :: npmms, ipmm
        character(len=:),   allocatable :: pmm_name

        npmms = registered_pmms%size()


        do ipmm = 1,npmms

            pmm_name = registered_pmms%data(ipmm)%pmm%get_name()
            call write_line(trim(pmm_name))

        end do ! ipmm


    end subroutine list_prescribed_mesh_motions
    !******************************************************************************************









end module mod_prescribed_mesh_motion
