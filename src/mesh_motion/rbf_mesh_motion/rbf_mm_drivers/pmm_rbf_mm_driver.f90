! This extended derived type implements procedures for updating RBF mesh motion
! source node values for displacement and velocity based on a prescribed mesh motion function.
!
!
module pmm_rbf_mm_driver
#include <messenger.h>
    use mod_kinds,                  only: ik, rk
    use type_rbf_mm_driver,         only: rbf_mm_driver_t
    use mod_prescribed_mesh_motion_function,    only: create_prescribed_mesh_motion_function
    use type_prescribed_mesh_motion_function,   only: prescribed_mesh_motion_function_t
    implicit none

    type, extends(rbf_mm_driver_t) :: rbf_mm_driver_pmm

        character(:), allocatable                                   :: pmmf_name
        class(prescribed_mesh_motion_function_t), allocatable       :: pmmf

    contains

        ! Deferred Procedures
        procedure :: init
        procedure :: compute_disp
        procedure :: compute_vel
 
        ! Specialized Procedures
        procedure       :: set_pmmf_name
        procedure       :: add_pmmf


    end type rbf_mm_driver_pmm

contains


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine init(self)
        class(rbf_mm_driver_pmm),       intent(inout) :: self
        
        call self%set_name('pmm')

    end subroutine init
    !**************************************************************



    !>
    !!
    !!
    !!
    !--------------------------------------------------------------
    function compute_disp(self,time,node) result(val)
        class(rbf_mm_driver_pmm),   intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)
        real(rk)                                    :: val(3)

        val = self%pmmf%compute_pos(time,node) - node

    end function compute_disp 
    !**************************************************************


    !>
    !!
    !!
    !!
    !--------------------------------------------------------------
    function compute_vel(self,time,node) result(val)
        class(rbf_mm_driver_pmm),   intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: node(3)
        real(rk)                                    :: val(3)

        val = self%pmmf%compute_vel(time,node) 
             
    end function compute_vel
    !**************************************************************

    
    !>
    !!
    !!
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine set_pmmf_name(self, pmmfstring)
        class(rbf_mm_driver_pmm),   intent(inout)   :: self
        character(*),               intent(in)      :: pmmfstring

        self%pmmf_name = pmmfstring

    end subroutine
    !********************************************************************************


    !>
    !!  Called from get_pmm_hdf
    !! 
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine add_pmmf(self, pmmf_in)
        class(rbf_mm_driver_pmm),        intent(inout)   :: self
        class(prescribed_mesh_motion_function_t), intent(inout)               :: pmmf_in
            
        integer(ik)     :: ierr
        if (allocated(self%pmmf)) deallocate(self%pmmf)
        allocate(self%pmmf, source = pmmf_in, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine
    !********************************************************************************



end module pmm_rbf_mm_driver
