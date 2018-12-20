module type_mesh_motion
#include <messenger.h>
    use mod_constants
    use mod_kinds,                              only: rk, ik
    use type_chidg_worker,                      only: chidg_worker_t
    use type_domain,                            only: domain_t
    use type_mesh,                              only: mesh_t
    implicit none


    !>  Mesh motion class. Abstract interface to enable treatment of prescribed mesh motion.
    !!
    !!  Data flow:
    !!  gridfile :-- readmeshmotion_hdf --> mm_group, mm_domain_data ...
    !!  :-- init_mm_group, init_mm_domain --> mesh_motion(:), mesh
    !!
    !!  @author Eric Wolf
    !!  @date   3/16/2017
    !!
    !---------------------------------------------------------------------------------------------------------------


    type, public, abstract     :: mesh_motion_t
        
        ! mesh motions are stored as a vector of mesh_motion_t's
        ! chidg%data%mesh%mesh_motion(1:nmm)

        !mm_ID is the unique ID for this pmm, so that it can be accessed via
        ! chidg%data%mesh%mesh_motion(mm_ID)

        integer(ik)                                                 :: mm_ID

        !mm_name is a name assigned to this mm
        ! This is used to identify the polymorphic type of the MM. 
        character(:), allocatable                                   :: mm_name
        character(:), allocatable                                   :: family_name


        real(rk)                                                    :: time = 0.0_rk

    contains

        ! Deferred procedures
        procedure(init_mm_name_interface),   deferred    :: init_name
        procedure(init_mm_interface),   deferred    :: init 
        procedure(update_mm_interface), deferred    :: update
        procedure(apply_mm_interface),  deferred    :: apply 

        ! Standard set/get procedures for pmm_name
        procedure   :: set_name
        procedure   :: get_name

        procedure   :: set_family_name
        procedure   :: get_family_name
        !These procedures read in the intermediate mm  data structures
        !created when reading in the grid file, which are
        ! mm_group, which contains a mm instance that is used as an
        ! allocation source for the present mm, and 
        ! mm_domain_data, which contains a domain name and a mm_ID, 
        ! which is used to store the mm_ID of the mesh in the appropriate domain
        procedure   :: init_mm_group 
        procedure   :: init_mm_domain


    end type mesh_motion_t

     abstract interface
         subroutine init_mm_name_interface(self)
            import mesh_motion_t

            class(mesh_motion_t),  intent(inout)   :: self
         end subroutine init_mm_name_interface



         subroutine init_mm_interface(self, mesh)
            use mod_kinds,  only: rk
            import mesh_motion_t
            import mesh_t

            class(mesh_motion_t),  intent(inout)   :: self
            type(mesh_t),         intent(inout)   :: mesh 
        end subroutine init_mm_interface


        subroutine update_mm_interface(self, mesh, time)
            use mod_kinds,  only: rk
            import mesh_motion_t
            import mesh_t

            class(mesh_motion_t),  intent(inout)   :: self
            type(mesh_t),         intent(inout)   :: mesh 
            real(rk),               intent(in)      :: time
        end subroutine update_mm_interface

        subroutine apply_mm_interface(self, mesh, time)
            use mod_kinds,  only: rk
            import mesh_motion_t
            import mesh_t

            class(mesh_motion_t),  intent(inout)   :: self
            type(mesh_t),         intent(inout)   :: mesh 
            real(rk),               intent(in)      :: time
        end subroutine apply_mm_interface 

    end interface

   
contains

    !>
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine set_name(self,mm_name)
        class(mesh_motion_t),            intent(inout)       :: self
        character(*),                               intent(in)          :: mm_name                        


        self%mm_name = mm_name


    end subroutine set_name
    !********************************************************************************

    !>
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    function get_name(self) result(mm_name)
        class(mesh_motion_t),            intent(inout)       :: self

        character(len=:), allocatable                                   :: mm_name                        


        mm_name =  self%mm_name

    end function get_name
    !********************************************************************************

    !>
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine set_family_name(self,mm_name)
        class(mesh_motion_t),            intent(inout)       :: self
        character(*),                               intent(in)          :: mm_name                        


        self%family_name = mm_name


    end subroutine set_family_name
    !********************************************************************************

    !>
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    function get_family_name(self) result(mm_name)
        class(mesh_motion_t),            intent(inout)       :: self

        character(len=:), allocatable                                   :: mm_name                        


        mm_name =  self%family_name

    end function get_family_name
    !********************************************************************************



    !>
    !!  This subroutine takes a an input mm (from a mm_group, which is generated by reading in the
    !!  grid file and initializing a mm instance according to a MM group),
    !!  and uses this mm instance as an allocation source for the present mm.
    !!
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine init_mm_group(self,mm_in)
        class(mesh_motion_t),               intent(inout)       :: self
        class(mesh_motion_t),               intent(inout)       :: mm_in


        integer(ik)     :: ierr
        

        call self%set_name(mm_in%get_name())

    end subroutine init_mm_group
    !********************************************************************************

    
    !>
    !!  This subroutine takes in a mesh (domain) that was selected from a
    !!  mm_domain_data_t, whch is generated by reading in the grid file,
    !!  and sets the mm_ID for all elements and faces in the mesh (domain)
    !!  to be equal to the present mm_ID.
    !!
    !!  Called from chidg%data%add_mm_domain(..)
    !!
    !!  @author Eric Wolf
    !!  @date 4/7/2017
    !--------------------------------------------------------------------------------
    subroutine init_mm_domain(self,domain)
        class(mesh_motion_t),               intent(inout)       :: self
        type(domain_t),                                       intent(inout)           ::domain 


        integer(ik)     :: ielem, nelem, iface
        
        ! Loop through the elements and faces in the mesh and assign the mm_ID


        domain%mm_ID = self%mm_ID
        nelem = domain%nelem
        do ielem = 1, nelem
            !Check if a MM has already been assigned
            !if (mesh%elems(ielem)%mm_ID /= NO_MM_ASSIGNED)  then
            !   bad bad bad
            !else
            domain%elems(ielem)%mm_ID = self%mm_ID
            do iface = 1, NFACES
                domain%faces(ielem, iface)%mm_ID = self%mm_ID
            end do
            !end if
        end do

    end subroutine init_mm_domain
    !********************************************************************************




end module type_mesh_motion
