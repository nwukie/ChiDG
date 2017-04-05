module type_mesh_new
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_chidg_mpi,              only: IRANK, NRANK, GLOBAL_MASTER
    use mpi_f08

    use type_point,                 only: point_t
    use type_domain,                only: domain_t
    use type_ivector,               only: ivector_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_element_connectivity,  only: element_connectivity_t
    implicit none
    private




    !> Data type for mesh information
    !!      - contains array of elements, array of faces for each element
    !!      - calls initialization procedure for elements and faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/27/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: mesh_new_t

        type(domain_t), allocatable :: domain(:)

    contains

        procedure   :: add_domain
        procedure   :: new_domain
        procedure   :: ndomains
        procedure   :: release

    end type mesh_new_t
    !****************************************************************************************





contains




    !>  Add a domain to the mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine add_domain(self, name, nodes, connectivity, nelements_g, spacedim, nterms_c, coord_system, eqn_ID)
        class(mesh_new_t),              intent(inout)   :: self
        character(*),                   intent(in)      :: name
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: nelements_g
        integer(ik),                    intent(in)      :: spacedim
        integer(ik),                    intent(in)      :: nterms_c
        character(*),                   intent(in)      :: coord_system
        integer(ik),                    intent(in)      :: eqn_ID

        integer(ik) :: idomain_l


        !
        ! Create new domain
        !
        idomain_l = self%new_domain()


        !
        ! Set name
        !
        self%domain(idomain_l)%name = trim(name)


        !
        ! Initialize new domain
        !
        call self%domain(idomain_l)%init_geom(idomain_l,    &
                                              nelements_g,  &
                                              spacedim,     &
                                              nterms_c,     &
                                              nodes,        &
                                              connectivity, &
                                              coord_system )

        
        !
        ! Set equation set identifier
        !
        call self%domain(idomain_l)%init_eqn(eqn_ID)

    end subroutine add_domain
    !***************************************************************************************














    !>  Extend domain storage by one and return the domain identifier of the newly 
    !!  allocated domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    function new_domain(self) result(idomain_l)
        class(mesh_new_t),   intent(inout)   :: self

        integer(ik)                 :: idomain_l, ierr
        type(domain_t), allocatable :: temp_domains(:)



        !
        ! Resize array storage
        !
        allocate(temp_domains(self%ndomains() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%ndomains() > 0) then
            temp_domains(1:size(self%domain)) = self%domain(1:size(self%domain))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_domains,self%domain)
        


        !
        ! Set domain identifier of newly allocated domain that will be returned
        !
        idomain_l = self%ndomains()


    end function new_domain
    !****************************************************************************************










    !> Return the number of domains in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function ndomains(self) result(ndomains_)
        class(mesh_new_t),    intent(in)      :: self

        integer :: ndomains_

        if (allocated(self%domain)) then
            ndomains_ = size(self%domain)
        else
            ndomains_ = 0
        end if

    end function ndomains
    !****************************************************************************************







    !>  Release any allocated resource on the object.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine release(self)
        class(mesh_new_t),  intent(inout)   :: self


        if (allocated(self%domain)) deallocate(self%domain)


    end subroutine release
    !***************************************************************************************










end module type_mesh_new
