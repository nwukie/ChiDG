module type_chidg_data
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use type_point,                 only: point_t
    use type_dict,                  only: dict_t

    ! Primary chidg_data_t components
    use type_mesh,                  only: mesh_t
    use type_bcset,                 only: bcset_t
    use type_equationset_wrapper,   only: equationset_wrapper_t
    use type_solverdata,            only: solverdata_t

    ! Factory methods
    use mod_equations,              only: create_equationset
    use mod_bc,                     only: create_bc

    ! Classes
    use type_bc,                    only: bc_t

    implicit none


    !> solver abstract type definition
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public  :: chidg_data_t

        integer(ik)                         :: ndomains = 0
        type(dict_t)                        :: domain_info      !< Dictionary of (domain_index, domain_name) pairs


        type(mesh_t),                   allocatable :: mesh(:)
        type(bcset_t),                  allocatable :: bcset(:)
        type(equationset_wrapper_t),    allocatable :: eqnset(:)


        type(solverdata_t)                          :: sdata


        logical                         :: solverInitialized = .false.


    contains

        procedure   :: init_sdata

        procedure   :: add_domain
        procedure   :: add_bc

        

    end type chidg_data_t
    !---------------------------------------------------------------------------------------




contains


    !> Initialize solution data storage structures
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_sdata(self)
        class(chidg_data_t),     intent(inout)   :: self


        call self%sdata%init(self%mesh)


    end subroutine








    !> Add a domain to ChiDG. Calls initialization routines for components which define a 
    !! domain. That is, a mesh_t, an equationset_t, and the number of terms in their 
    !! polynomial expansions.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  points      point_t matrix defining the mesh
    !!  @param[in]  nterms_c    Integer defining the number of terms in the element coordinate expansion
    !!  @param[in]  eqnset      Character string defining the equationset_t for the domain
    !!  @param[in]  nterms_s    Integer defining the number of terms in the solution expansion
    !!
    !---------------------------------------------------------------------------------------
    subroutine add_domain(self,name,points,nterms_c,eqnset,nterms_s)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: name
        type(point_t),          intent(in)      :: points(:,:,:)
        integer(ik),            intent(in)      :: nterms_c
        character(*),           intent(in)      :: eqnset
        integer(ik),            intent(in)      :: nterms_s

        integer(ik) :: idom, ierr


        type(mesh_t),                   allocatable :: temp_mesh(:)
        type(bcset_t),                  allocatable :: temp_bcset(:)
        type(equationset_wrapper_t),    allocatable :: temp_eqnset(:)


        !
        ! Increment number of domains by one
        !
        self%ndomains = self%ndomains + 1
        idom = self%ndomains


        !
        ! Add (name,index) pair to domain_info dictionary
        !
        call self%domain_info%set(name,idom)


        !
        ! Resize array storage
        !

        ! Allocate new storage arrays
        allocate(temp_mesh(self%ndomains),   &
                 temp_bcset(self%ndomains),  &
                 temp_eqnset(self%ndomains), stat=ierr) 
        if (ierr /= 0) call AllocationError


        ! Copy previously initialized instances to new array
        if (self%ndomains > 1) then
            temp_mesh(  1:size(self%mesh))   = self%mesh
            temp_bcset( 1:size(self%bcset))  = self%bcset
            temp_eqnset(1:size(self%eqnset)) = self%eqnset
        end if


        !
        ! Initialize new mesh
        !
        call temp_mesh(idom)%init_geom(idom,nterms_c,points)


        !
        ! Allocate equation set
        !
        call create_equationset(eqnset,temp_eqnset(idom)%item)


        !
        ! Initialize mesh numerics based on equation set and polynomial expansion order
        !
        call temp_mesh(idom)%init_sol(temp_eqnset(idom)%item%neqns,nterms_s)


        !
        ! Move rezied temp allocation back to chidg_data container
        !
        call move_alloc(temp_mesh,self%mesh)
        call move_alloc(temp_bcset,self%bcset)
        call move_alloc(temp_eqnset,self%eqnset)

    end subroutine add_domain











    subroutine add_bc(self,domain,bc,face,options)
        class(chidg_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: domain
        character(*),           intent(in)      :: bc
        integer(ik),            intent(in)      :: face
        type(dict_t), optional, intent(in)      :: options

        integer(ik)                 :: idom
        class(bc_t), allocatable    :: bc_instance


        !
        ! Get domain index from domain string
        !
        call self%domain_info%get(domain,idom)



        !
        ! Create boundary condition, specified in incoming string
        !
        call create_bc(bc,bc_instance)


        !
        ! Initialize new boundary condition from mesh data and face index
        !
        call bc_instance%init(self%mesh(idom),face,options)


        !
        ! Add initialized boundary condition to bcset_t for domain 'idom'
        !
        call self%bcset(idom)%add(bc_instance)


    end subroutine add_bc











end module type_chidg_data
