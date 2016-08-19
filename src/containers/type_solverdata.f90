module type_solverdata
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NFACES
    use type_chidgVector,               only: chidgVector_t
    use type_chidgMatrix,               only: chidgMatrix_t
    use type_mesh,                      only: mesh_t
    use type_function_status,           only: function_status_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_bcset_coupling,            only: bcset_coupling_t
    use type_BR2,                       only: BR2_t
    use type_element_info,              only: element_info_t
    use type_face_info,                 only: face_info_t
    use type_function_info,             only: function_info_t
    implicit none



    !> Container for solver data.
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   3/15/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!  @note   Added BR2 diffusion container
    !!
    !-------------------------------------------------------------------------------------------------------
    type, public  :: solverdata_t

        !
        ! Base solver data
        !
        type(chidgVector_t)             :: q                        !< Solution vector
        type(chidgVector_t)             :: dq                       !< Change in solution vector
        type(chidgVector_t)             :: rhs                      !< Residual of the spatial scheme
        type(chidgMatrix_t)             :: lhs                      !< Linearization of the spatial scheme


        !
        ! BR2 container
        !
        type(BR2_t)                     :: BR2



        !
        ! Time information
        !
        real(rk)                        :: t                        !< Global time
        real(rk),   allocatable         :: dt(:,:)                  !< Element-local time-step, (ndomains,maxelems)

        !
        ! Function registration
        !
        type(function_status_t)         :: function_status          !< Class for the status of a function residual and linearization


        logical                         :: solverInitialized = .false.




        ! NOTE: if one wanted to add specialized data, instead of deriving from chidgData, maybe you could add a
        !       chidgExtension class that could be specialized further which could contain non-standard data 
        !  class(chidgExtension_t)

    contains

        generic, public     :: init => init_base
        procedure, private  :: init_base

        procedure           :: update_diffusion

    end type solverdata_t
    !*******************************************************************************************************







contains







    !>  Initialize solver base data structures
    !!      - allocate and initialize q, dq, rhs, and linearization.
    !!      - Should be called by specialized 'init' procedure for derived solvers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh                Array of mesh_t instances which define storage requirements.
    !!  @param[in]  bcset_coupling      Array of bcset_coupling instances which describe the coupling of elements in bcs.
    !!  @param[in]  function_data       Array of containers that hold information on number of each function in eqnset.
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine init_base(self,mesh,bcset_coupling,function_data)
        class(solverdata_t),                intent(inout)           :: self
        type(mesh_t),                       intent(inout)           :: mesh(:)
        type(bcset_coupling_t),             intent(in)              :: bcset_coupling(:)
        type(equationset_function_data_t),  intent(in)              :: function_data(:)
        

        integer(ik) :: ierr, ndom, maxelems, idom
        logical     :: increase_maxelems = .false.


        ! Initialize and allocate storage
        call self%q%init(  mesh)
        call self%dq%init( mesh)
        call self%rhs%init(mesh)
        call self%lhs%init(mesh,bcset_coupling,'full')

        ! Initialize matrix parallel recv data
        call self%lhs%init_recv(self%rhs)





    
        !
        ! Find maximum number of elements in any domain
        !
        ndom = size(mesh)
        maxelems = 0
        do idom = 1,ndom

            increase_maxelems = ( mesh(idom)%nelem > maxelems )

            if (increase_maxelems) then
                maxelems = mesh(idom)%nelem
            end if

        end do



        !
        ! Allocate timestep storage
        !
        allocate(self%dt(ndom,maxelems),stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize storage on flux and linearization registration
        !
        call self%function_status%init( mesh, function_data)

        

        !
        ! Confirm solver initialization
        !
        self%solverInitialized = .true.

    end subroutine init_base
    !************************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine update_diffusion(self,mesh,elem_info)
        class(solverdata_t),    intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info


        !
        ! Update BR2 lifing operators
        !
        call self%BR2%update(mesh,elem_info,self%q)


    end subroutine update_diffusion
    !*************************************************************************************************************






















end module type_solverdata
