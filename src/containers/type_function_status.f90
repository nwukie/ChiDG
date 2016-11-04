module type_function_status
#include <messenger.h>
    use mod_constants,                  only: BOUNDARY_ADVECTIVE_FLUX
    use type_mesh,                      only: mesh_t
    use type_equationset_function_data, only: equationset_function_data_t
    use type_function_status_data,      only: function_status_data_t
    use type_face_info,                 only: face_info_t
    use type_function_info,             only: function_info_t
    implicit none




    !> Class is used to manage the contributions and linearizations of functions.
    !!
    !! Domain component contains arrays of logicals that indicate the status of a given function. 
    !!
    !! So, in this way, we can compute a flux on a face and add the contribution to both elements
    !! by registering that the function contribution was added. This container holds the registration
    !! information as arrays of logicals.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: function_status_t

       type(function_status_data_t), allocatable :: dom(:)

    contains

        procedure   :: init

        procedure   :: compute_function             !< Check if a function has been executed
        procedure   :: linearize_function           !< Check if a function has been linearized

        procedure   :: compute_function_equation    !< Check whether a function has contributed to a particular equation
        procedure   :: linearize_function_equation  !< Check whether a function contribution to a particualar equation has been linearized


        procedure   :: register_function_computed   !< Register a function as computed
        procedure   :: register_function_linearized !< Register a function as linearized

        procedure   :: clear                        !< Reset function status data

    end type function_status_t
    !*****************************************************************************************




contains



    !>  Initialize class storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh            Array of mesh_t instances.
    !!  @param[in]  function_data   Contains the number of each function type for every domain
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self, mesh, function_data)
        class(function_status_t),           intent(inout)   :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(equationset_function_data_t),  intent(in)      :: function_data(:)

        integer :: idom, ndom, ierr
        
        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%dom)) deallocate(self%dom)


        ndom = size(mesh)
        allocate(self%dom(ndom), stat=ierr)
        if (ierr /= 0) call AllocationError


        do idom = 1,ndom
            call self%dom(idom)%init( mesh(idom), function_data(idom) )
        end do



    end subroutine init
    !*****************************************************************************************









    !>  Call to determine if a function needs computed for a given face, element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info   Container for face information. Location in mesh, type, etc.
    !!  @param[in]  flux_info   Container for function information. iblk, idonor, iflux, type
    !!  @return     res         True if function needs computed. False if already computed.
    !!
    !-----------------------------------------------------------------------------------------
    function compute_function(self, face_info, function_info) result(res)
        class(function_status_t),   intent(in)              :: self
        type(face_info_t),          intent(in)              :: face_info
        type(function_info_t),      intent(in)              :: function_info
        
        logical :: function_already_computed
        logical :: res



        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn )


        function_already_computed = self%dom(idomain_l)%function_computed(type,ielement_l,iface,ifcn)


        if ( function_already_computed ) then
            res = .false.
        else
            res = .true.
        end if

        end associate

    end function compute_function
    !*****************************************************************************************














    !>  Call to determine if a particular equation has been contributed to from a given function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info   Container for face information. Location in mesh, type, etc.
    !!  @param[in]  flux_info   Container for function information. iblk, idonor, iflux, type
    !!  @param[in]  ieqn        Equation index within the function definition
    !!  @return     res         True if function needs computed. False if already computed.
    !!
    !-----------------------------------------------------------------------------------------
    function compute_function_equation(self, face_info, function_info, ieqn) result(res)
        class(function_status_t),   intent(in)              :: self
        type(face_info_t),          intent(in)              :: face_info
        type(function_info_t),      intent(in)              :: function_info
        integer(ik),                intent(in)              :: ieqn
        
        logical :: function_already_computed
        logical :: res



        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn,  iblk => function_info%idiff )


        function_already_computed = self%dom(idomain_l)%function_equation_computed(type,ielement_l,iface,ifcn,ieqn)


        if ( function_already_computed ) then
            res = .false.
        else
            res = .true.
        end if

        end associate

    end function compute_function_equation
    !*****************************************************************************************












    !>  Call to determine if a function needs linearized for a given face, element, seed.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info   Container for face information. Location in mesh, type, etc.
    !!  @param[in]  flux_info   Container for function information. iblk, idonor, iflux, type
    !!  @return     res         True if function needs linearized. False if already linearized.
    !!
    !-----------------------------------------------------------------------------------------
    function linearize_function(self, face_info, function_info) result(res)
        class(function_status_t),   intent(in)  :: self
        type(face_info_t),          intent(in)  :: face_info
        type(function_info_t),      intent(in)  :: function_info
        
        logical :: function_already_linearized
        logical :: res


        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn,  iblk => function_info%idiff )


        function_already_linearized = self%dom(idomain_l)%function_linearized(type,ielement_l,iface,ifcn,iblk)



        if ( function_already_linearized ) then
            res = .false.
        else
            res = .true.
        end if

        end associate

    end function linearize_function
    !*****************************************************************************************












    !>  Call to determine if a particular equation has been linearized from a given function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info   Container for face information. Location in mesh, type, etc.
    !!  @param[in]  flux_info   Container for function information. iblk, idonor, iflux, type
    !!  @return     res         True if function needs linearized. False if already linearized.
    !!
    !-----------------------------------------------------------------------------------------
    function linearize_function_equation(self, face_info, function_info, ieqn) result(res)
        class(function_status_t),   intent(in)  :: self
        type(face_info_t),          intent(in)  :: face_info
        type(function_info_t),      intent(in)  :: function_info
        integer(ik),                intent(in)  :: ieqn
        
        logical :: function_already_linearized
        logical :: res


        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn,  iblk => function_info%idiff )


        function_already_linearized = self%dom(idomain_l)%function_equation_linearized(type,ielement_l,iface,ifcn,iblk,ieqn)





        if ( function_already_linearized ) then
            res = .false.
        else
            res = .true.
        end if

        end associate

    end function linearize_function_equation
    !*****************************************************************************************












    !>  Register a function as having been computed for a given element/face.
    !!
    !!  This sets both the 'function_computed' and 'function_equation_computed' logicals.
    !!  So, in this way, because this gets called for every equation, the 'function_computed' value
    !!  will get set multiple times. Once for every equation. That is okay. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info       Contains face indices
    !!  @param[in]  function_info   Contains function indices, type, linearization index
    !!  @param[in]  ieqn            Index of equation being registered
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine register_function_computed(self,face_info,function_info,ieqn)
        class(function_status_t),   intent(inout)   :: self
        type(face_info_t),          intent(in)      :: face_info
        type(function_info_t),      intent(in)      :: function_info
        integer(ik),                intent(in)      :: ieqn
        


        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn,  iblk => function_info%idiff )


            self%dom(idomain_l)%function_computed(type,ielement_l,iface,ifcn) = .true.
            self%dom(idomain_l)%function_equation_computed(type,ielement_l,iface,ifcn,ieqn) = .true.

        end associate
    end subroutine register_function_computed
    !******************************************************************************************












    !>  Register a function as having been linearized for a given element/face/blk.
    !!
    !!  This sets both the 'function_linearized' and 'function_equation_linearized' logicals.
    !!  So, in this way, because this gets called for every equation, the 'function_linearized' value
    !!  will get set multiple times. Once for every equation. That is okay. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  face_info       Contains face indices
    !!  @param[in]  function_info   Contains function indices, type, linearization index
    !!  @param[in]  ieqn            Index of equation being registered
    !!
    !------------------------------------------------------------------------------------------
    subroutine register_function_linearized(self,face_info,function_info,ieqn)
        class(function_status_t),   intent(inout)   :: self
        type(face_info_t),          intent(in)      :: face_info
        type(function_info_t),      intent(in)      :: function_info
        integer(ik),                intent(in)      :: ieqn
        


        associate(idomain_l => face_info%idomain_l,  ielement_l => face_info%ielement_l, iface => face_info%iface, &
                  type => function_info%type, ifcn => function_info%ifcn,  iblk => function_info%idiff )


            self%dom(idomain_l)%function_linearized(type,ielement_l,iface,ifcn,iblk) = .true.
            self%dom(idomain_l)%function_equation_linearized(type,ielement_l,iface,ifcn,iblk,ieqn) = .true.

        end associate
    end subroutine register_function_linearized
    !*****************************************************************************************











    !>  Call clear on each domain local function status data. Should be called before every
    !!  evaluation of the right-hand side.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine clear(self)
        class(function_status_t),   intent(inout)   :: self

        integer :: idom


        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do

    end subroutine clear
    !****************************************************************************************








end module type_function_status
