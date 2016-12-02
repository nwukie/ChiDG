module SD_bc_operator
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: SD_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Scalar Diffusion BC Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Diffusive Operator")

        !
        ! Set operator equations
        !
        call self%add_primary_field("u")

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SD_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: iu


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            u, mu, dudx, dudy, dudz, flux_x, flux_y, flux_z, integrand


        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz


        !
        ! Get equation indices
        !
        iu = prop%get_primary_field_index("u")



        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        u     = worker%get_primary_field_face(iu, 'value'     , 'boundary')
        dudx  = worker%get_primary_field_face(iu, 'ddx + lift', 'boundary')
        dudy  = worker%get_primary_field_face(iu, 'ddy + lift', 'boundary')
        dudz  = worker%get_primary_field_face(iu, 'ddz + lift', 'boundary')


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Compute scalar coefficient
        !
        mu = prop%scalar%compute_mu(u,dudx,dudy,dudz)

        

        !=================================================
        ! Mass flux
        !=================================================
        flux_x = -mu*dudx
        flux_y = -mu*dudy
        flux_z = -mu*dudz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(iu, integrand)


    end subroutine compute
    !**********************************************************************************************










end module SD_bc_operator
