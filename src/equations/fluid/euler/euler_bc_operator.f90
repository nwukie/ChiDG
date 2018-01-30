module euler_bc_operator
    use mod_constants,      only: HALF,ONE, TWO
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: euler_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

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
        class(euler_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            density_bc,  mom1_bc, mom2_bc, mom3_bc, energy_bc,  &
            u_bc,    v_bc,    w_bc,                             &
            H_bc,    p_bc,                                      &
            flux_1,  flux_2,  flux_3

        real(rk),   allocatable, dimension(:) :: r


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        density_bc = worker%get_field("Density"   , 'value', 'boundary')
        mom1_bc    = worker%get_field("Momentum-1", 'value', 'boundary')
        mom2_bc    = worker%get_field("Momentum-2", 'value', 'boundary')
        mom3_bc    = worker%get_field("Momentum-3", 'value', 'boundary')
        energy_bc  = worker%get_field("Energy"    , 'value', 'boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary') 
            mom2_bc = mom2_bc / r
        end if


        !
        ! Compute gamma
        !
        p_bc = worker%get_field('Pressure', 'value', 'boundary')


        !
        ! Compute velocity components
        !
        u_bc = mom1_bc/density_bc
        v_bc = mom2_bc/density_bc
        w_bc = mom3_bc/density_bc



        !
        ! Compute boundary condition energy and enthalpy
        !
        H_bc = (energy_bc + p_bc)/density_bc




        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = (density_bc * u_bc)
        flux_2 = (density_bc * v_bc)
        flux_3 = (density_bc * w_bc)

        call worker%integrate_boundary_condition('Density','Advection',flux_1,flux_2,flux_3)

        !=================================================
        ! Momenum-1 flux
        !=================================================
        flux_1 = (mom1_bc * u_bc) + p_bc
        flux_2 = (mom1_bc * v_bc) 
        flux_3 = (mom1_bc * w_bc) 

        call worker%integrate_boundary_condition('Momentum-1','Advection',flux_1,flux_2,flux_3)

        !=================================================
        ! Momentum-2 flux
        !=================================================
        flux_1 = (mom2_bc * u_bc) 
        flux_2 = (mom2_bc * v_bc) + p_bc
        flux_3 = (mom2_bc * w_bc) 

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if

        call worker%integrate_boundary_condition('Momentum-2','Advection',flux_1,flux_2,flux_3)

        !=================================================
        ! Momentum-3 flux
        !=================================================
        flux_1 = (mom3_bc * u_bc) 
        flux_2 = (mom3_bc * v_bc) 
        flux_3 = (mom3_bc * w_bc) + p_bc

        call worker%integrate_boundary_condition('Momentum-3','Advection',flux_1,flux_2,flux_3)

        !=================================================
        ! Energy flux
        !=================================================
        flux_1 = (density_bc * H_bc * u_bc) 
        flux_2 = (density_bc * H_bc * v_bc) 
        flux_3 = (density_bc * H_bc * w_bc) 

        call worker%integrate_boundary_condition('Energy','Advection',flux_1,flux_2,flux_3)

    end subroutine compute
    !*******************************************************************************************























end module euler_bc_operator
