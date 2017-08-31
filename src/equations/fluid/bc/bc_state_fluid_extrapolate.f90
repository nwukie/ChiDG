module bc_state_fluid_extrapolate
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: fluid_extrapolate_t

    contains

        procedure   :: init                 !< Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     !< boundary condition function implementation

    end type fluid_extrapolate_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_extrapolate_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Fluid Extrapolate")
        call self%set_family("Symmetry")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("Momentum-1")
!        call self%set_equation("Momentum-2")
!        call self%set_equation("Momentum-3")
!        call self%set_equation("Energy"    )


        !
        ! Add functions
        !


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
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(fluid_extrapolate_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,                            &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m



        !
        ! Interpolate interior solution to face quadrature nodes
        !
      
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')

        grad1_density_m = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m    = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field('Momentum-3', 'grad3', 'face interior')

        grad1_energy_m  = worker%get_field('Energy'    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_field('Energy'    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_field('Energy'    , 'grad3', 'face interior')





        !
        ! Store boundary condition state
        !
        call worker%store_bc_state("Density"   , density_m, 'value')
        call worker%store_bc_state("Momentum-1", mom1_m,    'value')
        call worker%store_bc_state("Momentum-2", mom2_m,    'value')
        call worker%store_bc_state("Momentum-3", mom3_m,    'value')
        call worker%store_bc_state("Energy"    , energy_m,  'value')




        call worker%store_bc_state("Density"   , grad1_density_m, 'grad1')
        call worker%store_bc_state("Density"   , grad2_density_m, 'grad2')
        call worker%store_bc_state("Density"   , grad3_density_m, 'grad3')

        call worker%store_bc_state("Momentum-1", grad1_mom1_m,    'grad1')
        call worker%store_bc_state("Momentum-1", grad2_mom1_m,    'grad2')
        call worker%store_bc_state("Momentum-1", grad3_mom1_m,    'grad3')

        call worker%store_bc_state("Momentum-2", grad1_mom2_m,    'grad1')
        call worker%store_bc_state("Momentum-2", grad2_mom2_m,    'grad2')
        call worker%store_bc_state("Momentum-2", grad3_mom2_m,    'grad3')

        call worker%store_bc_state("Momentum-3", grad1_mom3_m,    'grad1')
        call worker%store_bc_state("Momentum-3", grad2_mom3_m,    'grad2')
        call worker%store_bc_state("Momentum-3", grad3_mom3_m,    'grad3')

        call worker%store_bc_state("Energy"    , grad1_energy_m,  'grad1')
        call worker%store_bc_state("Energy"    , grad2_energy_m,  'grad2')
        call worker%store_bc_state("Energy"    , grad3_energy_m,  'grad3')


    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_fluid_extrapolate
