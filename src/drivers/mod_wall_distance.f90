module mod_wall_distance
    use mod_constants,  only: ZERO
    use type_chidg,     only: chidg_t
    use type_bc_state,  only: bc_state_t
    use mod_bc,         only: create_bc
    use mpi_f08
    implicit none








contains


    !>  Solve for wall distance using a PDE-based approach.
    !!
    !!  Solve a p-Poisson equation for the approximate distance field. As 'p' goes to 
    !!  infinity, the scalar field satisfies the Eikonal equation, which gives the 
    !!  distance function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine compute_field(field_name,gridfile,solutionfile)
        character(*),   intent(in)      :: field_name
        character(*),   intent(in)      :: gridfile
        character(*),   intent(in)      :: solutionfile

!        type(chidg_t)                   :: chidg
!        type(dict_t)                    :: noptions, loptions
!        class(bc_state_t),  allocatable :: dirichlet_zero, neumann_zero
!
!
!        !
!        ! chidg%init('mpi') should have already been called by another ChiDG instance.
!        ! chidg%init('io')  should have already been called by another ChiDG instance.
!        !
!        ! Make sure this ChiDG environment is initialized.
!        !
!        call chidg%init('env',MPI_COMM_WORLD)
!
!
!
!        !
!        ! Set up boudary condition states to impose.
!        !
!        ! u = 0 on walls
!        !
!        ! dudn = 0 else where
!        !
!        call create_bc('Scalar Value', dirichlet_zero)
!        call create_bc('Scalar Derivative', neumann_zero)
!
!        call dirichlet_zero%set_fcn_option('Value','val', ZERO)
!        call neumann_zero%set_fcn_option('Derivative','val', ZERO)
!
!
!
!
!        !
!        ! Read grid, boundary conditions.
!        !
!        ! Override equation_set and boundary conditions with p-Poisson equation
!        ! for approximate wall distance field.
!        !
!        call chidg%read_grid(gridfile, equation_set='p-Poisson')
!        call chidg%read_boundaryconditions(gridfile,bc_wall=dirichlet_zero,     &
!                                                    bc_inlet=neumann_zero,      &
!                                                    bc_outlet=neumann_zero,     &
!                                                    bc_symmetry=neumann_zero,   &
!                                                    bc_farfield=neumann_zero)
!
!
!        !
!        ! Initialize options dictionaries
!        !
!        call noptions%set('tol',1.e-8_rk)   ! Set nonlinear solver options
!        call noptions%set('cfl0',2.0_rk)
!        call noptions%set('nsteps',50)
!
!        call loptions%set('tol',1.e-10_rk)  ! Set linear solver options
!
!
!
!        !
!        ! Set ChiDG components
!        !
!        call chidg%set('Time Integrator' , algorithm='Steady'                        )
!        call chidg%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
!        call chidg%set('Linear Solver'   , algorithm='FGMRES',       options=loptions)
!        call chidg%set('Preconditioner'  , algorithm='ILU0'                          )
!
!
!
!        !
!        ! Initialize solution
!        !
!        if (solutionfile_in == 'none') then
!
!        else
!
!            ! TODO: put in check that solutionfile actually contains solution
!            call chidg%read_solution(solutionfile_in, 'Wall Distance')
!
!        end if
!
!
!
!        do iorder = 1,order_max
!
!            !
!            ! Initialize domain storage, communication, matrix/vector storage
!            !
!            call chidg%set('Solution Order', integer_input=solution_order)
!            call chidg%initialize_solution_domains()
!            call chidg%init('communication')
!            call chidg%init('chimera')
!            call chidg%initialize_solution_solver()
!
!            !
!            ! Wrap-up initialization activities
!            !
!            call chidg%init('finalize')
!
!            !
!            ! Run ChiDG simulation
!            !
!            call chidg%report('before')
!            call chidg%run()
!            call chidg%report('after')
!        end do
!
!        
!        !
!        ! Write wall distance to auxiliary field
!        !
!        call chidg%write_solution(solutionfile_out)



    end subroutine compute_field
    !**************************************************************************************




end module mod_wall_distance
