module mod_wall_distance
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO
    use type_chidg,     only: chidg_t
    use type_bc_state,  only: bc_state_t
    use type_dict,      only: dict_t
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
    !!  @param[in]  gridfile        ChiDG-file containing a grid to be read and computed on.
    !!  @param[in]  solutionfile    ChiDG-file the field will be written to.
    !!  @param[in]  order           Polynomial order, the field will be computed with.
    !!
    !-------------------------------------------------------------------------------------
    subroutine compute_field(field_name,gridfile,solutionfile,order)
        character(*),   intent(in)      :: field_name
        character(*),   intent(in)      :: gridfile
        character(*),   intent(in)      :: solutionfile
        integer(ik),    intent(in)      :: order

        type(chidg_t)                   :: chidg
        type(dict_t)                    :: noptions, loptions
        class(bc_state_t),  allocatable :: dirichlet_zero, neumann_zero
        integer(ik)                     :: iorder, p


        !
        ! chidg%init('mpi') should have already been called by another ChiDG instance.
        ! chidg%init('io')  should have already been called by another ChiDG instance.
        !
        ! Make sure this ChiDG environment is initialized.
        !
        call chidg%init('env',MPI_COMM_WORLD)



        !
        ! Set up boudary condition states to impose.
        !
        ! u = 0 on walls
        !
        ! dudn = 0 else where
        !
        call create_bc('Scalar Value', dirichlet_zero)
        call create_bc('Scalar Derivative', neumann_zero)

        call dirichlet_zero%set_fcn_option( 'Value',     'val', ZERO)
        call neumann_zero%set_fcn_option(   'Derivative','val', ZERO)




        !
        ! Read grid, boundary conditions.
        !
        ! Override equation_set and boundary conditions with p-Poisson equation
        ! for approximate wall distance field.
        !
        ! Solid walls get dirichlet zero bc.
        ! All other families get neumann zero bc.
        !
        call chidg%read_grid(gridfile, equation_set='p-Poisson')
        call chidg%read_boundaryconditions(gridfile,bc_wall=dirichlet_zero,     &
                                                    bc_inlet=neumann_zero,      &
                                                    bc_outlet=neumann_zero,     &
                                                    bc_symmetry=neumann_zero,   &
                                                    bc_farfield=neumann_zero)


        !
        ! Initialize options dictionaries
        !
        call noptions%set('tol',1.e-8_rk)   ! Set nonlinear solver options
        call noptions%set('cfl0',2.0_rk)
        call noptions%set('nsteps',50)
        call loptions%set('tol',1.e-10_rk)  ! Set linear solver options



        !
        ! Set ChiDG components
        !
        call chidg%set('Time Integrator' , algorithm='Steady'                        )
        call chidg%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
        call chidg%set('Linear Solver'   , algorithm='FGMRES',       options=loptions)
        call chidg%set('Preconditioner'  , algorithm='ILU0'                          )




        ! Get wall-distance approximation for p-Poisson equation using a low-order
        ! polynomial expansion. We are going in steps of 'p' here to make sure
        ! we get good convergence of the Newton solver by having a good initial
        ! solution from the previous calculation.
        !
        !   Polynomial Order:
        !       P = 2 
        !
        !   p-Poisson Parameter:
        !       p = 2,4,6
        !
        iorder = 2
        do p = 2,2,6

            !
            ! (Re)Initialize domain storage, communication, matrix/vector storage
            !
            call chidg%set('Solution Order', integer_input=iorder)
            call chidg%initialize_solution_domains()
            call chidg%init('communication')
            call chidg%init('chimera')
            call chidg%initialize_solution_solver()
            call chidg%init('finalize')





            !
            ! Read solution if it exists.
            !
            call chidg%read_solution(solutionfile)


            !
            ! Run ChiDG simulation
            !
            call chidg%report('before')
            call chidg%run()
            call chidg%report('after')


            !
            ! Write wall distance to auxiliary field
            !
            call chidg%write_solution(solutionfile)


        end do





        ! Get wall-distance approximation for p-Poisson equation using a high 'p'
        ! and increasing polynomial expansion. We are going in steps of 'P' here to 
        ! make sure we get good convergence of the Newton solver by having a good initial
        ! solution from the previous calculation.
        !
        !   Polynomial Order:
        !       P = 2,3, ...
        !
        !   p-Poisson Parameter:
        !       p = 6
        !
        p = 6
        do iorder = 2,order

            !
            ! (Re)Initialize domain storage, communication, matrix/vector storage
            !
            call chidg%set('Solution Order', integer_input=iorder)
            call chidg%initialize_solution_domains()
            call chidg%init('communication')
            call chidg%init('chimera')
            call chidg%initialize_solution_solver()
            call chidg%init('finalize')





            !
            ! Read solution if it exists.
            !
            call chidg%read_solution(solutionfile)


            !
            ! Run ChiDG simulation
            !
            call chidg%report('before')
            call chidg%run()
            call chidg%report('after')


            !
            ! Write wall distance to auxiliary field
            !
            call chidg%write_solution(solutionfile)


        end do





        
        !
        ! Project normalization for better approximation of distance function
        !
!        call normalize_wall_distance(chidg)


        !
        ! Write wall distance to auxiliary field
        !
        call chidg%write_solution(solutionfile)



    end subroutine compute_field
    !**************************************************************************************




end module mod_wall_distance
