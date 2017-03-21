module mod_wall_distance
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    use eqn_wall_distance,  only: set_p_poisson_parameter
    use mod_function,       only: create_function
    use type_function,      only: function_t
    use type_chidg,         only: chidg_t
    use type_bc_state,      only: bc_state_t
    use type_dict,          only: dict_t
    use mod_chidg_post,     only: chidg_post,chidg_post_vtk
    use mod_bc,             only: create_bc
    use mod_io,             only: gridfile
    use mod_hdf_utilities,  only: get_properties_hdf, check_file_exists_hdf
    use mod_chidg_mpi,      only: IRANK, NRANK, ChiDG_COMM
    use type_file_properties,   only: file_properties_t
    implicit none





contains


    !>  Solve for wall distance using a PDE-based approach.
    !!
    !!  Solve a p-Poisson equation for the approximate distance field. As 'p' goes to 
    !!  infinity, the scalar field satisfies the Eikonal equation, which gives the 
    !!  distance function.
    !!
    !!  @date   11/15/2016
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  gridfile        ChiDG-file containing a grid to be read and computed on.
    !!  @param[in]  solutionfile    ChiDG-file the field will be written to.
    !!  @param[in]  order           Polynomial order, the field will be computed with.
    !!
    !-------------------------------------------------------------------------------------
    subroutine wall_distance_driver(chidg,fileout)
        type(chidg_t),  intent(inout)           :: chidg
        character(*),   intent(in), optional    :: fileout

        character(:), allocatable   :: user_msg
        integer(ik)                 :: order

        type(chidg_t)                   :: wall_distance
        type(dict_t)                    :: noptions, loptions
        class(bc_state_t),  allocatable :: dirichlet_zero, neumann_zero
        class(function_t),  allocatable :: constant
        integer(ik)                     :: iorder, p, aux_field_index, wd_nterms_s, ierr, iproc
        type(file_properties_t)         :: wd_props
        logical                         :: wd_file_exists, have_wd_field



        !
        ! chidg%init('mpi') should have already been called by another ChiDG instance.
        ! chidg%init('io')  should have already been called by another ChiDG instance.
        !
        ! Make sure this ChiDG environment is initialized.
        !
        call wall_distance%start_up('core')




        !
        ! Initialize options dictionaries
        !
        call noptions%set('tol',1.e-8_rk)   ! Set nonlinear solver options
        call noptions%set('cfl0',1.0_rk)
        call noptions%set('nsteps',50)
        call loptions%set('tol',1.e-10_rk)  ! Set linear solver options



        !
        ! Set ChiDG components
        !
        call wall_distance%set('Time Integrator' , algorithm='Steady'                        )
        call wall_distance%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
        call wall_distance%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)
        call wall_distance%set('Preconditioner'  , algorithm='RASILU0'                       )

        order = chidg%nterms_s_1d
        call wall_distance%set('Solution Order', integer_input=order)


        !
        ! Set up boudary condition states to impose.
        !
        ! u = 0 on walls
        !
        ! dudn = 0 else where
        !
        call create_bc('Scalar Value',      dirichlet_zero)
        call create_bc('Scalar Derivative', neumann_zero  )

        call dirichlet_zero%set_fcn_option( 'Value'     , 'val', ZERO)
        call neumann_zero%set_fcn_option(   'Derivative', 'val', ZERO)




        !
        ! Read grid, boundary conditions.
        !
        ! Override equation_set and boundary conditions with p-Poisson equation
        ! for approximate wall distance field.
        !
        ! Solid walls get dirichlet zero bc.
        ! All other families get neumann zero bc.
        !
        call wall_distance%read_grid(gridfile, equation_set='Wall Distance : p-Poisson')
        call wall_distance%read_boundaryconditions(gridfile, bc_wall     = dirichlet_zero, &
                                                             bc_inlet    = neumann_zero,   &
                                                             bc_outlet   = neumann_zero,   &
                                                             bc_symmetry = neumann_zero,   &
                                                             bc_farfield = neumann_zero)



        !
        ! Initialize wall_distance with chidg order in case we are going to read in a solution.
        !
        call wall_distance%init('all')




        !
        ! Check if we already have a wall distance solution
        !
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then

                !wd_file_exists = check_file_exists_hdf('wall_distance.h5')
                wd_file_exists = check_file_exists_hdf(fileout)
                if (wd_file_exists) then
                    !wd_props = get_properties_hdf('wall_distance.h5')
                    wd_props = get_properties_hdf(fileout)
                    wd_nterms_s = wd_props%nterms_s(1)

                    have_wd_field = (wd_nterms_s >= chidg%nterms_s)
                end if

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do




        !
        ! If we have a wall distance file and it has an accurate solution,
        ! just read from file.
        !
        if (wd_file_exists .and. have_wd_field) then

            !call wall_distance%read_solution('wall_distance.h5')
            call wall_distance%read_solution(fileout)
            wall_distance%data%sdata%q = wall_distance%data%sdata%q_in

        !
        ! If we don't have an accurate wall distance field in file, solve for a new one.
        !
        else

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
            call noptions%set('tol',1.e-4_rk)   ! Set nonlinear solver options
            call loptions%set('tol',1.e-5_rk)   ! Set linear solver options
            call wall_distance%set('Nonlinear Solver', algorithm='Newton', options=noptions)
            call wall_distance%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)
            do p = 2,6,2
                call write_line('Wall Distance Driver : Loop 1 : p = ', p)
                
                !
                ! Update p-Poisson fidelity
                !
                call set_p_poisson_parameter(real(p,rk))


                !
                ! (Re)Initialize domain storage, communication, matrix/vector storage
                !
                call wall_distance%set('Solution Order', integer_input=iorder)
                call wall_distance%init('all')


                !
                ! Read solution if it exists.
                !
                if (p == 2) then
                    call create_function(constant,'constant')
                    call constant%set_option('val',0.1_rk)
                    call wall_distance%data%sdata%q_in%project(wall_distance%data%mesh,constant,1)

                else
                    call wall_distance%read_solution(fileout)
                end if


                !
                ! Run ChiDG simulation
                !
                call wall_distance%report('before')
                call wall_distance%run(write_initial=.false., write_final=.false.)
                call wall_distance%report('after')


                !
                ! Write wall distance to auxiliary field
                !
                call wall_distance%write_solution(fileout)


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

            !
            ! Update p-Poisson fidelity
            !
            p = 6
            call set_p_poisson_parameter(real(p,rk))
            call noptions%set('tol',1.e-5_rk)   ! Set nonlinear solver options
            call loptions%set('tol',1.e-8_rk)   ! Set linear solver options
            call wall_distance%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
            call wall_distance%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)

            order = chidg%nterms_s_1d
            do iorder = 3,order
                call write_line('Wall Distance Driver : Loop 2 : order = ', iorder)


                !
                ! (Re)Initialize domain storage, communication, matrix/vector storage
                !
                call wall_distance%set('Solution Order', integer_input=iorder)
                call wall_distance%init('all')


                !
                ! Read solution if it exists.
                !
                call wall_distance%read_solution(fileout)



                !
                ! Run ChiDG simulation
                !
                call wall_distance%report('before')
                call wall_distance%run(write_initial=.false., write_final=.false.)
                call wall_distance%report('after')


                !
                ! Write wall distance to auxiliary field
                !
                call wall_distance%write_solution(fileout)


            end do


        end if ! have_wd_field .and. wd_file_exists

        call write_line('Storing Wall Distance field to Auxiliary field ChiDG Vector:', io_proc=GLOBAL_MASTER)


        !
        ! Try to find 'Wall Distance' auxiliary field storage.
        !
        aux_field_index = chidg%data%sdata%get_auxiliary_field_index('Wall Distance : p-Poisson')




        !
        ! If no 'Wall Distance' auxiliary field storage was not found, create one.
        !
        if (aux_field_index == 0) then

            call chidg%data%sdata%add_auxiliary_field('Wall Distance : p-Poisson', wall_distance%data%sdata%q)

        !
        ! If 'Wall Distance' auxiliary field storage was found, copy Wall Distance solution 
        ! to working ChiDG environment.
        !
        else
            chidg%data%sdata%auxiliary_field(aux_field_index) = wall_distance%data%sdata%q

        end if



    end subroutine wall_distance_driver
    !**************************************************************************************




end module mod_wall_distance
