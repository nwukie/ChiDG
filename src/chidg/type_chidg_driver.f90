!>  A submodule of type_chidg providing the 'auxiliary_driver' subroutine
!!  to drive chidg auxiliary_environment instances.
!!
!!  If you would like to provide another auxiliary driver routine: 
!!  --------------------------------------------------------------
!!      1: create a few file that is a submodule of type_chidg_driver.
!!         type_chidg_driver is itself a submodule of type_chidg so the
!!         first line should look something like:
!!  
!!         submodule (type_chidg:type_chidg_driver) my_driver
!!
!!      2: Make sure your routine in the submodule is declared as a
!!         'module subroutine'
!!
!!      3: In this file, declare your driver interface under the 
!!         "AUXILIARY DRIVER SUBMODULE INTERFACES" text.
!!
!!      4: Register your new driver under the 'auxiliary_driver' routine
!!         in this file.
!!
!!      NOTE: you can use wall_distance_driver as a template
!!
!!
!!  @author Nathan A. Wukie
!!  @date   6/25/2017
!!
!--------------------------------------------------------------------------------
submodule (type_chidg) type_chidg_driver
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO
    use mod_gridspace,          only: linspace
    use eqn_wall_distance,      only: set_p_poisson_parameter, p_max, p_min, p_increment, p_sub_rtol, p_sub_nmax
    use mod_function,           only: create_function
    use type_function,          only: function_t
    use type_bc_state,          only: bc_state_t
    use type_dict,              only: dict_t
    use mod_bc,                 only: create_bc
    use mod_hdf_utilities,      only: get_properties_hdf, check_file_exists_hdf
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_file_properties,   only: file_properties_t
    use precon_RASILU0,         only: precon_RASILU0_t
    implicit none


    !!-----------------------------------------------
    !!    AUXILIARY DRIVER SUBMODULE INTERFACES
    !!-----------------------------------------------
    !
    !! Provided by driver_wall_distance.f90 submodule
    !interface
    !    module subroutine wall_distance_driver(chidg,wall_distance,grid_file,aux_file)
    !        type(chidg_t),  intent(inout)   :: chidg
    !        type(chidg_t),  intent(inout)   :: wall_distance
    !        character(*),   intent(in)      :: grid_file
    !        character(*),   intent(in)      :: aux_file
    !    end subroutine wall_distance_driver
    !end interface



contains


    !>  An interface for calling auxiliary driver routines.
    !!
    !!  Called in chidg%prerun()
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2017
    !!
    !----------------------------------------------------------------------------------
    module subroutine auxiliary_driver(chidg,chidg_aux,case,solver_type,grid_file,aux_file)
        type(chidg_t),  intent(inout)   :: chidg
        type(chidg_t),  intent(inout)   :: chidg_aux
        character(*),   intent(in)      :: case
        character(*),   intent(in)      :: solver_type
        character(*),   intent(in)      :: grid_file
        character(*),   intent(in)      :: aux_file

        select case(trim(case))
            case('Wall Distance')

                select case(trim(solver_type))
                    case ('primal problem', 'primal solver', 'primal')
                        call wall_distance_primal_driver(chidg         = chidg,        &
                                                         wall_distance = chidg_aux,    &
                                                         grid_file     = grid_file,    &
                                                         aux_file      = aux_file)

                    case default
                       call chidg_signal_one(FATAL,"chidg_driver_t%auxiliary_driver: solver type found.", trim(solver_type)) 
                end select !Solver type

        end select !Problem type

    end subroutine auxiliary_driver
    !**********************************************************************************









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
    subroutine wall_distance_primal_driver(chidg,wall_distance,grid_file,aux_file)
        type(chidg_t),  intent(inout)   :: chidg
        type(chidg_t),  intent(inout)   :: wall_distance
        character(*),   intent(in)      :: grid_file
        character(*),   intent(in)      :: aux_file

        character(:), allocatable   :: user_msg
        integer(ik)                 :: order

        type(dict_t)                    :: noptions, loptions
        class(bc_state_t),  allocatable :: dirichlet_zero, neumann_zero
        class(function_t),  allocatable :: constant
        integer(ik)                     :: iorder, p_index, aux_field_index, wd_nterms_s, ierr, iproc, np_values
        real(rk), allocatable           :: p_values(:)
        real(rk)                        :: p
        type(file_properties_t)         :: wd_props
        logical                         :: wd_file_exists, have_wd_field
        type(functional_group_t)        :: functionals


        ! Make sure this ChiDG environment is initialized.
        call wall_distance%start_up('core')

        ! Initialize solver dictionaries
        call initialize_input_dictionaries(noptions,loptions)
        call noptions%set('tol',1.e-8_rk)
        call noptions%set('cfl0',1.0_rk)
        call noptions%set('nsteps',50)
        call noptions%set('search','Backtrack')
        call noptions%set('ptc',.false.)
        call noptions%set('smooth',.false.)

        call loptions%set('tol',1.e-10_rk)  
        call loptions%set('inner_fgmres',.true.)


        ! Set ChiDG components
        call wall_distance%set('Time Integrator' , algorithm='Steady')
        call wall_distance%set('Nonlinear Solver', algorithm='Newton', options=noptions)
        call wall_distance%set('Linear Solver'   , algorithm='fgmres', options=loptions)
        call wall_distance%set('Preconditioner'  , algorithm='RASILU0')


        ! Set RASILU settings for petsc
        select type (precon => wall_distance%preconditioner)
            type is (precon_RASILU0_t)
                precon%asm_overlap = 6
                precon%ilu_levels = 3
        end select


        order = chidg%nterms_s_1d
        call wall_distance%set('Solution Order', integer_input=order)



        ! Set up boudary condition states to impose.
        !   u    = 0      on walls
        !   dudn = 0      else where
        call create_bc('Scalar Value',      dirichlet_zero)
        call create_bc('Scalar Derivative', neumann_zero  )

        call dirichlet_zero%set_fcn_option('Value', 'val', ZERO)
        call neumann_zero%set_fcn_option('Normal Gradient', 'val', ZERO)


        ! Read grid, boundary conditions.
        !
        ! Override equation_set and boundary conditions with p-Poisson equation
        ! for approximate wall distance field.
        !
        ! Solid walls get dirichlet zero bc.
        ! All other families get neumann zero bc.
        call wall_distance%read_mesh(grid_file, storage      = 'primal storage',             &
                                                equation_set = 'Wall Distance : p-Poisson',  &
                                                bc_wall      = dirichlet_zero,               &
                                                bc_inlet     = neumann_zero,                 &
                                                bc_outlet    = neumann_zero,                 &
                                                bc_symmetry  = neumann_zero,                 &
                                                bc_farfield  = neumann_zero,                 &
                                                functionals  = functionals)


        ! Initialize wall_distance with chidg order in case we are going to read in a solution.
        call wall_distance%init('all')


        ! Check if we already have a wall distance solution
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then

                wd_file_exists = check_file_exists_hdf(aux_file)
                if (wd_file_exists) then
                    wd_props = get_properties_hdf(aux_file)
                    wd_nterms_s = wd_props%nterms_s(1)

                    have_wd_field = (wd_nterms_s >= chidg%nterms_s)
                end if

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do




        ! If we have a wall distance file and it has an accurate solution,
        ! just read from file.
        if (wd_file_exists .and. have_wd_field) then

            call wall_distance%read_fields(aux_file)
            associate (q_in => wall_distance%data%sdata%q_in, q => wall_distance%data%sdata%q)
                !wall_distance%data%sdata%q = wall_distance%data%sdata%q_in     ! For some reason, this causes a segfault, but is okay when wrapped inside associate
                q = q_in
            end associate


        ! If we don't have an accurate wall distance field in file, solve for a new one.
        else

            ! Store grid to file
            call wall_distance%write_mesh(aux_file)

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

            np_values = int(((p_max-p_min)/p_increment) + 1.)
            p_values = linspace(p_min,p_max, np_values)
            if (size(p_values) /= np_values) call chidg_signal(FATAL,"wall_distance_driver: np_values does not equal the size of actual array of p-values.")

            ! (Re)Initialize domain storage, communication, matrix/vector storage
            iorder = 2
            call wall_distance%set('Solution Order', integer_input=iorder)
            call wall_distance%init('all')

            do p_index = 1,np_values

                ! Update p-Poisson fidelity
                p = p_values(p_index)
                call set_p_poisson_parameter(p)
                call write_line('Wall Distance Driver : Loop 1 : p = ', p)
                
                if (p < p_max) then
                    call noptions%set('tol',1.e-3_rk)   
                    call noptions%set('rtol',p_sub_rtol)   
                    call noptions%set('ptc',.false.)    
                    call noptions%set('smooth',.false.)
                    call noptions%set('nmax',p_sub_nmax)
                    call loptions%set('tol',1.e-5_rk)   
                    call loptions%set('rtol',1.e-16_rk) 
                else
                    call noptions%set('tol',1.e-4_rk)   
                    call noptions%set('rtol',1.e-30_rk)   
                    call noptions%set('ptc',.false.)    
                    call noptions%set('smooth',.false.)
                    call noptions%set('nmax',-1)
                    call loptions%set('tol',1.e-5_rk)   
                    call loptions%set('rtol',1.e-16_rk) 
                end if

                call wall_distance%set('Nonlinear Solver', algorithm='Newton', options=noptions)
                call wall_distance%set('Linear Solver'   , algorithm='fgmres', options=loptions)


                ! Read solution if it exists.
                if (abs(p-TWO) < 1.e-6_rk) then
                    call create_function(constant,'constant')
                    call constant%set_option('val',0.001_rk)
                    !call create_function(constant,'radius')
                    call wall_distance%data%sdata%q_in%project(wall_distance%data%mesh,constant,1)

                else
                    call wall_distance%read_fields(aux_file)
                end if

                ! Run ChiDG simulation
                call wall_distance%reporter('before')
                call wall_distance%run(write_initial=.false., write_final=.false., write_tecio=.false.)
                call wall_distance%reporter('after')

                ! Write wall distance to auxiliary field
                call wall_distance%write_fields(aux_file)

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

            ! Update p-Poisson fidelity
            p = p_max
            call set_p_poisson_parameter(p)
            call noptions%set('tol',1.e-3_rk)   ! Set nonlinear solver options
            call noptions%set('rtol',1.e-16_rk) ! Set nonlinear solver options
            call noptions%set('ptc',.false.)     ! Set nonlinear solver options
            call noptions%set('smooth',.false.) ! Set nonlinear solver options
            call loptions%set('tol',1.e-7_rk)   ! Set linear solver options
            call loptions%set('rtol',1.e-16_rk) ! Set linear solver options
            call wall_distance%set('Nonlinear Solver', algorithm='Newton', options=noptions)
            call wall_distance%set('Linear Solver'   , algorithm='fgmres', options=loptions)

            order = chidg%nterms_s_1d
            do iorder = 3,order
                call write_line('Wall Distance Driver : Loop 2 : order = ', iorder)

                ! (Re)Initialize domain storage, communication, matrix/vector storage
                call wall_distance%set('Solution Order', integer_input=iorder)
                call wall_distance%init('all')

                ! Read solution if it exists.
                call wall_distance%read_fields(aux_file)

                ! Run ChiDG simulation
                call wall_distance%reporter('before')
                call wall_distance%run(write_initial=.false., write_final=.false., write_tecio=.false.)
                call wall_distance%reporter('after')

                ! Write wall distance to auxiliary field
                call wall_distance%write_fields(aux_file)

            end do

        end if ! have_wd_field .and. wd_file_exists


        ! Read scalar field 'u' to auxiliary field 'Wall Distance : p-Poisson'
        call write_line('Reading Wall Distance field from file to Auxiliary field ChiDG Vector:', io_proc=GLOBAL_MASTER)
        call chidg%read_auxiliary_field(aux_file,field='u',store_as='Wall Distance : p-Poisson')


    end subroutine wall_distance_primal_driver
    !**************************************************************************************



















end submodule type_chidg_driver
