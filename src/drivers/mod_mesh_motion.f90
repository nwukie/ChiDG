module mod_mesh_motion
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO
    use mod_function,           only: create_function
    use type_function,          only: function_t
    use type_chidg,             only: chidg_t
    use type_chidg_vector,             only: chidg_vector_t
    use type_bc_state,          only: bc_state_t
    use type_dict,              only: dict_t
    use mod_chidg_post,         only: chidg_post,chidg_post_vtk
    use mod_bc,                 only: create_bc
    use mod_io,                 only: gridfile
    use mod_hdf_utilities,      only: get_properties_hdf, check_file_exists_hdf
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_file_properties,   only: file_properties_t
    implicit none





contains


    !>  Solve for mesh motion using a PDE-based approach.
    !!
    !!  Solve a diffusion equation for the approximate mesh displacements. 
    !!
    !!  @date  4/26/2017 
    !!  @author Eric Wolf 
    !!
    !!  @param[in]  gridfile        ChiDG-file containing a grid to be read and computed on.
    !!  @param[in]  solutionfile    ChiDG-file the field will be written to.
    !!  @param[in]  order           Polynomial order, the field will be computed with.
    !!
    !-------------------------------------------------------------------------------------
    subroutine mesh_motion_driver(chidg,fileout)
        type(chidg_t),  intent(inout)           :: chidg
        character(*),   intent(in), optional    :: fileout

        character(:), allocatable   :: user_msg
        integer(ik)                 :: order

        type(chidg_t)                   :: mesh_motion,chidg_temp
        type(dict_t)                    :: noptions, loptions
        class(bc_state_t),  allocatable :: dirichlet_zero, neumann_zero, dirichlet_dspl
        class(function_t),  allocatable :: mmd_initial, mmd_analytical 
        integer(ik)                     :: iorder, p, wd_nterms_s, ierr, iproc
        integer(ik)                     :: aux_field_index1, aux_field_index2, aux_field_index3
        type(file_properties_t)         :: wd_props
        logical                         :: wd_file_exists, have_wd_field
        real(rk)                        :: error_val



        !
        ! chidg%init('mpi') should have already been called by another ChiDG instance.
        ! chidg%init('io')  should have already been called by another ChiDG instance.
        !
        ! Make sure this ChiDG environment is initialized.
        !
        call mesh_motion%start_up('core')




        !
        ! Initialize options dictionaries
        !
        call noptions%set('tol',1.e-4_rk)   ! Set nonlinear solver options
        call noptions%set('cfl0',0.1_rk)
        call noptions%set('nsteps',50)
        call loptions%set('tol',1.e-7_rk)  ! Set linear solver options



        !
        ! Set ChiDG components
        !
        call mesh_motion%set('Time Integrator' , algorithm='Steady'                        )
        call mesh_motion%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
        call mesh_motion%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)
        call mesh_motion%set('Preconditioner'  , algorithm='RASILU0'                       )

        order = chidg%nterms_s_1d
        call mesh_motion%set('Solution Order', integer_input=order)

        !mesh_motion%nonlinear_solver%search = .false.
        mesh_motion%nonlinear_solver%search = .true.


        !
        ! Set up boudary condition states to impose.
        !
        ! u = 0 on walls
        !
        ! dudn = 0 else where
        !
        call create_bc('Scalar Value',      dirichlet_zero)
        call create_bc('Scalar Value',      dirichlet_dspl)
        call create_bc('Scalar Derivative', neumann_zero  )

        call dirichlet_zero%set_fcn_option( 'Value'     , 'val', ZERO)
        call dirichlet_dspl%set_fcn_option( 'Value'     , 'val', 0.1_rk)
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
        call mesh_motion%read_grid(gridfile, equation_set='Mesh Motion : Diffusion')
        call mesh_motion%read_boundaryconditions(gridfile, bc_wall     = neumann_zero, &
                                                             bc_inlet    = dirichlet_dspl,   &
                                                             bc_outlet   = dirichlet_zero,   &
                                                             bc_symmetry = neumann_zero,   &
                                                             bc_farfield = neumann_zero)



        !
        ! Initialize mesh_motion with chidg order in case we are going to read in a solution.
        !
        call mesh_motion%init('all')




        !
        ! Check if we already have a wall distance solution
        !
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then

                wd_file_exists = check_file_exists_hdf(fileout)
                if (wd_file_exists) then
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

            call mesh_motion%read_solution(fileout)
            mesh_motion%data%sdata%q = mesh_motion%data%sdata%q_in

        !
        ! If we don't have an accurate wall distance field in file, solve for a new one.
        !
        else

            !
            ! Store grid to file
            !
            call mesh_motion%write_grid(fileout)

            
            iorder = 2
            call noptions%set('tol',1.e-4_rk)   ! Set nonlinear solver options
            call loptions%set('tol',1.e-5_rk)   ! Set linear solver options
            call mesh_motion%set('Nonlinear Solver', algorithm='Newton', options=noptions)
            call mesh_motion%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)
            call write_line('Mesh Motion Driver')

            !
            ! (Re)Initialize domain storage, communication, matrix/vector storage
            !
            call mesh_motion%set('Solution Order', integer_input=iorder)
            call mesh_motion%init('all')


            !
            ! Read solution if it exists.
            !
            call create_function(mmd_initial,'constant')
            call mmd_initial%set_option('val',0.0_rk)
            call mesh_motion%data%sdata%q_in%project(mesh_motion%data%mesh,mmd_initial,1)

            
            !
            ! Run ChiDG simulation
            !
            call mesh_motion%report('before')
            call mesh_motion%run(write_initial=.false., write_final=.false.)
            
!!
!            ! Read solution if it exists.
!            !
!            call create_function(mmd_analytical,'mmd_cdiff')
!            call mmd_analytical%set_option('dspl',0.1_rk) !Prescribed displacement of the moving wall (x1=0 face)
!            
!!            call create_function(mmd_analytical,'mmd_ldiff')
!!            call mmd_analytical%set_option('dspl',0.1_rk) !Prescribed displacement of the moving wall (x1=0 face)
!!            call mmd_analytical%set_option('frac',0.1_rk) !0<frac<1, rate of decrease of diffusivity away from moving wall (along x1)
!!            call mmd_analytical%set_option('xstart', ZERO)
!!            call mmd_analytical%set_option('xlength', ONE)
!
!
!            q_ref = mesh_motion%data%sdata%q
!            call q_ref%project(mesh_motion%data%mesh,mmd_analytical,1)
!
!            !Implement this!
!            ! Computes the L2 state error against a reference ChiDG vector.
!            error_val = mesh_motion%compute_l2_state_error(q_ref)

            call mesh_motion%report('after')


            !
            ! Write wall distance to auxiliary field
            !
            call mesh_motion%write_solution(fileout)



        end if ! have_wd_field .and. wd_file_exists

!        call write_line('Storing Mesh Motion field to Auxiliary field ChiDG Vector:', io_proc=GLOBAL_MASTER)
!
!        !Initialize chidg instance as temporary storage to copy grid displacement components into the 
!        call chidg_temp%start_up('core')
!        call chidg_temp%set('Time Integrator' , algorithm='Steady'                        )
!        call chidg_temp%set('Nonlinear Solver', algorithm='Quasi-Newton', options=noptions)
!        call chidg_temp%set('Linear Solver'   , algorithm='fgmres_cgs',   options=loptions)
!        call chidg_temp%set('Preconditioner'  , algorithm='RASILU0'                       )
!
!        call chidg_temp%set('Solution Order', integer_input=order)
!        call chidg_temp%init('all')
!
!        !chidg_temp%nonlinear_solver%search = .false.
!        chidg_temp%nonlinear_solver%search = .true.
!
!
!        !
!        ! Read grid, boundary conditions.
!        !
!        
!        call chidg_temp%read_grid(gridfile, equation_set='Scalar Diffusion')
!        
!        !
!        ! Try to find 'Mesh Motion' auxiliary field storage.
!        !
!        aux_field_index1 = chidg%data%sdata%get_auxiliary_field_index('Mesh Motion Grid Displacement 1')
!        aux_field_index2 = chidg%data%sdata%get_auxiliary_field_index('Mesh Motion Grid Displacement 2')
!        aux_field_index3 = chidg%data%sdata%get_auxiliary_field_index('Mesh Motion Grid Displacement 3')
!
!
!
!        !
!        ! Generalize the following code to copy the appropriate sections of q in the grid displacement aux fields
!        !   Or: Can an aux field consist of multiple equations/components?
!
!        !
!        ! If no 'Mesh Motion' auxiliary field storage was not found, create one.
!        !
!        !
!
!        ieqn = 1 
!        do idom = 1, ndom
!            do ielem = 1, nelems
!                chidg_temp%data%sdata%q%doms(idom)%vecs(ielem)%setvar( 1,1,& 
!                    mesh_motion%data%sdata%q%doms(idom)%vecs(ielem)%getvar(ieqn,1)) = 
!            end do
!        end do
!
!
!        if (aux_field_index1 == 0) then
!
!            
!            call chidg%data%sdata%add_auxiliary_field('Mesh Motion Grid Displacement 1', chidg_temp%data%sdata%q)
!
!        !
!        ! If 'Mesh Motion' auxiliary field storage was found, copy Wall Distance solution 
!        ! to working ChiDG environment.
!        !
!        else
!            chidg%data%sdata%auxiliary_field(aux_field_index1) = chidg_temp%data%sdata%q
!
!        end if
!
!        ieqn = 2 
!        do idom = 1, ndom
!            do ielem = 1, nelems
!                chidg_temp%data%sdata%q%doms(idom)%vecs(ielem)%setvar( 1,1,& 
!                    mesh_motion%data%sdata%q%doms(idom)%vecs(ielem)%getvar(ieqn,1)) = 
!            end do
!        end do
!
!
!        if (aux_field_index1 == 0) then
!
!            
!            call chidg%data%sdata%add_auxiliary_field('Mesh Motion Grid Displacement 2', chidg_temp%data%sdata%q)
!
!        !
!        ! If 'Mesh Motion' auxiliary field storage was found, copy Wall Distance solution 
!        ! to working ChiDG environment.
!        !
!        else
!            chidg%data%sdata%auxiliary_field(aux_field_index2) = chidg_temp%data%sdata%q
!
!        end if
!
!        ieqn = 3 
!        do idom = 1, ndom
!            do ielem = 1, nelems
!                chidg_temp%data%sdata%q%doms(idom)%vecs(ielem)%setvar( 1,1,& 
!                    mesh_motion%data%sdata%q%doms(idom)%vecs(ielem)%getvar(ieqn,1)) = 
!            end do
!        end do
!
!
!        if (aux_field_index1 == 0) then
!
!            
!            call chidg%data%sdata%add_auxiliary_field('Mesh Motion Grid Displacement 3', chidg_temp%data%sdata%q)
!
!        !
!        ! If 'Mesh Motion' auxiliary field storage was found, copy Wall Distance solution 
!        ! to working ChiDG environment.
!        !
!        else
!            chidg%data%sdata%auxiliary_field(aux_field_index3) = chidg_temp%data%sdata%q
!
!        end if



    end subroutine mesh_motion_driver
    !**************************************************************************************




end module mod_mesh_motion
