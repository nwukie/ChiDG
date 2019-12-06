module mod_tutorial_smoothbump_new
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_string,             only: string_t
    use mod_file_utilities,     only: delete_file
    use type_tutorial,          only: tutorial_t
    use type_tutorial_step,     only: tutorial_step_t
    use type_bc_state_group,    only: bc_state_group_t
    use mod_test_utilities,     only: create_mesh_file
    use type_bc_state,          only: bc_state_t
    use mod_bc,                 only: create_bc
    use mod_io
    implicit none


    !>
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    type, extends(tutorial_t), public   :: smooth_bump_tutorial_t

    contains

        procedure   :: print_configuration
        procedure   :: initialize
        procedure   :: clean_up

    end type smooth_bump_tutorial_t
    !*****************************************************************************




    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step1
    contains
        procedure   :: execute => execute_step1
    end type bump_step1
    !*******************************************************


    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step2
    contains
        procedure   :: execute => execute_step2
    end type bump_step2
    !*******************************************************


    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step3
    contains
        procedure   :: execute => execute_step3
    end type bump_step3
    !*******************************************************


    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step4
    contains
        procedure   :: execute => execute_step4
    end type bump_step4
    !*******************************************************

    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step5
    contains
        procedure   :: execute => execute_step5
    end type bump_step5
    !*******************************************************

    !-------------------------------------------------------
    type, extends(tutorial_step_t), public  :: bump_step6
    contains
        procedure   :: execute => execute_step6
    end type bump_step6
    !*******************************************************


contains




    !>  Initialize the smooth bump tutorial steps.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine initialize(self)
        class(smooth_bump_tutorial_t),  intent(inout)   :: self
        
        type(bump_step1)   :: step1
        type(bump_step2)   :: step2
        type(bump_step3)   :: step3
        type(bump_step4)   :: step4
        type(bump_step5)   :: step5
        type(bump_step6)   :: step6

        !
        ! Set tutorial text data
        !
        call self%set_title('Flow over smooth bump')
        call self%set_tags('Euler')
        call self%set_files('smooth_bump.x, smooth_bump_set.h5, chidg.nml')


        !
        ! Set step titles
        !
        call step1%set_title('Overview')
        call step2%set_title('Generate grid + boundary conditions')
        call step3%set_title('Run calculation')
        call step4%set_title('Post-processing')
        call step5%set_title('Other resources')
        call step6%set_title('Tutorial complete')
        

        !
        ! Add steps
        !
        call self%add_step(step1)
        call self%add_step(step2)
        call self%add_step(step3)
        call self%add_step(step4)
        call self%add_step(step5)
        call self%add_step(step6)

        !stop
    end subroutine initialize
    !*********************************************************************************




    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine print_configuration(self)
        class(smooth_bump_tutorial_t),  intent(inout)   :: self

        call write_line('                                                                           ',ltrim=.false.)
        call write_line('Configuration:                                                             ',ltrim=.false., bold=.true.)
        call write_line('                                   Wall                                    ',ltrim=.false.)
        call write_line('                       ___________________________                         ',ltrim=.false.)
        call write_line('                      |                           |                        ',ltrim=.false.)
        call write_line('   Total Pressure     |         Flow --->         |  Constant Pressure     ',ltrim=.false.)
        call write_line('   Total Temperature  |           _____           |                        ',ltrim=.false.)
        call write_line('                      |__________/     \__________|                        ',ltrim=.false.)
        call write_line('                                   Wall                                    ',ltrim=.false.)
        call write_line('                                                                           ',ltrim=.false.)


    end subroutine print_configuration
    !*********************************************************************************



    !>
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine clean_up(self)
        class(smooth_bump_tutorial_t),  intent(inout)   :: self

        call delete_file('smooth_bump.x')
        call delete_file('smooth_bump_setup.h5')


    end subroutine clean_up
    !*********************************************************************************






    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step1(self)
        class(bump_step1), intent(inout)   :: self

        character(:),   allocatable :: overview


        overview = "This tutorial investigates the problem of inviscid flow over a smooth &
                    bump in a channel. Since the flow is inviscid and subsonic, there should &
                    be no entropy generated in the problem. However, the solution of the &
                    discretized problem will have some entropy generated. We can use the &
                    ammount of entropy generated in the numerical solution as an error metric &
                    for determining the order of accuracy for the numerical scheme."

        call write_line(overview, width=100, color='blue')


    end subroutine execute_step1
    !*********************************************************************************


    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step2(self)
        class(bump_step2), intent(inout)   :: self

        type(string_t)                  :: group_names(1,6)
        type(bc_state_group_t)          :: bc_state_groups(3)
        class(bc_state_t),  allocatable :: bc_state
        character(:),       allocatable :: step2_1, step2_2

        step2_1 = "A multi-block, double-precision, unformatted Plot3D grid for this problem &
                 has been created and is stored in the file, 'smooth_bump.x'. The Plot3D &
                 grid was generated with extra resolution so that the points can be agglomerated &
                 by ChiDG to create quartic, curved elements. The Plot3D grid generally needs &
                 to be generated by the user. The Plot3D grid must be converted to the ChiDG format &
                 using the 'convert' utility as 'chidg convert smooth_bump.x' and then the boundary &
                 conditions can be added and set in the 'edit' utility as 'chidg edit smooth_bump.h5'. & 
                 This has been done for you here and the ChiDG grid is stored in 'smooth_bump_setup.h5' &
                 with the boundary condition functions and patch groups already set up. You may get a &
                 feel for the conversion and edit process by running the 'convert' action on the Plot3D &
                 grid(chidg convert smooth_bump.x) and then 'edit' the generated file to set &
                 boundary conditions(chidg edit smooth_bump.h5). For now we will use the file &
                 'smooth_bump_setup.h5' that has already been converted and set up. You can use this &
                 file as a reference and view the setup by running 'chidg edit smooth_bump_setup.h5'."

        step2_2 = "You may investigate the contents of the generated HDF5 file by running the &
                  action 'chidg edit smooth_bump_setup.h5' in another terminal. Just make sure to &
                  exit the edit session before proceeding with this tutorial, since only one process &
                  may access the file at any given time."



        call write_line(step2_1, width=100, color='blue')
        call write_line(' ', ltrim=.false.)


        !
        ! Setup boundary conditions
        !
        bc_state_groups(1)%name = 'Inlet'
        bc_state_groups(2)%name = 'Outlet'
        bc_state_groups(3)%name = 'Walls'


        call create_bc('Inlet - Total', bc_state)
        call bc_state%set_fcn_option('Total Pressure',   'val', 110000._rk)
        call bc_state%set_fcn_option('Total Temperature','val', 300._rk   )
        call bc_state_groups(1)%add_bc_state(bc_state)

        call create_bc('Outlet - Constant Pressure', bc_state)
        call bc_state%set_fcn_option('Static Pressure', 'val', 93000._rk )
        call bc_state_groups(2)%add_bc_state(bc_state)

        call create_bc('Wall', bc_state)
        call bc_state_groups(3)%add_bc_state(bc_state)


        ! Define patch group names
        group_names(1,:) = [string_t('Inlet') , &
                            string_t('Outlet'), &
                            string_t('Walls') , &
                            string_t('Walls') , &
                            string_t('Walls') , &
                            string_t('Walls') ]

        !
        ! Create grid file
        !
        call write_line("Creating grid file:", width=90, color='blue', bold=.true.)
        call create_mesh_file(selector        = "Smooth Bump",          &
                              filename        = 'smooth_bump_setup.h5', &
                              equation_sets   = [string_t("Euler")],    &
                              group_names     = group_names,            &
                              bc_state_groups = bc_state_groups,        &
                              nelem_xi        = 15,                     &
                              nelem_eta       = 8,                      &
                              nelem_zeta      = 1,                      &
                              save_intermediate_files = .true.)
        call write_line(" ", width=90, color='blue', ltrim=.false.)



        !
        ! Set namelist values and write to chidg.nml                
        !
        gridfile            = 'smooth_bump_setup.h5'
        solutionfile_out    = 'smooth_bump_setup.h5'
        solution_order      = 1
        gq_rule             = 3
        time_integrator     = 'steady'
        nonlinear_solver    = 'newton'
        cfl0                = 1.0
        ntol                = 1.e-5
        linear_solver       = 'fgmres'
        ltol                = 1.e-9
        preconditioner      = 'RASILU0'
        initial_fields(1:5) = [1.19, 100.0, 0.0, 0.0, 250000.0]


        call write_line("Creating namelist input file: 'chidg.nml'", width=90, color='blue', bold=.true.)
        call write_line(" ", width=90, color='blue', bold=.true., ltrim=.false.)
        call write_namelist()


        call write_line("Equivalent steps:", color='blue', bold=.true.)
        call write_line("   - Create Plot3D grid: smooth_bump.x.",                                        color='blue', ltrim=.false.)
        call write_line("   - Convert Plot3D grid to HDF5: 'chidg convert grid.x'",                       color='blue', ltrim=.false.)
        call write_line("   - Edit boundary state functions/boundary patch groups: 'chidg edit grid.h5'", color='blue', ltrim=.false.)
        call write_line(" ", width=90, color='blue', ltrim=.false.)

        call write_line("Note:", color='blue', bold=.true.)
        call write_line(step2_2, width=100, color='blue')


    end subroutine execute_step2
    !*********************************************************************************





    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step3(self)
        class(bump_step3), intent(inout)   :: self

        character(:),   allocatable :: step3_1, step3_2, step3_3, step3_4

        step3_1 = "Run the ChiDG calculation. In a separate terminal, navigate to the directory that the &
                 tutorial is running in. Run the ChiDG executable. If the executable is in your path: &
                 ( user@:~$ chidg )."
        step3_2 = "1) ChiDG searches for 'chidg.nml' and obtains the name of the grid files. Also initializes solver settings."
        step3_3 = "2) From the file specified in 'chidg.nml', ChiDG reads the grid and initializes data."
        step3_4 = "3) ChiDG runs the calculation."


        call write_line(step3_1, width=120, color='blue')
        call write_line(' ', ltrim=.false.)
        call write_line('During execution:', width=120, color='blue', bold=.true.)
        call write_line(step3_2, width=120, color='blue')
        call write_line(step3_3, width=120, color='blue')
        call write_line(step3_4, width=120, color='blue')


    end subroutine execute_step3
    !*********************************************************************************



    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step4(self)
        class(bump_step4), intent(inout)   :: self

        character(:),   allocatable :: step4_1, step4_2


        step4_1 = "Once the calculation is finished, in the same terminal and in the same tutorial directory, &
                 we want to post-process the solution for visualization. ChiDG can generate visualization &
                 files for either Tecplot360(.szplt) or Paraview(.pvd, VTK). These are generated using &
                 the '2tec' and '2vtk' actions respectively as: 'chidg 2tec smooth_bump_setup.h5' or &
                 'chidg 2vtk smooth_bump_setup.h5'."
        step4_2 = "The ChiDG post-processing action('chidg 2tec' or 'chidg 2vtk') expects a ChiDG file containing &
                   a solution."


        call write_line(step4_1, width=120, color='blue')
        call write_line(' ', ltrim=.false.)
        call write_line(step4_2, width=120, color='blue')


    end subroutine execute_step4
    !*********************************************************************************




    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step5(self)
        class(bump_step5), intent(inout)   :: self

        character(:),   allocatable :: step5_1, step5_2, step5_3, step5_4

        step5_1 = "International Workshop on High-Order CFD Methods: https://how4.cenaero.be/content/bi2-inviscid-flow-over-bump"
        step5_2 = "www.tecplot.com"
        step5_3 = "www.paraview.org"
        step5_4 = "Reach out with questions to Nathan Wukie(nwukie@gmail.com)."


        call write_line(' ', ltrim=.false.)
        call write_line(' ', ltrim=.false.)
        call write_line('General problem specification: ',color='blue', bold=.true.)
        call write_line(step5_1, width=150, color='blue')
        call write_line(' ', ltrim=.false.)
        call write_line(' ', ltrim=.false.)
        call write_line('Visualization applications: ', color='blue', bold=.true.)
        call write_line(step5_2, width=150, color='blue')
        call write_line(step5_3, width=150, color='blue')
        call write_line(' ', ltrim=.false.)
        call write_line(' ', ltrim=.false.)
        call write_line('Developer contact: ', color='blue', bold=.true.)
        call write_line(step5_4, width=150, color='blue')

    end subroutine execute_step5
    !*********************************************************************************





    !>
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine execute_step6(self)
        class(bump_step6), intent(inout)   :: self

        character(:),   allocatable :: step6



    end subroutine execute_step6
    !*********************************************************************************






end module mod_tutorial_smoothbump_new
