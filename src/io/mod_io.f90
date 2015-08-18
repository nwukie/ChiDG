module mod_io
    !!          Module Includes             !!
    use mod_kinds,     only: rk,ik
    use mod_constants, only: MAXBLOCKS

    !!          Variable Declarations       !!
    implicit none

    integer(kind=ik),   save    :: sol_poly_order = 2
    integer(kind=ik),   save    :: mesh_poly_order = 1

    integer(kind=ik),   save    :: gq_rule = 2          !> 1: Collocation, 2: Over-integration
    integer(kind=ik),   save    :: nterms_sol1d = 2
    integer(kind=ik),   save    :: nterms_sol2d = 2
    integer(kind=ik),   save    :: nterms_sol3d = 2
    integer(kind=ik),   save    :: nterms_mesh1d = 1
    integer(kind=ik),   save    :: nterms_mesh2d = 4
    integer(kind=ik),   save    :: nterms_mesh3d = 8

    integer(kind=ik),   save    :: ntime_instances = 1
    character(len=20),  save    :: solver
    real(kind=rk),      save    :: dt = 0.001_rk
    integer(kind=ik),   save    :: nsteps = 100
    integer(kind=ik),   save    :: nwrite = 100

    character(len=20),  save    :: gridfile,gridtype
    integer(kind=ik),   save    :: bc_ximin(MAXBLOCKS),   bc_ximax(MAXBLOCKS)
    integer(kind=ik),   save    :: bc_etamin(MAXBLOCKS),  bc_etamax(MAXBLOCKS)
    integer(kind=ik),   save    :: bc_zetamin(MAXBLOCKS), bc_zetamax(MAXBLOCKS)
    real(kind=rk),      save    :: bcpar1(6), bcpar2(6), bcpar3(6), bcpar4(6)

    integer(kind=ik),   save    :: output_res = 5


    contains
    !--------------------------------------------------

    subroutine read_input()
    !   read case variables from namelist input file
    !   Edit list:
    !       Nathan Wukie - 2/9/2015
        use mod_ordering,   only: MAX_POLY_ORDER
        implicit none

        namelist /global/   sol_poly_order,     &
                            mesh_poly_order,    &
                            gq_rule,            &
                            ntime_instances

        namelist /time/     solver, dt, nsteps, nwrite


        namelist /grid/     gridfile,   gridtype,   &
                            bc_ximin,   bc_ximax,   &
                            bc_etamin,  bc_etamax,  &
                            bc_zetamin, bc_zetamax, &
                            bcpar1, bcpar2, bcpar3, bcpar4

        namelist /io/       output_res

        open(unit=7,form='formatted',file="case.inp")
        read(7,nml=global)
        read(7,nml=time)
        read(7,nml=grid)
        read(7,nml=io)

        ! solution terms
        if(sol_poly_order > MAX_POLY_ORDER) stop "Error: maximum solution order exceeded"
        nterms_sol1d = sol_poly_order+1
        nterms_sol2d = nterms_sol1d*nterms_sol1d
        nterms_sol3d = nterms_sol1d*nterms_sol1d*nterms_sol1d

        ! mesh terms
        nterms_mesh1d = mesh_poly_order+1
        nterms_mesh2d = nterms_mesh1d*nterms_mesh1d
        nterms_mesh3d = nterms_mesh1d*nterms_mesh1d*nterms_mesh1d

    end subroutine read_input




end module mod_io
