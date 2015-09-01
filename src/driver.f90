!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie



program driver
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: XI_MIN, ETA_MIN, ZETA_MIN
    use type_chidg,             only: chidg_t
    use type_domain,            only: domain_t
    use atype_solver,           only: solver_t
    use mod_solver,             only: create_solver
    use atype_matrixsolver,     only: matrixsolver_t
    use mod_matrixsolver,       only: create_matrixsolver
    use mod_hdfio,              only: read_grid_hdf
    use mod_grid_operators,     only: initialize_variable
    use atype_function,         only: function_t
    use mod_function,           only: create_function
    use mod_tecio,              only: write_tecio_variables
    use mod_io
    
    !-----------------------------------------------------------
    ! Variable declarations
    !-----------------------------------------------------------
    implicit none
    type(chidg_t)                       :: chidg
    type(domain_t),         allocatable :: domains(:)
    class(matrixsolver_t),  allocatable :: matrixsolver
    class(solver_t),        allocatable :: solver
    class(function_t),      allocatable :: fcn
   

    !
    ! Initialize ChiDG environment
    !
    call chidg%init('full')


    !
    ! Allocate solution-advancement routine
    !
    call create_solver(temporal_scheme,solver)


    !
    ! Initialize grid and numerics
    !
    call read_grid_hdf(gridfile,domains)




    call domains(1)%init_bc('periodic',XI_MIN)
    call domains(1)%init_bc('periodic',ETA_MIN)
    call domains(1)%init_bc('periodic',ZETA_MIN)

    call domains(1)%init_sol(eqnset,nterms_s)


    !
    ! Initialize solution
    !
    call create_function(fcn,'gaussian')
    call initialize_variable(domains(1),1,fcn)



    !
    ! Write initial solution
    !
    !call write_tecio_variables(domains(1),'0.plt',1)



    !
    ! Initialize solver
    !
    call solver%init(domains(1))



    !
    ! Create matrix solver
    !
    call create_matrixsolver(msolver,matrixsolver)



    !
    ! Execute solver routine
    !
    call solver%solve(domains(1),matrixsolver)



end program driver
