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
    use type_chidg,             only: chidg_t
    use type_domain,            only: domain_t
    use atype_solver,           only: solver_t
    use mod_solver,             only: create_solver
    use mod_hdfio,              only: read_grid_hdf
    use mod_grid_operators,     only: initialize_variable
    use atype_function,         only: function_t
    use mod_function,           only: create_function
    use mod_tecio,              only: write_tecio_variables

    !-----------------------------------------------------------
    ! Variable declarations
    !-----------------------------------------------------------
    implicit none
    type(chidg_t)                   :: chidg
    type(domain_t),     allocatable :: domains(:)
    class(solver_t),    allocatable :: solver
    class(function_t),  allocatable :: fcn
    integer(ik)     :: ielem


    ! Initialize ChiDG environment
    print*, 'calling chidg%init'
    call chidg%init()

    ! Allocate solver
    print*, 'creating solver'
    call create_solver('fe',solver)


    ! Initialize grid and solution
    print*, 'reading grid'
    call read_grid_hdf('4x4x4.h5',domains)

    print*, 'initializing domain numerics'
    call domains(1)%init_sol('scalar',8)


    print*, 'create_function'
    call create_function(fcn,'gaussian')


    print*, 'initialize_variable'
    call initialize_variable(domains(1),1,fcn)


    call write_tecio_variables(domains(1),'0.plt',1)










    call solver%init(domains(1))

    print*, 'solve'
    call solver%solve(domains(1))



end program driver
