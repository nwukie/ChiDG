!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie



program driver
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


    ! Initialize ChiDG environment
    call chidg%init()

    ! Allocate solver
    call create_solver('fe',solver)


    ! Initialize grid and solution
    call read_grid_hdf('9x9x9.h5',domains)
    call domains(1)%init_sol('scalar',64)


    call create_function(fcn,'gaussian')
    call initialize_variable(domains(1),1,fcn)


    call write_tecio_variables(domains(1),'test.plt',1)



end program driver
