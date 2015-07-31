!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie



program driver
#include <messenger.h>
    use mod_kinds,              only: ik,rk
    use mod_constants,          only: ZERO, ONE, FIVE
    use mod_tecio,              only: write_tecio_variables
    use mod_hdfio,              only: read_grid_hdf
    use type_chidg,             only: chidg_t
    use type_point,             only: point_t
    use type_domain,            only: domain_t
    use type_element,           only: element_t
    use atype_function,         only: function_t
    use mod_grid_operators,     only: initialize_variable
    use mod_function,           only: assign_function

    !-----------------------------------------------------------
    ! Variable declarations
    !-----------------------------------------------------------
    implicit none


    type(chidg_t)                   :: chidg
    type(domain_t),     allocatable :: dom(:)
    class(function_t),  allocatable :: fcn


    call chidg%init()



!    call warn(1,"hi")



!    call signal(WARN,'hi')

    call AllocationError



    call read_grid_hdf('D1_E27_M1.h5',dom)
    call dom(1)%init_sol('scalar',27)



    ! Initialize function
    call assign_function(fcn,'gaussian')

    ! Initialize solution
    call initialize_variable(dom(1),1,fcn)

    call write_tecio_variables(dom(1),'test.plt',1)



end program driver
