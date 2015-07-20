!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie


program driver
    use mod_kinds,              only: ik,rk
    use mod_constants,          only: ZERO, ONE, FIVE
    use mod_tecio,              only: write_tecio_variables
    use mod_hdfio,              only: read_grid_hdf5
    use type_chidg,             only: chidg_t
    use type_point,             only: point_t
    use type_domain,            only: domain_t
    use type_element,           only: element_t
    use atype_function,         only: function_t
    use mod_grid_operators,     only: initialize_variable
    use mod_function,           only: assign_function
    use messenger

    !-----------------------------------------------------------
    ! Variable declarations
    !-----------------------------------------------------------
    implicit none
    type(chidg_t)                   :: chidg
!    type(domain_t), allocatable     :: dom(:)
!    type(domain_t), target              :: dom(1)
    type(domain_t), target              :: dom
    integer(ik)                     :: ipt, ixi, ieta, izeta
    real(rk)                        :: x(12), y(12), z(12)
    type(point_t)                   :: pts(3,2,2)
    class(function_t), allocatable  :: fcn
    type(element_t), pointer                :: elems(:)
    type(element_t), pointer                :: elems_3d(:,:,:)
    type(element_t), allocatable, target    :: elems_1d(:)


    call chidg%init()

    allocate(elems_1d(2))

    !> Initialize coredg static data
    x = [ZERO, FIVE, 10._rk, ZERO, FIVE, 10._rk, ZERO, FIVE, 10._rk, ZERO, FIVE, 10._rk]
    y = [ZERO, ZERO, ZERO,   ONE,  ONE,  ONE,    ZERO, ZERO, ZERO,   ONE,  ONE,  ONE]
    z = [ZERO, ZERO, ZERO,   ZERO, ZERO, ZERO,   ONE,  ONE,  ONE,    ONE,  ONE,  ONE]

    ipt = 1
    do izeta = 1,2
        do ieta = 1,2
            do ixi = 1,3
                call pts(ixi,ieta,izeta)%set(x(ipt),y(ipt),z(ipt))
                ipt = ipt + 1
            end do
        end do
    end do

!    call read_grid_hdf5('D1_E27_M1.h5',dom)

!    call dom(1)%init_geom(8,pts)
!    call dom(1)%init_sol('scalar',27)

    call dom%init_geom(8,pts)
    call dom%init_sol('scalar',27)

    elems_1d(1:2) = dom%mesh%elems(1:2)

    elems_3d(1:2,1:1,1:1) => elems_1d(1:2)
    elems(1:2)            => dom%mesh%elems

    print*, "driver - elems_3d"
    print*, elems_3d(1,1,1)%ielem
    print*, elems_3d(2,1,1)%ielem

    print*, "driver - elems"
    print*, elems(1)%ielem
    print*, elems(2)%ielem

    print*, "reference - dom%mesh%elems"
    print*, dom%mesh%elems(2)%ielem
    print*, dom%mesh%elems(1)%ielem


!    ! Initialize function
!    call assign_function(fcn,'gaussian')
!
!    ! Initialize solution
!    call initialize_variable(dom(1),1,fcn)
!
!    call write_tecio_variables(dom(1),'test.plt',1)



end program driver
