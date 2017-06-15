module type_support
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    implicit none






    !>  
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/14/2017
    !!
    !-------------------------------------------------------------------
    type, public support_t

    contains

        procedure   :: 

    end type support_t
    !*******************************************************************



contains





    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/14/2017
    !!
    !-------------------------------------------------------------------
    function produce(self,geometry,set,orientation,level) result(nodes)
        class(support_t),   intent(in)  :: self
        character(*),       intent(in)  :: geometry
        character(*),       intent(in)  :: set
        integer(ik),        intent(in)  :: orientation
        integer(ik),        intent(in)  :: level

        real(rk)    :: nodes(:,:)



    end function produce
    !*******************************************************************




    call geometry%set_geometry(nodes, connectivity)
    call geometry%set_basis(basis, order)
    call geometry%set_interpolation(rule, level)
    call geometry%construct()


    points = cloud%produce('Equally Spaced', 
    


    call geometry%construct(nodes,connectivity)
    call geometry%support(rule='Quadrature', level=1)
    call geometry%enrich(basis='Legendre',   order=2)


    nodes = geometry%node_set('Hexahedral',    'Quadrature', face=1,      level=1)
    nodes = geometry%node_set('Quadrilateral', 'Quadrature', component=1, level=1)
    nodes = geometry%node_set('Line',          'Quadrature',              level=1)


    call interpolator%construct('Tensor-Product',nodes)































end module type_support
