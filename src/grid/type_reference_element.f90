module type_reference_element
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/26/2017
    !!
    !-------------------------------------------------------
    type, public :: reference_element_t



    contains

        procedure   :: init

        procedure   :: nodes
        procedure   :: weights

        procedure   :: interpolator


    end type reference_element_t
    !*******************************************************






contains





    call reference_element%nodes(  'uniform',    level=1, dim=3)
    call reference_element%nodes(  'quadrature', level=2, dim=3)
    call reference_element%weights('quadrature', level=2, dim=3)


    call reference_element%nodes(  'uniform',    level=1, dim=2)
    call reference_element%nodes(  'quadrature', level=2, dim=2)
    call reference_element%weights('quadrature', level=2, dim=2)


    
    call reference_element%interpolator('Value')
    call reference_element%interpolator('Grad1')
    call reference_element%interpolator('Grad2')
    call reference_element%interpolator('Grad3')














end module type_reference_element
