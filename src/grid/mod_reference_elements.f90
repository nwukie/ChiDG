module mod_reference_elements
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: NO_ID, ZERO
    use type_reference_element, only: reference_element_t
    implicit none

    ! We don't want this to be allocatable, target, because everytime
    ! it gets extended all the memory locations change and then any
    ! elements/faces that are pointing to those locations in memory
    ! no longer have valid pointers.
    !type(reference_element_t),  allocatable, target :: ref_elems(:)
    type(reference_element_t),  target  :: ref_elems(100)
    integer(ik)                         :: nref_elems = 0   ! number of ref_elem objects that have been initialized

contains


    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !----------------------------------------------------------------------------
    function get_reference_element(element_type,polynomial,nterms,node_set,level,nterms_rule) result(ref_elem_ID)
        integer(ik),    intent(in)              :: element_type
        character(*),   intent(in), optional    :: polynomial
        integer(ik),    intent(in), optional    :: nterms
        character(*),   intent(in), optional    :: node_set
        integer(ik),    intent(in), optional    :: level
        integer(ik),    intent(in), optional    :: nterms_rule

        character(:),   allocatable :: user_msg
        integer(ik)                 :: iref_elem, ref_elem_ID
        logical                     :: polynomial_matches, nterms_matches, node_set_matches, &
                                       level_matches, element_matches, rule_matches

        !
        ! First, try to find an existing reference object to provide
        !
        ref_elem_ID = NO_ID
        do iref_elem = 1,nref_elems

            !
            ! Try to find fully initialized reference element: 
            ! (reference nodes + interpolation nodes)
            !
            if (present(polynomial) .and. &
                present(nterms)     .and. &
                present(node_set)   .and. &
                present(level)      .and. &
                present(nterms_rule) ) then
                if (ref_elems(iref_elem)%interpolation_initialized) then
                    element_matches    = (ref_elems(iref_elem)%element_type == element_type    )
                    polynomial_matches = (ref_elems(iref_elem)%polynomial   == trim(polynomial))
                    nterms_matches     = (ref_elems(iref_elem)%nterms_i()   == nterms          )
                    node_set_matches   = (ref_elems(iref_elem)%node_set     == trim(node_set)  )
                    level_matches      = (ref_elems(iref_elem)%level        == level           )
                    rule_matches       = (ref_elems(iref_elem)%nterms_rule  == nterms_rule     )

                    if (element_matches     .and. &
                        polynomial_matches  .and. &
                        nterms_matches      .and. &
                        node_set_matches    .and. &
                        level_matches       .and. &
                        rule_matches) then
                        ref_elem_ID = iref_elem
                        exit
                    end if
                end if !interpolation_initialized

            !
            ! Try to find reference element only: (reference nodes)
            !
            else if ( (.not. present(polynomial)) .and. &
                      (.not. present(nterms))     .and. &
                      (.not. present(node_set))   .and. &
                      (.not. present(level))      .and. &
                      (.not. present(nterms_rule)) ) then

                element_matches = ref_elems(iref_elem)%element_type == element_type
                if ( element_matches ) then
                    ref_elem_ID = iref_elem
                    exit
                end if

            else
                user_msg = "get_reference_element: Invalid combination of optional parameters &
                            was passed to try and retrieve a reference element object."
                call chidg_signal(FATAL,user_msg)
            end if

        end do !iref_elem



        !
        ! If existing reference element was not found that matched the requirements, 
        ! create a new one.
        !
        if (ref_elem_ID == NO_ID) ref_elem_ID = new_reference_element(element_type, polynomial, nterms, node_set, level, nterms_rule)


    end function get_reference_element
    !****************************************************************************






    !>  Extend the 'ref_elems' module variable allocation and return the ID of 
    !!  the new object.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2017
    !!
    !----------------------------------------------------------------------------
    function new_reference_element(element_type,polynomial,nterms,node_set,level, nterms_rule) result(ref_elem_ID)
        integer(ik),    intent(in)              :: element_type
        character(*),   intent(in), optional    :: polynomial
        integer(ik),    intent(in), optional    :: nterms
        character(*),   intent(in), optional    :: node_set
        integer(ik),    intent(in), optional    :: level
        integer(ik),    intent(in), optional    :: nterms_rule


        character(:),               allocatable :: user_msg
        integer(ik)                             :: ref_elem_ID, ierr


        !
        ! Initialize new object(last 
        !
        ref_elem_ID = nref_elems + 1
        nref_elems  = nref_elems + 1


        !
        ! All new elements are initialized with their reference nodes and projector
        !
        call ref_elems(ref_elem_ID)%init_element(element_type)



        !
        ! If specified, initialize the new object with interpolators: 
        ! (reference nodes + interpolation nodes)
        !
        if (present(polynomial) .and. &
            present(nterms)     .and. &
            present(node_set)   .and. &
            present(level)      .and. &
            present(nterms_rule)) then
            call ref_elems(ref_elem_ID)%init_interpolator(polynomial,nterms,node_set,level,nterms_rule)
        end if !iref_elem


    end function new_reference_element
    !****************************************************************************





end module mod_reference_elements
