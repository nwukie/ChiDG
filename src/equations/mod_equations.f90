module mod_equations
    use mod_kinds,                      only: rk,ik
    use atype_equationset,              only: equationset_t
    use type_testeq,                    only: testeq_t
    implicit none


contains


    subroutine AssignEquationSet(eqn_string,eqnset)
        character(*),                      intent(in)      :: eqn_string
        class(equationset_t), allocatable, intent(inout)   :: eqnset

        type(testeq_t)  :: testeq




        if (eqn_string == 'testeq') then
            call testeq%init()
            allocate(eqnset, source=testeq)
        end if

    end subroutine



end module mod_equations
