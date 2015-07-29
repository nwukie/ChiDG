module atype_equationset
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: XI_MIN,XI_MAX,ETA_MIN,ETA_MAX, &
                              ZETA_MIN,ZETA_MAX,NFACES

    use type_equation,  only: equation_t
    use type_mesh,      only: mesh_t
    use type_expansion, only: expansion_t


    implicit none
    private

    type, public, abstract :: equationset_t
        character(100)              :: name
        integer(ik)                 :: neqns

        ! List equation variables with indices
        type(equation_t), allocatable  :: eqns(:)


    contains
        ! Must define these procedures in the extended type
        procedure(self_interface),     deferred  :: init
        procedure(boundary_interface), deferred  :: compute_boundary_average_flux
        procedure(boundary_interface), deferred  :: compute_boundary_upwind_flux
        procedure(volume_interface),   deferred  :: compute_volume_flux
        procedure(volume_interface),   deferred  :: compute_volume_source


        procedure   :: get_var

    end type equationset_t

    abstract interface
        subroutine self_interface(self)
            import equationset_t
            class(equationset_t), intent(inout) :: self
        end subroutine
    end interface


    abstract interface
        subroutine boundary_interface(self,mesh,q,rhs,ielem,iface,iblk)
            use mod_kinds,  only: ik
            import equationset_t
            import mesh_t
            import expansion_t

            class(equationset_t),   intent(in)          :: self
            class(mesh_t),          intent(in)          :: mesh
            class(expansion_t),     intent(inout)       :: q(:)
            class(expansion_t),     intent(inout)       :: rhs(:)
            integer(ik),            intent(in)          :: ielem
            integer(ik),            intent(in)          :: iface
            integer(ik),            intent(in)          :: iblk
        end subroutine
    end interface


    abstract interface
        subroutine volume_interface(self,mesh,q,rhs,ielem,iblk)
            use mod_kinds,  only: ik
            import equationset_t
            import mesh_t
            import expansion_t

            class(equationset_t),   intent(in)          :: self
            class(mesh_t),          intent(in)          :: mesh
            class(expansion_t),     intent(inout)       :: q(:)
            class(expansion_t),     intent(inout)       :: rhs(:)
            integer(ik),            intent(in)          :: ielem
            integer(ik),            intent(in)          :: iblk
        end subroutine
    end interface

contains

    !===================================================
    !
    !   Given a character string for the variable name,
    !   returns the variable index
    !
    !===================================================
    function get_var(self,varstring) result(varindex)
        class(equationset_t),   intent(in)  :: self
        character(*),           intent(in)  :: varstring

        integer(ik)  :: varindex,ieq

        varindex = 123456789


        ! Search for character string in equation, if found set index
        do ieq = 1,self%neqns
            if (varstring == self%eqns(ieq)%name) then
                varindex = self%eqns(ieq)%ind
                exit
            end if
        end do

        ! Check if index was found
        if (varindex == 123456789) stop "Error: invalid equation string for get_var"

    end function


    !========================================================
    !
    !   Call the boundary flux calculation on each element face
    !
    !========================================================
    subroutine compute_boundary_flux_all(self)
        class(equationset_t),   intent(in)      :: self

    end subroutine

    !========================================================
    !
    !   Call the local lift calculation on each element face
    !
    !========================================================
    subroutine compute_lift_all(self)
        class(equationset_t),   intent(in)      :: self



    end subroutine






end module atype_equationset
