module atype_equationset
    use mod_kinds,          only: rk,ik
    use type_equation,      only: equation_t
    use type_mesh,          only: mesh_t
    use atype_solverdata,   only: solverdata_t


    implicit none
    private

    type, public, abstract :: equationset_t
        character(100)              :: name
        integer(ik)                 :: neqns

        ! List equation variables with indices
        type(equation_t), allocatable  :: eqns(:)


        ! Array of boundary flux functions
        !type(boundary_flux_t),  allocatable   :: bnd_adv_flux(:)
        !type(boundary_flux_t),  allocatable   :: bnd_diff_flux(:)

        ! Array of volume flux functions
        !type(volume_flux_t),    allocatable   :: vol_adv_flux(:)
        !type(volume_flux_t),    allocatable   :: vol_diff_flux(:) 

        ! Array of volume source functions
        !type(source_t),         allocatable   :: v_source(:)


    contains
        ! Must define these procedures in the extended type
        procedure(self_interface),     deferred  :: init
        procedure(boundary_interface), deferred  :: compute_boundary_average_flux
        procedure(boundary_interface), deferred  :: compute_boundary_upwind_flux
        procedure(volume_interface),   deferred  :: compute_volume_flux
        procedure(volume_interface),   deferred  :: compute_volume_source


        procedure   :: get_var
        procedure   :: add

    end type equationset_t







    !> Interface definitions
    abstract interface
        subroutine self_interface(self)
            import equationset_t
            class(equationset_t), intent(inout) :: self
        end subroutine
    end interface


    abstract interface
        subroutine boundary_interface(self,mesh,sdata,ielem,iface,iblk)
            use mod_kinds,  only: ik
            import equationset_t
            import mesh_t
            import solverdata_t

            class(equationset_t),   intent(in)          :: self
            class(mesh_t),          intent(in)          :: mesh
            class(solverdata_t),    intent(inout)       :: sdata
            integer(ik),            intent(in)          :: ielem
            integer(ik),            intent(in)          :: iface
            integer(ik),            intent(in)          :: iblk
        end subroutine
    end interface


    abstract interface
        subroutine volume_interface(self,mesh,sdata,ielem,iblk)
            use mod_kinds,  only: ik
            import equationset_t
            import mesh_t
            import solverdata_t

            class(equationset_t),   intent(in)          :: self
            class(mesh_t),          intent(in)          :: mesh
            class(solverdata_t),    intent(inout)       :: sdata
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







    subroutine add(self,varstring,varindex)
        class(equationset_t),   intent(inout)  :: self
        character(*),           intent(in)     :: varstring
        integer(ik),            intent(in)     :: varindex

        type(equation_t), allocatable    :: temp(:)
        integer(ik) :: ieq



        if (allocated(self%eqns)) then
            ! Get allocate temp eqn array with one extra slot for new eqn
            allocate(temp(size(self%eqns) + 1))

            ! Copy current eqns to first temp slots
            do ieq = 1,size(self%eqns)
                temp(ieq) = self%eqns(ieq)
            end do

            ! Add new eqn to last slot
            temp(size(temp))%name = varstring
            temp(size(temp))%ind  = varindex
        else

            allocate(self%eqns(1))
            self%eqns(1)%name = varstring
            self%eqns(1)%ind  = varindex

        end if


    end subroutine









end module atype_equationset
