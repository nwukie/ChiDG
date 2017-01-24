module type_pseudo_timestep
    use mod_kinds,          only: ik, rk
    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_properties,    only: properties_t
    implicit none


    !>  A class that handles computing a pseudo timestep for Quasi-Newton
    !!  methods.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !------------------------------------------------------------------------
    type, abstract, public :: pseudo_timestep_t

    contains

        procedure(compute_interface), deferred :: compute

    end type pseudo_timestep_t
    !************************************************************************



    abstract interface
        subroutine compute_interface(self,idomain,mesh,prop,sdata,cfl,itime)
            import ik, rk
            import pseudo_timestep_t
            import mesh_t
            import properties_t
            import solverdata_t
            class(pseudo_timestep_t),   intent(in)      :: self
            integer(ik),                intent(in)      :: idomain
            type(mesh_t),               intent(inout)   :: mesh(:)
            type(properties_t),         intent(in)      :: prop
            type(solverdata_t),         intent(inout)   :: sdata
            real(rk),                   intent(in)      :: cfl(:)
            integer(ik),                intent(in)      :: itime
        end subroutine
    end interface









    !>  Default pseudo_timestep implementation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(pseudo_timestep_t), public :: default_pseudo_timestep_t


    contains

        procedure   :: compute

    end type default_pseudo_timestep_t
    !*************************************************************************








contains




    !>  Default procedure for computing a pseudo timestep.
    !!
    !!  Default is based strictly
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------
    subroutine compute(self,idomain,mesh,prop,sdata,cfl,itime)
        class(default_pseudo_timestep_t),   intent(in)      :: self
        integer(ik),                        intent(in)      :: idomain
        type(mesh_t),                       intent(inout)   :: mesh(:)
        type(properties_t),                 intent(in)      :: prop
        type(solverdata_t),                 intent(inout)   :: sdata
        real(rk),                           intent(in)      :: cfl(:)
        integer(ik),                        intent(in)      :: itime

        integer(ik) :: ielem, ieqn
        real(rk)    ::  h

        
        !
        ! Loop through elements and compute time-step function
        !
       do ielem = 1,mesh(idomain)%nelem

            !
            ! Compute element spacing parameter
            !
            h = mesh(idomain)%elems(ielem)%vol**(1._rk / 3._rk)


            !
            ! Compute elemen-local timestep
            !
            !sdata%dt(idomain,ielem) = cfl*h
            do ieqn = 1,size(cfl)
                mesh(idomain)%elems(ielem)%dtau(ieqn) = cfl(ieqn)*h
            end do


        end do  ! ielem

    end subroutine compute
    !***********************************************************************

















end module type_pseudo_timestep
