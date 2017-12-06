module bc_state_P_neumannpoint
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    use ieee_arithmetic
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: P_neumannpoint_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type P_neumannpoint_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(P_neumannpoint_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('P Neumann Point')
        call self%set_family('Scalar')

    end subroutine init
    !******************************************************************************************












    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(P_neumannpoint_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM


        type(AD_D), allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc


        !
        ! Get 'u' value from face interior to extrapolate
        !
        u_bc = worker%get_field('u', 'value', 'face interior')


        !
        ! Initialize derivative arrays
        !
        dudx_bc = ZERO*worker%get_field('u', 'grad1','face interior')
        dudy_bc = ZERO*dudx_bc
        dudz_bc = ZERO*dudx_bc


        if ((worker%element_info%ielement_g == 1) .and. (worker%iface == 1)) then
        !if (worker%iface == 1) then
           u_bc(1) = 100000._rk
        end if 



        !
        ! Get derivative value
        !
        dudx_bc = self%bcproperties%compute("Grad1",worker%time(),worker%coords())
        dudy_bc = self%bcproperties%compute("Grad2",worker%time(),worker%coords())
        dudz_bc = self%bcproperties%compute("Grad3",worker%time(),worker%coords())


        call worker%store_bc_state('u', u_bc,    'value')
        call worker%store_bc_state('u', dudx_bc, 'grad1')
        call worker%store_bc_state('u', dudy_bc, 'grad2')
        call worker%store_bc_state('u', dudz_bc, 'grad3')






    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_P_neumannpoint
