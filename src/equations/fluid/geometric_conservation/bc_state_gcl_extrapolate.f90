module bc_state_gcl_extrapolate
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_point,         only: point_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/10/2017
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: gcl_extrapolate_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type gcl_extrapolate_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(gcl_extrapolate_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Geometric Conservation Extrapolate')
        call self%set_family('Extrapolation')


    end subroutine init
    !******************************************************************************************










    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      eqnset  Equation Set type governing the current domain
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !---------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(gcl_extrapolate_t),   intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM

        type(AD_D), allocatable, dimension(:)   ::  &
                g_bar_bc, dgdx_bc, dgdy_bc, dgdz_bc


        !
        ! Get u and grad(u) from face interior to extrapolate
        !
        g_bar_bc = worker%get_field('g_bar','value', 'face interior')
        dgdx_bc  = worker%get_field('g_bar','grad1', 'face interior')
        dgdy_bc  = worker%get_field('g_bar','grad2', 'face interior')
        dgdz_bc  = worker%get_field('g_bar','grad3', 'face interior')


        !
        ! Store as extrpolated boundary condition value and gradient
        !
        call worker%store_bc_state('g_bar',g_bar_bc, 'value')
        call worker%store_bc_state('g_bar',dgdx_bc,  'grad1')
        call worker%store_bc_state('g_bar',dgdy_bc,  'grad2')
        call worker%store_bc_state('g_bar',dgdz_bc,  'grad3')


    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_gcl_extrapolate
