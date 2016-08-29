module bc_lineardiffusion_value
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ME
    use type_bc,            only: bc_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_point,         only: point_t
    use DNAD_D
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: lineardiffusion_value_t



    contains

        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type lineardiffusion_value_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(lineardiffusion_value_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('lineardiffusion_value')


        !
        ! Add functions
        !
        call self%bcproperties%add('Value','Required')         ! add StaticPressure


        !
        ! Add parameters
        !


    end subroutine add_options
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
    subroutine compute(self,worker,prop)
        class(lineardiffusion_value_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: iu

        real(rk)        :: time

        type(AD_D),     allocatable, dimension(:)   :: u_m, flux, integrand, lift, dudx, tmp, tmp2, lift_gq
        real(rk),       allocatable, dimension(:)   :: u_bc, normx
        type(point_t),  allocatable, dimension(:)   :: coords


        !
        ! Get equation index
        !
        iu = prop%get_eqn_index("u")


        !
        ! Get derivative value
        !
        coords = worker%coords()
        time   = worker%time()
        u_bc   = self%bcproperties%compute("Value",time,coords)


        !
        ! Initialize derivatives
        !
        u_m  = worker%interpolate(iu, 'value', ME)
        dudx = worker%interpolate(iu, 'ddx',   ME)


!        tmp = -(u_bc - u_m)
!        tmp2 = matmul(transpose(worker%mesh(worker%face_info%idomain_l)%elems(worker%face_info%ielement_l)%gq%face%val(:,:,worker%face_info%iface)), tmp)
!        lift = matmul(worker%mesh(worker%face_info%idomain_l)%elems(worker%face_info%ielement_l)%invmass, tmp2)
!
!        lift_gq = matmul(worker%mesh(worker%face_info%idomain_l)%elems(worker%face_info%ielement_l)%gq%face%val(:,:,worker%face_info%iface), lift)
!        dudx = dudx + real(6,rk)*lift_gq


        !
        ! Compute diffusive flux
        !
        flux = -dudx


        !
        ! Compute integrand
        !
        normx = worker%normal(1)
        integrand = flux*normx


        call worker%integrate_boundary(iu, integrand)



    end subroutine compute
    !*********************************************************************************************






end module bc_lineardiffusion_value
