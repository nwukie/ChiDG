module bc_lineardiffusion_derivative
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
    type, public, extends(bc_t) :: lineardiffusion_derivative_t



    contains

        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type lineardiffusion_derivative_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(lineardiffusion_derivative_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('lineardiffusion_derivative')


        !
        ! Add functions
        !
        call self%bcproperties%add('Derivative','Required')         ! add StaticPressure


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
        class(lineardiffusion_derivative_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: iu

        real(rk)        :: time

        type(AD_D),     allocatable, dimension(:)   :: flux, integrand
        real(rk),       allocatable, dimension(:)   :: uderiv, normx
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
        uderiv = self%bcproperties%compute("Derivative",time,coords)


        !
        ! Initialize derivatives
        !
        flux = worker%interpolate(iu, 'value', ME)

        normx = worker%normal(1)

        !
        ! Compute diffusive flux
        !
        flux = -uderiv


        !
        ! Compute integrand
        !
        integrand = flux*normx


        call worker%integrate_boundary(iu, integrand)



    end subroutine compute
    !*********************************************************************************************






end module bc_lineardiffusion_derivative
