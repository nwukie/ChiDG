module sst_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_source_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_source_operator_t
    !******************************************************************************










contains

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_source_operator_t),   intent(inout)      :: self

        ! Set operator name.
        call self%set_name('SST Source Operator')

        ! Set operator type.
        call self%set_operator_type('Volume Diffusive Operator')
        call self%add_auxiliary_field('Wall Distance : p-Poisson')
        ! Set operator equations being integrated.
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')

    end subroutine init
    !********************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_source_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                            &
                                    omega_src, k_src
        real(rk)    :: const, epsilon_vorticity, eps
        

        !
        ! Interpolate solution to quadrature nodes
        !

        omega_src   = worker%get_field('SST Omega Source Term', 'value', 'element')
        k_src       = worker%get_field('SST k Source Term', 'value', 'element')

        


        !========================================================================
        !                       Omega Source Term
        !========================================================================
        call worker%integrate_volume_source('Density * Omega',omega_src)
        
        !========================================================================
        !                       k Source Term
        !========================================================================

        call worker%integrate_volume_source('Density * k',k_src)



    end subroutine compute
    !*********************************************************************************************************






end module sst_source
