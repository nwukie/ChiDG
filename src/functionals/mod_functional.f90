module mod_functional
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_fclvector,     only: fclvector_t
    use type_evaluator,     only: evaluator_t
    use mod_string,         only: string_to_lower
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D

    ! INTERNAL FLOW FUNCTIONALS 
    use fcl_mass_averaged_entropy,          only: mass_averaged_entropy_t
    use fcl_mass_averaged_P0_ratio,         only: mass_averaged_P0_ratio_t
    use fcl_mass_averaged_P_ratio,          only: mass_averaged_P_ratio_t
    use fcl_mass_averaged_flowangle12,      only: mass_averaged_flowangle12_t
    use fcl_mass_averaged_flowangle13,      only: mass_averaged_flowangle13_t
    use fcl_mass_averaged_flowangle32,      only: mass_averaged_flowangle32_t
    use fcl_mass_flux,                      only: mass_flux_t
    use fcl_mass_flux_balance,              only: mass_flux_balance_t
    use fcl_mass_averaged_total_pressure,   only: mass_averaged_total_pressure_t
    use fcl_mass_averaged_total_temperature,only: mass_averaged_total_temperature_t

    ! EXTERNAL FLOW FUNCTIONALS
    !use fcl_kinetic_energy,         only: kinetic_energy_t
    use fcl_xforce,                 only: xforce_t
    use fcl_yforce,                 only: yforce_t
    use fcl_zforce,                 only: zforce_t
    use fcl_xforce_coeff,           only: xforce_coeff_t
    use fcl_yforce_coeff,           only: yforce_coeff_t
    use fcl_zforce_coeff,           only: zforce_coeff_t


    ! SCALAR FUNCTIONALS
    use fcl_test_1D_energy,         only: test_1D_energy_t
    use fcl_test_1D_integral,       only: test_1D_integral_t
    use fcl_test_1D_face_integral,  only: test_1D_face_integral_t
    
    implicit none


    ! Register of functionals
    type(fclvector_t)   :: registered_functional
    logical             :: initialized = .false.


contains



    !>  Register functionals in a module vector
    !!  
    !!  This gets called by chidg%start_up('core')
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/10/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine register_functionals()
        integer(ik)     :: nfcls
        integer(ik)     :: ifcl
        
        ! Instantiate functionals
        type(mass_averaged_entropy_t)           :: MASS_AVERAGED_ENTROPY
        type(mass_averaged_P0_ratio_t)          :: MASS_AVERAGED_P0_RATIO
        type(mass_averaged_P_ratio_t)           :: MASS_AVERAGED_P_RATIO
        type(mass_averaged_flowangle12_t)       :: MASS_AVERAGED_FLOWANGLE_12
        type(mass_averaged_flowangle13_t)       :: MASS_AVERAGED_FLOWANGLE_13
        type(mass_averaged_flowangle32_t)       :: MASS_AVERAGED_FLOWANGLE_32
        type(mass_averaged_total_pressure_t)    :: MASS_AVERAGED_TOTAL_PRESSURE
        type(mass_averaged_total_temperature_t) :: MASS_AVERAGED_TOTAL_TEMPERATURE
        !type(kinetic_energy_t)                  :: KINETIC_ENERGY
        type(test_1D_energy_t)                  :: TEST_1D_ENERGY
        type(test_1D_integral_t)                :: TEST_1D_INTEGRAL
        type(test_1D_face_integral_t)           :: TEST_1D_FACE_INTEGRAL
        type(xforce_t)                          :: XFORCE
        type(yforce_t)                          :: YFORCE
        type(zforce_t)                          :: ZFORCE
        type(xforce_coeff_t)                    :: XFORCE_COEFF
        type(yforce_coeff_t)                    :: YFORCE_COEFF
        type(zforce_coeff_t)                    :: ZFORCE_COEFF
        type(mass_flux_t)                       :: MASS_FLUX
        type(mass_flux_balance_t)               :: MASS_FLUX_BALANCE


        if ( .not. initialized ) then
            
            call registered_functional%push_back(MASS_AVERAGED_ENTROPY)
            call registered_functional%push_back(MASS_AVERAGED_P0_RATIO)
            call registered_functional%push_back(MASS_AVERAGED_P_RATIO)
            call registered_functional%push_back(MASS_AVERAGED_FLOWANGLE_12)
            call registered_functional%push_back(MASS_AVERAGED_FLOWANGLE_13)
            call registered_functional%push_back(MASS_AVERAGED_FLOWANGLE_32)
            call registered_functional%push_back(MASS_AVERAGED_TOTAL_PRESSURE)
            call registered_functional%push_back(MASS_AVERAGED_TOTAL_TEMPERATURE)
            !call registered_functional%push_back(KINETIC_ENERGY)
            call registered_functional%push_back(XFORCE)
            call registered_functional%push_back(YFORCE)
            call registered_functional%push_back(ZFORCE)
            call registered_functional%push_back(XFORCE_COEFF)
            call registered_functional%push_back(YFORCE_COEFF)
            call registered_functional%push_back(ZFORCE_COEFF)
            call registered_functional%push_back(ZFORCE_COEFF)
            call registered_functional%push_back(MASS_FLUX)
            call registered_functional%push_back(MASS_FLUX_BALANCE)
            call registered_functional%push_back(TEST_1D_ENERGY)
            call registered_functional%push_back(TEST_1D_INTEGRAL)
            call registered_functional%push_back(TEST_1D_FACE_INTEGRAL)

        end if

        nfcls = registered_functional%size()

        do ifcl = 1,nfcls

            call registered_functional%data(ifcl)%func%init()

        end do

        ! Confirm initialization
        initialized = .true.


    end subroutine register_functionals
    !******************************************************************************************


    

    !>  List all the functionals registered
    !!
    !!  This get called by chidg_edit_adjoint
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/10/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_functionals()
        
        integer(ik)                 :: ifcl, nfcls
        character(:), allocatable   :: name_

        nfcls = registered_functional%size()

        do ifcl = 1,nfcls
            name_=registered_functional%data(ifcl)%func%get_name()
            call write_line(trim(name_))
        end do

    end subroutine list_functionals
    !******************************************************************************************






    !>  Verify that an input functional is registered
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/10/2017
    !!
    !------------------------------------------------------------------------------------------
    function check_functional_existence(func_name) result(exists) 
        character(*),   intent(in)  :: func_name

        logical                     :: exists
        integer(ik)                 :: ifcl, nfcls
        character(:), allocatable   :: name_

        nfcls = registered_functional%size()

        exists = .false.
        do ifcl = 1,nfcls
            name_=registered_functional%data(ifcl)%func%get_name()
            if (name_ == trim(func_name)) then
                exists = .true.
                exit
            end if
        end do

    end function check_functional_existence
    !******************************************************************************************






    !>  Crate a concrete functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/14/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine create_functional(fcl_name,func)
        character(*),                       intent(in)      :: fcl_name
        class(evaluator_t),    allocatable, intent(inout)   :: func

        integer(ik)     :: ierr, ofindex

        if (allocated(func)) deallocate(func)
        ! Find the functional in the register "registered_functional"
        ofindex = registered_functional%index_by_name(trim(fcl_name))
        if (ofindex == 0) call chidg_signal_one(FATAL,"create_functional: functional not recognized", trim(fcl_name))

        ! Allocate concrete functional
        allocate(func, source=registered_functional%data(ofindex)%func, stat=ierr)
        if (ierr/=0) call chidg_signal(FATAL,"create_functional: error allocating functional from register")

        ! Chekc if functional was allocated
        if (.not. allocated(func)) call chidg_signal(FATAL,"create_functional: error allocating the concrete functional")

    end subroutine create_functional
    !******************************************************************************************




end module mod_functional
