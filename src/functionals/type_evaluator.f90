module type_evaluator
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_io,                 only: backend
    use type_chidg_worker,      only: chidg_worker_t
    use type_svector,           only: svector_t
    use mod_string,             only: string_t
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
    use type_functional_cache,  only: functional_cache_t
    
    use DNAD_D

    implicit none

    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/10/2017
    !!
    !!  Tested in io/test_hdfio_functional.pf
    !!
    !---------------------------------------------------------------------------------------------
    type, public, abstract  :: evaluator_t

        type(string_t)      :: name
        type(string_t)      :: evaluator_type           ! functional (May include bcs in the future)
        type(string_t)      :: integration_type         ! "VOLUME INTEGRAL" or "SURFACE INTEGRAL"
        type(svector_t)     :: reference_geom           ! List of reference geometries
        type(svector_t)     :: auxiliary_geom           ! List of auxiliary geometries (if any)
        integer(ik)         :: ID                       ! Identification number

        type(svector_t)     :: intermediate_integrals
        
        !type(zone_t),               :: zone_ID

    contains

       
        procedure, public   :: get_name
        procedure, public   :: set_name
       
        procedure, public   :: set_int_type
        procedure, public   :: get_int_type
        
        procedure, public   :: set_eval_type
        procedure, public   :: get_eval_type
        
        procedure, public   :: set_ref_geom
        procedure, public   :: get_ref_geom
        
        procedure, public   :: set_aux_geom
        procedure, public   :: get_aux_geom
        
        procedure, public   :: n_ref_geom
        procedure, public   :: n_aux_geom
       

        ! Functional overwritten procedures
        procedure, public   :: init
        procedure, public   :: check
        procedure, public   :: compute_functional
        procedure, public   :: compute_auxiliary
        procedure, public   :: finalize_functional
        procedure, public   :: finalize_auxiliary
        procedure, public   :: store_value 
        procedure, public   :: store_deriv

        procedure, public   :: add_integral
        

    end type evaluator_t
    !*********************************************************************************************




contains

    
    !> Set the evaluator name
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine set_name(self,eval_name)
        class(evaluator_t),     intent(inout)   :: self
        character(*),           intent(in)      :: eval_name

        call self%name%set(trim(eval_name))

    end subroutine set_name
    !*********************************************************************************************






    !> Get the evaluator name
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function get_name(self) result(eval_name)
        class(evaluator_t),    intent(in)  :: self
        
        character(:), allocatable  :: eval_name

        eval_name = trim(self%name%get())

    end function get_name
    !*********************************************************************************************





    !>  Set the evaluator type
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !!
    !------------------------------------------------------------------------------
    subroutine set_eval_type(self,name_e)
        class(evaluator_t), intent(inout)   :: self
        character(*),       intent(in)      :: name_e

        call self%evaluator_type%set(trim(name_e))

    end subroutine set_eval_type
    !------------------------------------------------------------------------------





    !>  Get evaluator type
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function get_eval_type(self) result(type_e)
        class(evaluator_t),             intent(in)  :: self
        character(:),   allocatable :: type_e

        type_e = self%evaluator_type%get()

    end function get_eval_type
    !------------------------------------------------------------------------------






    !>  Set type of integration used by the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    subroutine set_int_type(self,integ_type)
        class(evaluator_t), intent(inout)   :: self
        character(*),       intent(in)      :: integ_type

        character(:),   allocatable     :: usr_msg

        select case (trim(integ_type))
            case ('VOLUME INTEGRAL','volume','VOL','volume integral','Volume Integral')

                call self%integration_type%set('VOLUME INTEGRAL')


            case ('SURFACE INTEGRAL','surface','SURF','surface integral','Surface Integral',&
                   'Face Integral','FACE INTEGRAL')

                call self%integration_type%set('FACE INTEGRAL')

            case default

                usr_msg  = "Wrong integration type. Please, check the implementation of the current evaluator"
                
                call chidg_signal(FATAL,usr_msg)

        end select

    end subroutine set_int_type
    !------------------------------------------------------------------------------






    !>  Get type of integration used by the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function get_int_type(self) result(type_e)
        class(evaluator_t),             intent(in)  :: self

        character(:),   allocatable     :: type_e
        character(:),   allocatable     :: user_msg

        type_e = self%integration_type%get()
        

    end function get_int_type
    !------------------------------------------------------------------------------






    !>  Set reference geometry of the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !!  @param[IN]  input_geom    Comma-separated list of reference geometries coming from
    !!                      'chidg edit'
    !!
    !------------------------------------------------------------------------------
    subroutine set_ref_geom(self,input_geom)
        class(evaluator_t), intent(inout)   :: self
        character(*),       intent(in)      :: input_geom

        type(string_t)                  :: full_string
        type(string_t), allocatable     :: geom_list(:)
        integer(ik)                     :: igeom

        ! Strip the comma-separated list
        call full_string%set(trim(input_geom))
        geom_list = full_string%strip(',')


        ! Add each geoemtry to self%ref_geom
        do igeom = 1,size(geom_list)
            
            ! Check that "none" is not added to reference geometries
            if (geom_list(igeom)%get() /= "none") then
                call self%reference_geom%push_back_unique(geom_list(igeom))
            end if

        end do


    end subroutine set_ref_geom
    !------------------------------------------------------------------------------



    


    !>  Get reference geometry of the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function get_ref_geom(self,geom_index) result(geom_name)
        class(evaluator_t), intent(in)              :: self
        integer(ik),        intent(in), optional    :: geom_index

        type(string_t)                :: str_name
        character(len=:), allocatable :: geom_name

        if (present(geom_index)) then
        
            str_name  = self%reference_geom%at(geom_index)
            geom_name = str_name%get()

        else
            
            str_name  = self%reference_geom%concatenate(",") 
            geom_name = str_name%get()
            ! If no reference is presenet, return 'none'
            if (geom_name == '') then
                geom_name = "none"
            else
            end if

        end if

    end function get_ref_geom
    !------------------------------------------------------------------------------




    
    
    
    !>  Set auxiliary geometry of the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !!  @param[IN]  input_geom    Comma-separated list of auxiliary geometries coming from
    !!                      'chidg edit'
    !!
    !------------------------------------------------------------------------------
    subroutine set_aux_geom(self,input_geom)
        class(evaluator_t), intent(inout)   :: self
        character(*),       intent(in)      :: input_geom

        type(string_t)                  :: full_string
        type(string_t), allocatable     :: geom_list(:)
        integer(ik)                     :: igeom

        ! Strip the comm:a-separated list
        call full_string%set(trim(input_geom))
        geom_list = full_string%strip(',')


        ! Add each geoemtry to self%ref_geom
        do igeom = 1,size(geom_list)

            ! Check that "none" is not added to reference geometries
            if (geom_list(igeom)%get() /= "none") then
                call self%auxiliary_geom%push_back_unique(geom_list(igeom))
            end if
            
        end do


    end subroutine set_aux_geom
    !------------------------------------------------------------------------------






    !>  Get auxiliary geometry of the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function get_aux_geom(self,geom_index) result(geom_name)
        class(evaluator_t), intent(in)           :: self
        integer(ik),        intent(in), optional :: geom_index
        
        type(string_t)                :: str_name
        character(len=:), allocatable :: geom_name

        if (present(geom_index)) then
        
            str_name  = self%auxiliary_geom%at(geom_index)
            geom_name = str_name%get()

        else
            
            str_name  = self%auxiliary_geom%concatenate(",") 
            geom_name = str_name%get()
            ! If no auxiliary is presenet, return 'none'
            if (geom_name == '') then
                geom_name = "none"
            else
            end if

        end if

    end function get_aux_geom
    !------------------------------------------------------------------------------






    !>  Get number of reference geometries registered with the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function n_ref_geom(self) result(n_geom)
        class(evaluator_t), intent(in)  :: self
        
        integer(ik) :: n_geom

        n_geom = self%reference_geom%size()

    end function n_ref_geom
    !------------------------------------------------------------------------------






    !>  Get number of auxiliary geometries registered with the evaluator
    !!
    !!  @author Matteo Ugolotti
    !!  @date   02/12/2017
    !!
    !------------------------------------------------------------------------------
    function n_aux_geom(self) result(n_geom)
        class(evaluator_t), intent(in)  :: self
        
        integer(ik) :: n_geom

        n_geom = self%auxiliary_geom%size()

    end function n_aux_geom
    !------------------------------------------------------------------------------







    !>  This is overwritten by individual functional and is meant to initialize the functional
    !!  with fundamental general information
    !!
    !!  It does nothing if it is not overwritten.
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   12/2/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self)
        class(evaluator_t),    intent(inout)  :: self


    end subroutine init
    !*********************************************************************************************






    !>  Check the correctness of reference and auxiliary geometries. Verify that the user's input
    !!  satisfies the requirements 
    !!
    !!  It does nothing if it is not overwritten.
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine check(self)
        class(evaluator_t),    intent(inout)  :: self


    end subroutine check
    !*********************************************************************************************





    !>  Computation of the funcitonal on the reference geometry. Each functional will overwrite this procedure 
    !!  This is used to compute the necessary integrals on each element/face belonging to the reference
    !!  geometry
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(evaluator_t),         intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(functional_cache_t),   intent(inout)   :: cache


    end subroutine compute_functional
    !*********************************************************************************************





    !>  Computation of the functional on the auxiliary geometry. Some functional might need to overwrite
    !!  this procedure 
    !!
    !!  This might not be overwritten if the functional does not require any auxiliary computation
    !!
    !!  This is used to compute the necessary integrals on each element/face belonging to the auxiliary
    !!  geometry
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine compute_auxiliary(self,worker,cache)
        class(evaluator_t),         intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(functional_cache_t),   intent(inout)   :: cache


    end subroutine compute_auxiliary
    !*********************************************************************************************





    !>  This procedure is complmentary to "compute_auxiliary". 
    !!  Once all the singluar integral on each face/element are computed, this subroutine gives the 
    !!  user the opportunity to compute the operations on the overall integral on the auxiliary
    !!  geometry 
    !!  
    !! 
    !!  For instance mass-averaged pressure needs to be computed on the auxiliary geoemtry to eventually
    !!  compute a pressure ratio. The pressure flux can be divided here by the mass flux over the 
    !!  entire auxiliary geometry:
    !!
    !!          ma_pressure = (A*rho*vn*P)/(A*rho*vn)
    !!
    !!  This might not be overwritten if the auxiliary geoemtry is not needed.
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_auxiliary(self,worker,cache)
        class(evaluator_t),         intent(inout)   :: self
        type(chidg_worker_t),       intent(in)      :: worker
        type(functional_cache_t),   intent(inout)   :: cache


    end subroutine finalize_auxiliary
    !*********************************************************************************************






    !>  This procedure is complmentary to "compute_functional". 
    !!  Once all the singluar integral on each face/element are computed, this subroutine gives the 
    !!  user the opportunity to compute the operations on the overall integral on the reference
    !!  geometry 
    !!  
    !! 
    !!  For instance mass-averaged entropy functional will need to divide the overall entropy
    !!  flux by the mass flux over the entire reference geometry:
    !!
    !!          ma_entropy = (A*rho*vn*S)/(A*rho*vn)
    !!  
    !!  Since the functional cache is available here, one could also use integrals evaluated on the
    !!  auxiliary geometry for further reduction of the functional computed on the reference geometry
    !!  For instance, for pressure ratio functional, you might want here to divide the mass-averaged
    !!  pressure computed on the reference geometry (obtained in this same procedure) and divide it 
    !!  by the mass averaged pressure computed on the auxiliary geometry (obtained in finalize auxiliary)  
    !!
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_functional(self,worker,cache) 
        class(evaluator_t),         intent(inout)   :: self
        type(chidg_worker_t),       intent(in)      :: worker
        type(functional_cache_t),   intent(inout)   :: cache


    end subroutine finalize_functional
    !*********************************************************************************************

    


    
    !>  Used to store the final functional value
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function store_value(self,cache) result(res) 
        class(evaluator_t),         intent(inout)   :: self
        type(functional_cache_t),   intent(inout)   :: cache
        
        real(rk)    :: res 

    end function store_value
    !*********************************************************************************************





    !>  Used to store the final functional derivatives
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function store_deriv(self,cache) result(res) 
        class(evaluator_t),         intent(inout)   :: self
        type(functional_cache_t),   intent(inout)   :: cache
        
        type(chidg_vector_t)    :: res 

        res = chidg_vector(trim(backend))

    end function store_deriv
    !*********************************************************************************************







    !>  Add intermediate integral to the functional. This will trigger allocation in the
    !!  integral_cache.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/16/2020
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_integral(self,integral_name)
        class(evaluator_t), intent(inout)   :: self
        character(*),       intent(in)      :: integral_name

        call self%intermediate_integrals%push_back_unique(string_t(trim(integral_name)))

    end subroutine add_integral
    !*********************************************************************************************



end module type_evaluator
