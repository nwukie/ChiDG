module mod_update_functionals
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NO_DIFF
    use type_timer,                     only: timer_t
    use type_chidg_data,                only: chidg_data_t
    use type_chidg_worker,              only: chidg_worker_t
    use type_chidg_cache,               only: chidg_cache_t
    use type_functional_cache_handler,  only: functional_cache_handler_t
    use type_functional_cache,          only: functional_cache_t
    use mod_io,                         only: verbosity
    implicit none


contains


    !>  Update functionals (i.e. Lift, Drag, Mass averaged entropy, ecc)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   6/17/2017
    !!
    !!  Allowed functionals to be linearized wrt grid-node. Added timer.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/5/2018
    !!
    !!
    !!  TODO: add unsteady and HB (steady only as of now)
    !!
    !!  @param[inout]      data             chidg data
    !!  @param[in]         differentiate    flag for functional derivative computation
    !!
    !---------------------------------------------------------------------------------
    subroutine update_functionals(data,differentiate)
        type(chidg_data_t),     intent(inout),      target      :: data
        integer(ik),            intent(in)                      :: differentiate

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        integer(ik)                 :: ifunc, nfunc
        character(:),   allocatable :: func_name
        type(timer_t)               :: comm_timer, func_timer
        

        ! Initialize Chidg Worker references
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, data%time_manager, cache)


        ! Communicate solution vector
        call comm_timer%start()
        call data%sdata%q%assemble()
        call comm_timer%stop()

        ! Set time info to worker
        worker%itime = data%time_manager%itime
        worker%t     = data%time_manager%t

        
        ! Get overall number of functionals
        nfunc = data%functional_group%n_functionals()


        ! Compute functionals one by one
        do ifunc = 1,nfunc
            
            ! Print functional name being processed
            func_name = data%functional_group%fcl_entities(ifunc)%func%get_name()
            call write_line('-  Updating functional: ',func_name,delimiter=' ',io_proc=GLOBAL_MASTER,silence=(verbosity<3))
            
            ! Update functional
            call func_timer%reset()
            call func_timer%start()
            call update_functional(worker,data,ifunc,differentiate=differentiate) 
            call func_timer%stop()

            ! Print functional computational time
            call write_line('Functionals compute time: ',func_timer%elapsed(),delimiter=' ',io_proc=GLOBAL_MASTER,silence=(verbosity<3))

        end do !i_func


    end subroutine update_functionals
    !********************************************************************************





    
    !>  For each functional/auxiliary initialize a functional_cache_handler and proceede 
    !!  to compute and store the functional/auxiliary
    !!
    !!  N.B. If a functional does not have an auxiliary geometry, this procedure is   
    !!  called anyway and the 'update' subroutine will not do anything.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/3/2017
    !!
    !!
    !!  @param[inout]      worker           chidg worker
    !!  @param[inout]      data             chidg data
    !!  @param[in]         ifunc            functional ID 
    !!  @param[in]         geom_feature     geometry on which the functional needs to 
    !!                                      be compute: reference/auxiliary
    !!  @param[in]         differentiate    flag for functional derivative computation
    !!
    !!
    !!  TODO: add unsteady and HB (steady only as of now)
    !!
    !---------------------------------------------------------------------------------
    subroutine update_functional(worker,data,ifunc,differentiate)
        type(chidg_worker_t),   intent(inout)   :: worker
        type(chidg_data_t),     intent(inout)   :: data
        integer(ik),            intent(in)      :: ifunc
        integer(ik),            intent(in)      :: differentiate

        type(functional_cache_handler_t)     :: functional_cache_handler
        type(functional_cache_t)             :: functional_cache
        
        
        ! Update auxiliary        
        call functional_cache_handler%update(worker,functional_cache,data,ifunc,'auxiliary',differentiate)

        ! Update reference        
        call functional_cache_handler%update(worker,functional_cache,data,ifunc,'reference',differentiate)


    end subroutine update_functional
    !********************************************************************************






end module mod_update_functionals
