!*************************************************************************
!* Dual Number Automatic Differentiation (DNAD) of Fortran Codes 
!*----------------------------------------------------------------
!* COPYRIGHT (c) Wenbin Yu, All rights reserved, you are free to copy, 
!* modify or translate this code to other languages such as c/c++. If 
!* you find a bug please let me know through wenbinyu.heaven@gmail.com. If 
!* you added new functions and want to share with others, please let me know too.
!* You are welcome to share your successful stories with us through 
!* http://groups.google.com/group/hifi-comp. 
!******************************************************************************
!* Simple Instruction:
!*---------------------
!* This is the procedure file: define the procedures (functions/subroutines) 
!* needed for overloading intrinsic functions and operators. If a function 
!* or operation is not defined (unlikely), you need to create a corresponding 
!* interface in header file (DNADHeaders.f90) and a corresponding 
!* function/subroutine in DNAD.f90.
!* If you are uncertain whether I did things in the correct way or not, you can 
!* create a new module with the same definition of the type AD_D, but delete
!* all the interfaces and procedures first. Then compile it with your analysis 
!* codes after you changed all the definitions of real numbers to be dual numbers. 
!* The compiler will complain about some operations and functions are not defined. 
!* Then you do it one by one by looking at my original module, copying the needed
!* interfaces/procedures into your module. This way, you can also get a leaner
!* model and use those only needed for your analysis code.
!*********************************************************************************
!* Acknowledgements
!*-------------------
!* The development of DNAD is supported, in part, by the Chief Scientist 
!* Innovative Research Fund at AFRL/RB WPAFB, and by Department of Army 
!* SBIR (Topic A08-022) through Advanced Dynamics Inc. The views and 
!* conclusions contained herein are those of the authors and should not be 
!* interpreted as necessarily representing the official policies or endorsement,
!* either expressed or implied, of the funding agency.    
!**********************************************************************************
!* Citation
!*-------------------
!* Your citation of the following two papers are appreciated: 
!* Yu, W. and Blair, M.: "DNAD, a Simple Tool for Automatic Differentiation of 
!* Fortran Codes Using Dual Numbers," Computer Physics Communications, vol. 184, 
!* 2013, pp. 1446-1452. 
!* 
!* Spall, R. and Yu, W.: "Imbedded Dual-Number Automatic Differentiation for CFD 
!* Sensitivity Analysis," Journal of Fluids Engineering, vol. 135, 2013, 014501.
!**********************************************************************************
!* Basic idea
!*-------------------------------------
!* Carry out Automatic Differentiation using the arithmetic of dual numbers. 
!* It is supposed to be more efficient and more accurate than complex-step 
!* approximation. What one needs to do is to create a new data type according  
!* to the rules of dual numbers, which is equivalent to create a quasi-complex
!* numbers. And overload all the intrinsic operators and functions for
!* such a type.
!* For simplicity and memoric convenience, we call a dual number as a dual (D), 
!* a double precision real number as a real (R) and a single precision real 
!* number as a single (S), and an integer number as an integer (I). We also 
!* use M to denote a 2D array, V to denote a vector, M3 to denote a 3D arry. 
!* By defining functions/subroutines to the ELEMENTAL, the inputs/outputs 
!* can be automatically overloaded for vectors and matrices of the same shape
!**********************************************************************************
!* To AD a Fortran code use DNAD
!* Step 0: compile DNAD.f90 to be DNAD.o
!* Step 1: Replace all the definitions of real numbers with dual numbers
!*        For example
!*              replace REAL(8) :: x     
!*               with  TYPE(AD_D) :: x
!*              replace REAL(8), PARAMETER:: ONE=1.0D0
!*               with TYPE(AD_D),PARAMETER::ONE=AD_D(1.0D0,0.D0)
!* Step 2: Insert USE Dual_Num_Auto_Diff right after Module/Function/Subroutine/Program
!*         statements used TYPE(AD_D)
!* Step 3: Change IO commands correspondingly if the code does not use free formatting 
!*         read and write (can be automated by written some general-purpose utility subroutines)
!* Step 4: Recompile the source along with DNAD.o
!* The whole process can be automated, and even manually it only takes just a few minutes 
!* for most real analysis codes, although step 3 is code dependent.
!*****************************************************************************************
!* Example
!* The analysis version of the code computing the area of a circle with input radius
!* PROGRAM CircleArea
!*    REAL(8),PARAMETER:: PI=3.141592653589793D0
!*    REAL(8):: radius, area
!*    READ(*,*) radius
!*    Area=PI*radius**2
!*    WRITE(*,*) "AREA=", Area
!* END PROGRAM CircleArea
!*Input: 5
!*Output: AREA=78.5398163397448
!*---------------------------------------------------------------------------------------
!* The AD version of the code computing the area of a circle and sensitivity with input radius
!*PROGRAM CircleArea
!*    USE DNAD
!*    TYPE (AD_D),PARAMETER:: PI=AD_D(3.141592653589793D0,0.D0)
!*    TYPE (AD_D):: radius,area
!*    READ(*,*) radius
!*    Area=PI*radius**2
!*    WRITE(*,*) "AREA=",Area
!*END PROGRAM CircleArea
!* Input 5, 1
!* Output: AREA=78.5398163397448, 31.4159265358979
!*--------------------------------------------------------
!*****************************************************************************************
!* Possible Mistakes to Avoid
!* 11/23/2011: DNAD should always compute the same functional value. However, in very rare situations, 
!* sensitivities are not computed as expected. One such case has been identified with 
!* dgemv in the blas package. The main reason is that dgemv skip some calculations due to
!* the fact of multiplying zero. However, sometimes, the functional value is zero, but 
!* the sensitivity is not zero. In DNAD the comparision is done only between functional values
!* hence is sensitivity calculation is also avoid altogether. For example, the following piece of code
!*                IF( X( JX ).NE.ZERO )THEN
!*                  TEMP = ALPHA*X( JX )
!*                  DO 50, I = 1, M
!*                     Y( I ) = Y( I ) + TEMP*A( I, J )
!*   50             CONTINUE
!*              END IF
!* the solution is to either to rewrite dgemv or assign epsilon(0.0d0) to x(jx) if the functional value
!* while the corresponding derivatives are not zero. 
!************************************************************************************************************

MODULE DNAD_D
    use mod_kinds,  only: rsingle,rdouble,ishort,ilong,rk
IMPLICIT NONE

! Dynamic allocation, so this is no longer used
!INTEGER(2), PUBLIC,PARAMETER:: NDV_AD=4 ! number of design variables

PRIVATE
! Imported kinds from common file for consistency
integer, parameter :: SNG_AD = rsingle
integer, parameter :: DBL_AD = rdouble
!INTEGER, PARAMETER:: SNG_AD=SELECTED_REAL_KIND(6)  ! single precision
!INTEGER, PARAMETER:: DBL_AD=SELECTED_REAL_KIND(15) ! double precision
REAL(DBL_AD)      ::negative_one=-1.0d0

    TYPE,PUBLIC:: AD_D
    !TYPE :: AD_D
                                    ! make this private will create difficulty to use the original write/read commands,
                                   ! hence x_ad_ and xp_ad_ are variables which can be accessed using D%x_ad_ and
                                    ! D%xp_ad_ in other units using this module in which D is defined as TYPE DUAL_NUM.
        REAL(DBL_AD)              :: x_ad_      ! functional value
        REAL(DBL_AD), allocatable :: xp_ad_(:)  ! derivative


    END TYPE AD_D


    interface AD_D
        module procedure new_dual_num_s
        module procedure new_dual_num_l
    end interface


INCLUDE 'DNADHeaders_D.f90'

 
CONTAINS
!=======================
!
!   Contstructors
!
!=======================
    ! SHORT INTEGER
    type(AD_D) function new_dual_num_s(nd)
        integer(ishort), intent(in) :: nd
        integer(ishort)    :: ierr

        allocate(new_dual_num_s%xp_ad_(nd),stat=ierr)
        if (ierr /= 0) stop "Allocation error: DNAD_D - new_dual_num"

        new_dual_num_s%xp_ad_ = real(0.,DBL_AD)
    end function

    ! LONG INTEGER
    type(AD_D) function new_dual_num_l(nd)
        integer(ilong), intent(in) :: nd
        integer(ishort)    :: ierr

        allocate(new_dual_num_l%xp_ad_(nd),stat=ierr)
        if (ierr /= 0) stop "Allocation error: DNAD_D - new_dual_num"

        new_dual_num_l%xp_ad_ = real(0.,DBL_AD)
    end function




!*********Begin: functions/subroutines for overloading operators    

!******* Begin: (=)
!--------------------- 
    !** dual=integer:     ! <u,up>=<a,0.D0>
    ELEMENTAL SUBROUTINE ASSIGN_DI_D(u,n)
         TYPE (AD_D), INTENT(INOUT) :: u
         INTEGER, INTENT(IN)::n

 !         u=AD_D(n,0.d0)    ! It is shown that this approach is much slower than the new approach
         u%x_ad_= REAL(n,DBL_AD) ! REAL(n,DBL_AD) is faster than let the code do the conversion.
         u%xp_ad_=0.0D0

    END SUBROUTINE ASSIGN_DI_D

    !** dual=real: <u,up>=<a,0.D0>
    ELEMENTAL SUBROUTINE ASSIGN_DR_D(u,n)
         TYPE (AD_D), INTENT(INOUT)::u
         REAL(DBL_AD), INTENT(IN)::n

          u%x_ad_=n
          u%xp_ad_=0.0D0

    END SUBROUTINE ASSIGN_DR_D

    !** dual=single
    ! <u,up>=<n,0.D0>
    ELEMENTAL SUBROUTINE ASSIGN_DS_D(u,n)
         TYPE (AD_D), INTENT(INOUT)::u
         REAL(SNG_AD), INTENT(IN)::n

         u%x_ad_=REAL(n,DBL_AD)
         u%xp_ad_=0.0D0

    END SUBROUTINE ASSIGN_DS_D



    !** integer=dual
    ! n=u%x_ad_
    ELEMENTAL SUBROUTINE ASSIGN_ID_D(n,u)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER, INTENT(OUT)::n

         ! Wukie
         n=int(u%x_ad_)
         !n=u%x_ad_

    END SUBROUTINE ASSIGN_ID_D

!******* End: (=)
!--------------------- 


!******* Begin: (+)
!--------------------- 


    !dual=+<v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res=u ! It is faster than assigning component wise
      
    END FUNCTION ADD_D_D

 
    ! <u,up>+<v,vp>=<u+v,up+vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u,v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = u%x_ad_+v%x_ad_ 
         res%xp_ad_ = u%xp_ad_+v%xp_ad_

    END FUNCTION ADD_DD_D

    ! dual+integer
    ! <v,vp>+<n,0>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_DI_D(v,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::v
         INTEGER,INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
            

         res%x_ad_  = REAL(n,DBL_AD)+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_DI_D
    
    ! dual +real
    ! <v,vp>+<n,0>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_DR_D(v,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::v
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = n+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_DR_D

    !-----------------------------------------
    ! dual+single
    ! <v,vp>+<n,0>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_DS_D(v,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::v
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = REAL(n,DBL_AD)+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_DS_D
    !-----------------------------------------


    !-----------------------------------------
    !an integer+ dual number
    ! <n,0>+<v,vp>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_ID_D(n,v) RESULT(res)
         INTEGER,INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = REAL(n,DBL_AD)+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_ID_D

   
    !-----------------------------------------
    ! real + dual 
    ! <n,0>+<v,vp>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_RD_D(n,v) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = n+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_RD_D

     ! single + dual 
    ! <n,0>+<v,vp>=<n+v,vp>
    !-----------------------------------------
    ELEMENTAL FUNCTION ADD_SD_D(n,v) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         res%x_ad_  = REAL(n,DBL_AD)+v%x_ad_ 
         res%xp_ad_ = v%xp_ad_

    END FUNCTION ADD_SD_D

!******* End: (+)
!--------------------- 


!******* Begin: (-)
!--------------------- 


    !-----------------------------------------
    ! negate a dual number
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = -u%x_ad_ 
         res%xp_ad_= -u%xp_ad_
    
    END FUNCTION MINUS_D_D

    !-----------------------------------------
    ! <u,up>-<v,vp>=<u-v,up-vp>
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u,v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = u%x_ad_-v%x_ad_ 
         res%xp_ad_= u%xp_ad_-v%xp_ad_
    
    END FUNCTION MINUS_DD_D

    !-----------------------------------------
    ! dual number - integer
    ! <u,up>-<n,0>=<u-n,up>
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_DI_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER,INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = u%x_ad_-REAL(n,DBL_AD) 
         res%xp_ad_= u%xp_ad_
    
    END FUNCTION MINUS_DI_D

    !-----------------------------------------
    ! dual number - real
    ! <u,up>-<n,0>=<u-n,up>
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_DR_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = u%x_ad_-n
         res%xp_ad_= u%xp_ad_
    
    END FUNCTION MINUS_DR_D

    !-----------------------------------------
    ! dual number - single
    ! <u,up>-<n,0>=<u-n,up>
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_DS_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = u%x_ad_-REAL(n,DBL_AD) 
         res%xp_ad_= u%xp_ad_
    
    END FUNCTION MINUS_DS_D


    !-----------------------------------------
    ! integer-dual number 
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_ID_D(n,u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER,INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = REAL(n,DBL_AD)-u%x_ad_ 
         res%xp_ad_= -u%xp_ad_
    
    END FUNCTION MINUS_ID_D

    !-----------------------------------------
    ! real-dual number
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_RD_D(n,u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ = n-u%x_ad_ 
         res%xp_ad_= -u%xp_ad_
    
    END FUNCTION MINUS_RD_D

    !-----------------------------------------
    ! single-dual number  
    ! <n,0>-<u,up>=<n-u,-up>
    !-------------------------------------------------
    ELEMENTAL FUNCTION MINUS_SD_D(n,u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
    
         res%x_ad_ =REAL(n,DBL_AD) - u%x_ad_
         res%xp_ad_=- u%xp_ad_
    
    END FUNCTION MINUS_SD_D

!******* END: (-)
!--------------------- 


!******* BEGIN: (*)
!--------------------- 

    !-----------------------------------------
    ! <u,up>*<v,vp>=<u*v,up*v+u*vp>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u,v
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = u%x_ad_*v%x_ad_ 
         res%xp_ad_= u%xp_ad_*v%x_ad_ + u%x_ad_*v%xp_ad_
 
    END FUNCTION MULT_DD_D


    !-----------------------------------------
    !  dual*integer
    ! <u,up>*<n,0>=<u*n,up*n>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_DI_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER,INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = REAL(n,DBL_AD)*u%x_ad_ 
         res%xp_ad_= REAL(n,DBL_AD)*u%xp_ad_

    END FUNCTION MULT_DI_D

    !-----------------------------------------
    !  dual*real
    ! <u,up>*<n,0>=<u*n,up*n>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_DR_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = n*u%x_ad_ 
         res%xp_ad_= n*u%xp_ad_

    END FUNCTION MULT_DR_D


    !-----------------------------------------
    !  dual*single
    ! <u,up>*<n,0>=<u*n,up*n>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_DS_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = REAL(n,DBL_AD)*u%x_ad_ 
         res%xp_ad_= REAL(n,DBL_AD)*u%xp_ad_

    END FUNCTION MULT_DS_D

    
    !-----------------------------------------
    ! integer*dual
    ! <n,0>*<v,vp>=<n*v,n*vp>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_ID_D(n,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::v
         INTEGER,INTENT(IN)::n
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))

  
         res%x_ad_ = REAL(n,DBL_AD)*v%x_ad_ 
         res%xp_ad_= REAL(n,DBL_AD)*v%xp_ad_

    END FUNCTION MULT_ID_D

    !-----------------------------------------
    ! real* dual 
    ! <n,0>*<u,up>=<u*n,up*n>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_RD_D(n,u) RESULT(res)
!    FUNCTION MULT_RD_D(n,u) RESULT(res)
         REAL(DBL_AD),  INTENT(IN)  ::n
         TYPE(AD_D),    INTENT(IN)  ::u
         TYPE(AD_D) ::res

         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = n*u%x_ad_
         res%xp_ad_= n*u%xp_ad_

    END FUNCTION MULT_RD_D


!
!    FUNCTION MULT_RD_D(n,u) RESULT(res)
!         REAL(DBL_AD),  INTENT(IN)  ::n
!         TYPE(AD_D),    INTENT(IN)  ::u(:)
!         TYPE(AD_D), allocatable                 ::res(:)
!
!        integer :: i
!!         allocate(res(size(u)))
!         res = u
!         res = 0._DBL_AD
!!         allocate(res%xp_ad_(size(u%xp_ad_)))
!
!         res%x_ad_ = n*u%x_ad_
!
!         do i = 1,size(res)
!            res(i)%xp_ad_= n*u(i)%xp_ad_
!        end do
!
!    END FUNCTION MULT_RD_D



    !-----------------------------------------
    ! MULTIPLY a dual number with REAL number
    ! <n,0>*<u,up>=<u*n,up*n>
    !----------------------------------------
    ELEMENTAL FUNCTION MULT_SD_D(n,u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL(SNG_AD),INTENT(IN)::n
         allocate(res%xp_ad_(size(u%xp_ad_)))
        
         res%x_ad_ = REAL(n,DBL_AD)*u%x_ad_ 
         res%xp_ad_= REAL(n,DBL_AD)*u%xp_ad_
        
    END FUNCTION MULT_SD_D


!******* END: (*)
!--------------------- 


!******* BEGIN: (/)
!--------------------- 

    !-----------------------------------------
    ! <u,up>/<v,vp>=<u/v,(up-u vp/v)/v>
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u,v
         REAL(DBL_AD)::tmp 
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         tmp=1.D0/v%x_ad_
         res%x_ad_ = u%x_ad_*tmp
         res%xp_ad_ =(u%xp_ad_- res%x_ad_*v%xp_ad_)*tmp
    
    END FUNCTION DIV_DD_D

    !-----------------------------------------
    ! <u,up>/<n,0>=<u/n,up/n>
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_DI_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER,INTENT(IN)::n
         REAL(DBL_AD)::tmp 
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         tmp=1.D0/REAL(n,DBL_AD)
         res%x_ad_ = u%x_ad_*tmp
         res%xp_ad_ =u%xp_ad_*tmp
     
    END FUNCTION DIV_DI_D

    !-----------------------------------------
    ! DIVIDE dual number with respect to real numbers
    ! <u,up>/<n,0>=<u/n,up/v>
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_DR_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D):: res
         REAL(DBL_AD)::tmp 
         allocate(res%xp_ad_(size(u%xp_ad_)))
         
         tmp=1.0D0/n
         res%x_ad_ = u%x_ad_*tmp
         res%xp_ad_ =u%xp_ad_*tmp
     
    END FUNCTION DIV_DR_D

  
    !-----------------------------------------
    ! DIVIDE dual number with respect to single
    ! <u,up>/<n,0>=<u/n,up/v>
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_DS_D(u,n) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D):: res
         REAL(DBL_AD)::tmp 
         allocate(res%xp_ad_(size(u%xp_ad_)))
         
         tmp=1.0D0/REAL(n,DBL_AD)
         res%x_ad_ = u%x_ad_*tmp
         res%xp_ad_ =u%xp_ad_*tmp
     
    END FUNCTION DIV_DS_D
    
    !-----------------------------------------
    ! integer/<v,vp>
    ! <n,0>/<v,vp>=<n/v,-n vp/v/v>=n*<1,-vp/v>/v
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_ID_D(n,v) RESULT(res)
         INTEGER,INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         tmp=1.D0/v%x_ad_
         res%x_ad_=REAL(n,DBL_AD)*tmp
         res%xp_ad_=-res%x_ad_*tmp*v%xp_ad_
     
    END FUNCTION DIV_ID_D


    !-----------------------------------------
    ! real/<v,vp>
    ! <n,0>/<v,vp>=<n/v,-n vp/v/v>=n*<1,-vp/v>/v
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_RD_D(n,v) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         tmp=1.D0/v%x_ad_
         res%x_ad_=n*tmp
         res%xp_ad_=-res%x_ad_*tmp*v%xp_ad_
     
    END FUNCTION DIV_RD_D

    !-----------------------------------------
    ! single/<v,vp>
    ! <n,0>/<v,vp>=<n/v,-n vp/v/v>=n*<1,-vp/v>/v
    !----------------------------------------
    ELEMENTAL FUNCTION DIV_SD_D(n,v) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::n
         TYPE (AD_D), INTENT(IN)::v
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(v%xp_ad_)))
         
         tmp=1.D0/v%x_ad_
         res%x_ad_=REAL(n,DBL_AD)*tmp
         res%xp_ad_=-res%x_ad_*tmp*v%xp_ad_
     
    END FUNCTION DIV_SD_D

!******* END: (/)
!--------------------- 


!******* BEGIN: (**)
!--------------------- 

    !-----------------------------------------
    ! POWER dual numbers
    ! <u,up>^k=<u^k,k u^{k-1} up>
    !----------------------------------------
    ELEMENTAL FUNCTION POW_I_D(u,k) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         INTEGER,INTENT(IN)::k
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
         
         tmp=u%x_ad_**(k-1)
         res%x_ad_ = u%x_ad_*tmp 
         res%xp_ad_=REAL(k,DBL_AD)*tmp*u%xp_ad_

    END FUNCTION POW_I_D

    !-----------------------------------------
    ! POWER dual numbers
    ! <u,up>^k=<u^k,k u^{k-1} up>
    !----------------------------------------
    ELEMENTAL FUNCTION POW_R_D(u,k) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD),INTENT(IN)::k
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
         
         tmp=u%x_ad_**(k-1)
         res%x_ad_ = u%x_ad_*tmp 
         res%xp_ad_=k*tmp*u%xp_ad_

    END FUNCTION POW_R_D

    !-----------------------------------------
    ! POWER dual numbers
    ! <u,up>^k=<u^k,k u^{k-1} up>
    !----------------------------------------
    ELEMENTAL FUNCTION POW_S_D(u,k) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(SNG_AD),INTENT(IN)::k
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
         
         tmp=u%x_ad_**(k-1)
         res%x_ad_ = u%x_ad_*tmp 
         res%xp_ad_=k*tmp*u%xp_ad_

    END FUNCTION POW_S_D
    
    !-----------------------------------------
    ! POWER dual numbers to a dual power
    ! <u,up>^(v,vp)=<u^v,u^v (v/u*up+Log(u)*vp>
    !----------------------------------------
    ELEMENTAL FUNCTION POW_D_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D), INTENT(IN)::v
         REAL(DBL_AD)::uf,vf
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
         uf=u%x_ad_
         vf=v%x_ad_
         
         res%x_ad_ =uf**vf
         res%xp_ad_=res%x_ad_*(vf/uf*u%xp_ad_+LOG(uf)*v%xp_ad_)

    END FUNCTION POW_D_D
!******* END: (**)
!--------------------- 


!******* BEGIN: (==)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION EQ_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ == rhs%x_ad_)
    
    END FUNCTION EQ_DD_D
  
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION EQ_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ == REAL(rhs,DBL_AD))
    
    END FUNCTION EQ_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION EQ_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ == rhs)
    
    END FUNCTION EQ_DR_D

   ! compare a dual with a single
   !-----------------------------------------
    ELEMENTAL FUNCTION EQ_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ == REAL(rhs,DBL_AD))
    
    END FUNCTION EQ_DS_D

    !-----------------------------------------
    ! compare an integer with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION EQ_ID_D(lhs, rhs) RESULT(res)
         INTEGER,INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)==rhs%x_ad_)
    
    END FUNCTION EQ_ID_D

    !-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION EQ_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs==rhs%x_ad_)
    
    END FUNCTION EQ_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION EQ_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)==rhs%x_ad_)
    
    END FUNCTION EQ_SD_D

!******* END: (==)
!--------------------- 


!******* BEGIN: (<=)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION LE_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ <= rhs%x_ad_)
    
    END FUNCTION LE_DD_D
   
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION LE_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ <= REAL(rhs,DBL_AD))
    
    END FUNCTION LE_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION LE_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ <= rhs)
    
    END FUNCTION LE_DR_D

  ! compare a dual with a single
   !----------------------------------------
    ELEMENTAL FUNCTION LE_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ <= REAL(rhs,DBL_AD))
    
    END FUNCTION LE_DS_D
 
    !-----------------------------------------
    ! compare a dual number with an integer
    !----------------------------------------
    ELEMENTAL FUNCTION LE_ID_D(n, rhs) RESULT(res)
         INTEGER,INTENT(IN)::n 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(n,DBL_AD)<=rhs%x_ad_)
    
    END FUNCTION LE_ID_D

!-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION LE_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs<=rhs%x_ad_)
    
    END FUNCTION LE_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION LE_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)<=rhs%x_ad_)
    
    END FUNCTION LE_SD_D

!******* END: (<=)
!--------------------- 

!******* BEGIN: (<)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION LT_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ < rhs%x_ad_)
    
    END FUNCTION LT_DD_D
  
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION LT_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ < REAL(rhs,DBL_AD))
    
    END FUNCTION LT_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION LT_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ < rhs)
    
    END FUNCTION LT_DR_D

   ! compare a dual with a single
   !----------------------------------------
    ELEMENTAL FUNCTION LT_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ < REAL(rhs,DBL_AD))
    
    END FUNCTION LT_DS_D

 !-----------------------------------------
    ! compare a dual number with an integer
    !----------------------------------------
    ELEMENTAL FUNCTION LT_ID_D(n, rhs) RESULT(res)
         INTEGER,INTENT(IN)::n 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(n,DBL_AD)<rhs%x_ad_)
    
    END FUNCTION LT_ID_D

    !-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION LT_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs < rhs%x_ad_)
    
    END FUNCTION LT_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION LT_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD) < rhs%x_ad_)
    
    END FUNCTION LT_SD_D


!******* END: (<)
!--------------------- 


!******* BEGIN: (>=)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION GE_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ >= rhs%x_ad_)
    
    END FUNCTION GE_DD_D
  
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION GE_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ >= REAL(rhs,DBL_AD))
    
    END FUNCTION GE_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION GE_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ >= rhs)
    
    END FUNCTION GE_DR_D

   ! compare a dual with a single
   !----------------------------------------
    ELEMENTAL FUNCTION GE_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ >=REAL(rhs,DBL_AD))
    
    END FUNCTION GE_DS_D
 !-----------------------------------------
    ! compare a dual number with an integer
    !----------------------------------------
    ELEMENTAL FUNCTION GE_ID_D(n, rhs) RESULT(res)
         INTEGER,INTENT(IN)::n 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(n,DBL_AD)>=rhs%x_ad_)
    
    END FUNCTION GE_ID_D

!-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION GE_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs>=rhs%x_ad_)
    
    END FUNCTION GE_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION GE_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)>=rhs%x_ad_)
    
    END FUNCTION GE_SD_D


!******* END: (>=)
!--------------------- 


!******* BEGIN: (>)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION GT_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ > rhs%x_ad_)
    
    END FUNCTION GT_DD_D
  
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION GT_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ > REAL(rhs,DBL_AD))
    
    END FUNCTION GT_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION GT_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ > rhs)
    
    END FUNCTION GT_DR_D

   ! compare a dual with a single
   !----------------------------------------
    ELEMENTAL FUNCTION GT_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ >REAL(rhs,DBL_AD))
    
    END FUNCTION GT_DS_D

  !-----------------------------------------
    ! compare a dual number with an integer
    !----------------------------------------
    ELEMENTAL FUNCTION GT_ID_D(n, rhs) RESULT(res)
         INTEGER,INTENT(IN)::n 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(n,DBL_AD)>rhs%x_ad_)
    
    END FUNCTION GT_ID_D

!-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION GT_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs>rhs%x_ad_)
    
    END FUNCTION GT_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION GT_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)>rhs%x_ad_)
    
    END FUNCTION GT_SD_D

    
!******* END: (>)
!--------------------- 



!******* BEGIN: (/=)
!--------------------- 
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION NE_DD_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs,rhs
         LOGICAL::res
      
         res = (lhs%x_ad_ /= rhs%x_ad_)
    
    END FUNCTION NE_DD_D
  
   ! compare a dual with an integer
   !-----------------------------------------
    ELEMENTAL FUNCTION NE_DI_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         INTEGER,INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ /= REAL(rhs,DBL_AD))
    
    END FUNCTION NE_DI_D

    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    ELEMENTAL FUNCTION NE_DR_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(DBL_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ /= rhs)
    
    END FUNCTION NE_DR_D

   ! compare a dual with a single
   !----------------------------------------
    ELEMENTAL FUNCTION NE_DS_D(lhs, rhs) RESULT(res)
         TYPE (AD_D), INTENT(IN):: lhs
         REAL(SNG_AD),INTENT(IN)::rhs 
         LOGICAL::res 

         res = (lhs%x_ad_ /= REAL(rhs,DBL_AD))
    
    END FUNCTION NE_DS_D

 !-----------------------------------------
    ! compare a dual number with an integer
    !----------------------------------------
    ELEMENTAL FUNCTION NE_ID_D(n, rhs) RESULT(res)
         INTEGER,INTENT(IN)::n 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(n,DBL_AD) /= rhs%x_ad_)
    
    END FUNCTION NE_ID_D

!-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION NE_RD_D(lhs, rhs) RESULT(res)
         REAL(DBL_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (lhs/=rhs%x_ad_)
    
    END FUNCTION NE_RD_D

    !-----------------------------------------
    ! compare a single with a dual
    !----------------------------------------
    ELEMENTAL FUNCTION NE_SD_D(lhs, rhs) RESULT(res)
         REAL(SNG_AD),INTENT(IN)::lhs 
         TYPE (AD_D), INTENT(IN):: rhs
         LOGICAL::res
         
         res = (REAL(lhs,DBL_AD)/=rhs%x_ad_)
    
    END FUNCTION NE_SD_D

!******* END: (/=)
!--------------------- 

    !---------------------------------------------------
    ! ABS of dual numbers
    ! ABS<u,up>=<ABS(u),up SIGN(u)>
    !----------------------------------------------------
    ELEMENTAL FUNCTION ABS_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))
                   
         IF(u%x_ad_>0) THEN
            res%x_ad_ = u%x_ad_
            res%xp_ad_= u%xp_ad_
         ELSE IF (u%x_ad_<0) THEN
            res%x_ad_ = -u%x_ad_
            res%xp_ad_= -u%xp_ad_
         ELSE 
            res%x_ad_ = 0.0D0
!            res%xp_ad_= 0.0 !Set_NaN_D() ! Indicating an undefined derivative, however, for some codes it will cause problem.

            res%xp_ad_ = u%xp_ad_   !> NEED THIS TO ENFORCE ALLOCATION OF res%xp_ad_
            res%xp_ad_ = 0.0D0      !> THEN, SET TO ZERO

         ENDIF

    END FUNCTION ABS_D_D


    !-----------------------------------------
    ! ACOS of dual numbers
    ! ACOS<u,up>=<ACOS(u),-up/sqrt(1-(u%x_ad_)**2)>
    !----------------------------------------
    ELEMENTAL FUNCTION ACOS_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL(DBL_AD)::tmp 
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = ACOS(u%x_ad_)
         IF(u%x_ad_==1.0D0.OR.u%x_ad_==-1.0D0) THEN  
            res%xp_ad_=Set_NaN_D() ! Indicating an undefined derivative
         ELSE 
            tmp= -1.0d0/SQRT(1.0D0-(u%x_ad_)**2)
            res%xp_ad_= u%xp_ad_*tmp
         ENDIF
    END FUNCTION ACOS_D_D


    !-----------------------------------------
    ! ASIN of dual numbers
    ! ASIN<u,up>=<ASIN(u),up 1/SQRT(1-u^2)>
    !----------------------------------------
    ELEMENTAL FUNCTION ASIN_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL (DBL_AD):: tmp
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = ASIN(u%x_ad_) 
         IF(u%x_ad_==1.0D0.OR.u%x_ad_==-1.0D0) THEN  
            res%xp_ad_=Set_NaN_D() ! Indicating an undefined derivative
         ELSE 
            tmp= 1.0d0/SQRT(1.0D0-(u%x_ad_)**2)
            res%xp_ad_= u%xp_ad_*tmp
         ENDIF
    END FUNCTION ASIN_D_D

        
      !-----------------------------------------
    ! COS of dual numbers
    ! COS<u,up>=<COS(u),-up*SIN(u)>
    !----------------------------------------
    ELEMENTAL FUNCTION COS_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL(DBL_AD):: tmp
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = COS(u%x_ad_) 
         tmp=-SIN(u%x_ad_)
         res%xp_ad_= u%xp_ad_*tmp

    END FUNCTION COS_D_D

    !-----------------------------------------
    ! DOT PRODUCT two dual number vectors
    ! <u,up>.<v,vp>=<u.v,u.vp+up.v>
    !----------------------------------------
    FUNCTION DOT_PRODUCT_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u(:),v(:)
         TYPE (AD_D)::res
!         INTEGER:: i
!
         stop "DOT_PRODUCT not implemented"
!         res = AD_D(size(u(1)%xp_ad_))
!         res = 0.
!
!         res%x_ad_ = DOT_PRODUCT(u%x_ad_,v%x_ad_)
!
!         DO i=1,size(u(1)%xp_ad_)
!           res%xp_ad_(i) =DOT_PRODUCT(u%x_ad_,v%xp_ad_(i))+DOT_PRODUCT(u%xp_ad_(i),v%x_ad_)
!         ENDDO
    END FUNCTION DOT_PRODUCT_DD_D

   
   
    !-----------------------------------------
    ! EXPONENTIAL OF dual numbers
    ! EXP<u,up>=<EXP(u),up EXP(u)>
    !----------------------------------------
    ELEMENTAL FUNCTION EXP_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD)::tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         tmp=EXP(u%x_ad_)
         res%x_ad_ = tmp
         res%xp_ad_ =u%xp_ad_*tmp
              
    END FUNCTION EXP_D_D

    !-----------------------------------------
    ! EXPONENTIAL OF dual numbers
    ! EXP<u,up>=<EXP(u),up EXP(u)>
    !----------------------------------------
    ELEMENTAL FUNCTION INT_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD)::tmp
         INTEGER::res

         tmp=u%x_ad_
         res = INT(tmp)
              
    END FUNCTION INT_D_D

    
    !-----------------------------------------
    ! LOG OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    ELEMENTAL FUNCTION LOG_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL(DBL_AD)::tmp
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = LOG(u%x_ad_)
         tmp=1.0d0/u%x_ad_
         res%xp_ad_ =u%xp_ad_*tmp
                       
    END FUNCTION LOG_D_D

!-----------------------------------------
    ! LOG OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    ELEMENTAL FUNCTION LOG10_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL(DBL_AD)::tmp
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = LOG10(u%x_ad_)
         tmp=1.0d0/u%x_ad_/LOG(10.0D0)
         res%xp_ad_ =u%xp_ad_*tmp
                       
    END FUNCTION LOG10_D_D


    !-----------------------------------------
    ! MULTIPLY two dual number matrices
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    FUNCTION MATMUL_DD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u(:,:),v(:,:)
         TYPE (AD_D)::res(SIZE(u,1),SIZE(v,2))
         REAL(DBL_AD) :: xp_ad_u(size(u,1),size(u,2)), xp_ad_v(size(v,1),size(v,2)), xp_ad(size(res,1),size(res,2))
         INTEGER:: i,j,k

         do k = 1,size(v,2)
            do j = 1,size(u,1)
                allocate(res(j,k)%xp_ad_(size(u(1,1)%xp_ad_)))
            end do
        end do


         res%x_ad_ = MATMUL(u%x_ad_,v%x_ad_)
         DO i=1,size(u(1,1)%xp_ad_)
            ! Assemble derivative components for U
            do k = 1,size(u,2)
                do j = 1,size(u,1)
                    xp_ad_u(j,k) = u(j,k)%xp_ad_(i)
                end do
            end do

            ! Assemble derivative components for V
            do k = 1,size(v,2)
                do j = 1,size(v,1)
                    xp_ad_v(j,k) = v(j,k)%xp_ad_(i)
                end do
            end do

            xp_ad = MATMUL(xp_ad_u,v%x_ad_) + MATMUL(u%x_ad_,xp_ad_v)
!            res%xp_ad_(i)= MATMUL(xp_ad_u,v%x_ad_) + MATMUL(u%x_ad_,xp_ad_v)

            do k = 1,size(res,2)
                do j = 1,size(res,1)
                    res(j,k)%xp_ad_(i) = xp_ad(j,k)
                end do
            end do
         ENDDO

!         DO i=1,NDV_AD
!            res%xp_ad_(i)= MATMUL(u%xp_ad_(i),v%x_ad_) + MATMUL(u%x_ad_,v%xp_ad_(i))
!         ENDDO
    END FUNCTION MATMUL_DD_D

    !-----------------------------------------
    ! MULTIPLY a dual number matrix with a dual number 
    ! vector 
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    FUNCTION MATMUL_DV_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u(:,:),v(:)
         TYPE (AD_D)::res(SIZE(u,1))
         REAL(DBL_AD) :: xp_ad_u(size(u,1),size(u,2)), xp_ad_v(size(v)), xp_ad(size(res))
         INTEGER:: i,j,k

         do j = 1,size(u,1)
                allocate(res(j)%xp_ad_(size(u(1,1)%xp_ad_)))
         end do

         res%x_ad_ = MATMUL(u%x_ad_,v%x_ad_)
         DO i=1,size(u(1,1)%xp_ad_)
            ! Assemble derivative components for U
            do k = 1,size(u,2)
                do j = 1,size(u,1)
                    xp_ad_u(j,k) = u(j,k)%xp_ad_(i)
                end do
            end do

            ! Assemble derivative components for V
            do j = 1,size(v)
                    xp_ad_v(j) = v(j)%xp_ad_(i)
            end do


           xp_ad = MATMUL(xp_ad_u,v%x_ad_) + MATMUL(u%x_ad_,xp_ad_v)

           do j = 1,size(res)
               res(j)%xp_ad_(i) = xp_ad(j)
           end do
         ENDDO

!         DO i=1,NDV_AD
!            res%xp_ad_(i)= MATMUL(u%xp_ad_(i),v%x_ad_) + MATMUL(u%x_ad_,v%xp_ad_(i))
!         ENDDO
    END FUNCTION MATMUL_DV_D



    !-----------------------------------------
    ! MULTIPLY Matrix(Real)-Vector(Dual)
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    FUNCTION MATMUL_MV_RD_D(u,v) RESULT(res)
        REAL(DBL_AD), intent(in)     :: u(:,:)
        TYPE(AD_D), intent(in)       :: v(:)
!        TYPE(AD_D)                      :: res(size(v))
!        TYPE(AD_D)                      :: res(size(u,1))  !> Causes memory leak in ifort 15.0.3. Doesn't seem to be deallocating everything correctly
        TYPE(AD_D), allocatable, dimension(:) :: res        !> Declaring as allocatable to fix memory leak
!        REAL(DBL_AD)                    :: xp_ad_v(size(v)), &
!                                           xp_ad(size(u,1))
        REAL(DBL_AD), dimension(size(v),size(v(1)%xp_ad_))      :: xp_ad_vm
        REAL(DBL_AD), dimension(size(u,1),size(v(1)%xp_ad_))    :: res_xp_m
        INTEGER:: i,j



        allocate(res(size(u,1)))
        do j = 1,size(res)
               allocate(res(j)%xp_ad_(size(v(1)%xp_ad_)))
        end do

        !
        ! Standard matrix multiplication of function values
        !
        res%x_ad_ = MATMUL(u,v%x_ad_)

        !
        ! Assemble derivative components as a matrix
        !
        do i = 1,size(v)
            xp_ad_vm(i,:) = v(i)%xp_ad_
        end do


        !print*, size(u,1), size(u,2), size(xp_ad_vm,1), size(xp_ad_vm,2)
        res_xp_m = matmul(u,xp_ad_vm)



        !
        ! Distribute derivatives
        !
        do i = 1,size(res)
            res(i)%xp_ad_ = res_xp_m(i,:)
        end do




!        do i=1,size(v(1)%xp_ad_)
!
!            ! Assemble derivative components for V
!            do j = 1,size(v)
!                xp_ad_v(j) = v(j)%xp_ad_(i)
!            end do
!
!            xp_ad = MATMUL(u,xp_ad_v)
!
!            do j = 1,size(res)
!                res(j)%xp_ad_(i) = xp_ad(j)
!            end do
!
!       end do

    END FUNCTION MATMUL_MV_RD_D




    !-----------------------------------------
    ! MULTIPLY a dual vector with a  dual matrix
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    FUNCTION MATMUL_VD_D(u,v) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u(:),v(:,:)
         TYPE (AD_D)::res(SIZE(v,2))
         REAL(DBL_AD) :: xp_ad_u(size(u)), xp_ad_v(size(v,1),size(v,2)), xp_ad(size(res))
         INTEGER::i,j,k

        ! Disable vector-matrix multiplication
        stop "Vector-Matrix multiplication disabled"

         do j = 1,size(v,2)
                allocate(res(j)%xp_ad_(size(u(1)%xp_ad_)))
         end do

         res%x_ad_= MATMUL(u%x_ad_,v%x_ad_)
         DO i=1,size(u(1)%xp_ad_)
            ! Assemble derivative components for U
            do j = 1,size(u)
                xp_ad_u(j) = u(j)%xp_ad_(i)
            end do

            ! Assemble derivative components for V
            do k = 1,size(v,2)
                do j = 1,size(v,1)
                    xp_ad_v(j,k) = v(j,k)%xp_ad_(i)
                end do
            end do

            xp_ad = MATMUL(xp_ad_u,v%x_ad_) + MATMUL(u%x_ad_,xp_ad_v)

            do j = 1,size(res)
                res(j)%xp_ad_(i) = xp_ad(j)
            end do
!            res%xp_ad_(i)= MATMUL(xp_ad_u,v%x_ad_) + MATMUL(u%x_ad_,xp_ad_v)
         ENDDO

!         DO i=1,size(u(1)%xp_ad_)
!            res%xp_ad_(i)= MATMUL(u%xp_ad_(i),v%x_ad_) + MATMUL(u%x_ad_,v%xp_ad_(i))
!         ENDDO
    END FUNCTION MATMUL_VD_D

    !-----------------------------------------
    ! Obtain the max of 2 to 5 dual numbers
    !----------------------------------------    
    ELEMENTAL FUNCTION MAX_DD_D(val1, val2, val3, val4,val5) RESULT(res)
        TYPE (AD_D), INTENT(IN)::val1, val2
        TYPE (AD_D), INTENT(IN),OPTIONAL:: val3, val4,val5
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(val1%xp_ad_)))
   
        IF (val1%x_ad_ > val2%x_ad_) THEN
            res = val1
        ELSE 
            res = val2
        ENDIF
        IF(PRESENT(val3))THEN
           IF(res%x_ad_<val3%x_ad_) res=val3
        ENDIF
        IF(PRESENT(val4))THEN
           IF(res%x_ad_<val4%x_ad_) res=val4
        ENDIF
        IF(PRESENT(val5))THEN
           IF(res%x_ad_<val5%x_ad_) res=val5
        ENDIF

    END FUNCTION MAX_DD_D


    !-----------------------------------------
    ! Obtain the max of a dual number and an integer
    !----------------------------------------    
    ELEMENTAL FUNCTION MAX_DI_D(u, n) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u
        INTEGER,INTENT(IN)::n
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(u%xp_ad_)))
   
        IF (u%x_ad_ > n) THEN
            res = u
        ELSE 
            res = n
        ENDIF
        
    END FUNCTION MAX_DI_D

    !-----------------------------------------
    ! Obtain the max of a dual number and a real number
    !----------------------------------------    
    ELEMENTAL FUNCTION MAX_DR_D(u, n) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u
        REAL(DBL_AD),INTENT(IN)::n
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(u%xp_ad_)))
   
        IF (u%x_ad_ > n) THEN
            res = u
        ELSE 
            res = n
        ENDIF
        
    END FUNCTION MAX_DR_D


 !-----------------------------------------
    ! Obtain the max of a dual number and a single
    !----------------------------------------    
    ELEMENTAL FUNCTION MAX_DS_D(u, n) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u
        REAL(SNG_AD),INTENT(IN)::n
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(u%xp_ad_)))
   
        IF (u%x_ad_ > n) THEN
            res = u
        ELSE 
            res = n
        ENDIF
        
    END FUNCTION MAX_DS_D

   
    !---------------------------------------------------
    ! Obtain the max of a real and a dual
    ! note the real argument is renamed as r to avoid
    ! the ambiguity due to positional association using
    ! keywords when functions/subroutines are called
    !---------------------------------------------------    
     ELEMENTAL FUNCTION MAX_RD_D(r,u) RESULT(res)
        REAL(DBL_AD),INTENT(IN)::r
        TYPE (AD_D), INTENT(IN)::u
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(u%xp_ad_)))
   
        IF (u%x_ad_ > r) THEN
            res = u
        ELSE 
            res = r
        ENDIF
        
    END FUNCTION MAX_RD_D

    !-----------------------------------------
    ! Obtain the max value of vector u
    !----------------------------------------    
    FUNCTION MAXVAL_D_D(u) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u(:)
        INTEGER::iloc(1)
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(u(1)%xp_ad_)))

        stop "MAXVAL - might not work properly"

        iloc=MAXLOC(u%x_ad_)
        res=u(iloc(1))
        
    END FUNCTION MAXVAL_D_D

   
    !-----------------------------------------
    ! Obtain the min of 2 to 4 dual numbers
    !----------------------------------------    
    ELEMENTAL FUNCTION MIN_DD_D(val1, val2, val3, val4) RESULT(res)
        TYPE (AD_D), INTENT(IN)::val1, val2
        TYPE (AD_D), INTENT(IN),OPTIONAL:: val3, val4
        TYPE (AD_D)::res
   
        IF (val1%x_ad_ < val2%x_ad_) THEN
            res = val1
        ELSE 
            res = val2
        ENDIF
        IF(PRESENT(val3))THEN
           IF(res%x_ad_>val3%x_ad_) res=val3
        ENDIF
        IF(PRESENT(val4))THEN
           IF(res%x_ad_>val4%x_ad_) res=val4
        ENDIF

    END FUNCTION MIN_DD_D

    !-----------------------------------------
    ! Obtain the min of a dual number and a single
    !----------------------------------------    
    ELEMENTAL FUNCTION MIN_DR_D(u, n) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u
        REAL(DBL_AD),INTENT(IN)::n
        TYPE (AD_D)::res
   
        IF (u%x_ad_ < n) THEN
            res = u
        ELSE 
            res = n
        ENDIF
        
    END FUNCTION MIN_DR_D

    !-----------------------------------------
    ! Obtain the min of a dual number and a single
    !----------------------------------------    
    ELEMENTAL FUNCTION MIN_DS_D(u, n) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u
        REAL(SNG_AD),INTENT(IN)::n
        TYPE (AD_D)::res

        IF (u%x_ad_ < n) THEN
            res = u
        ELSE 
            res = n
        ENDIF
        
    END FUNCTION MIN_DS_D

    !-----------------------------------------
    ! Obtain the min value of vector u
    !----------------------------------------    
    FUNCTION MINVAL_D_D(u) RESULT(res)
        TYPE (AD_D), INTENT(IN)::u(:)
        INTEGER::iloc(1)
        TYPE (AD_D)::res

        stop "MINVAL - might not work properly"

        iloc=MINLOC(u%x_ad_)
        res=u(iloc(1))
        
    END FUNCTION MINVAL_D_D

    !------------------------------------------------------
    !Returns the nearest integer to u%x, ELEMENTAL
    !------------------------------------------------------
    ELEMENTAL FUNCTION NINT_D_D(u) RESULT(res)
        TYPE (AD_D), INTENT(IN):: u
        INTEGER::res
        
        res=NINT(u%x_ad_)
  
    END FUNCTION NINT_D_D



    !----------------------------------------------------------------
    ! SIGN(a,b) with two dual numbers as inputs, 
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    ELEMENTAL FUNCTION SIGN_DD_D(val1, val2) RESULT(res)
        TYPE (AD_D), INTENT(IN) :: val1, val2
        TYPE (AD_D)::res

   
        IF (val2%x_ad_ < 0.D0) THEN
            res = -ABS(val1)
        ELSE
            res =  ABS(val1)
        ENDIF
        
     END FUNCTION SIGN_DD_D
  

    !----------------------------------------------------------------
    ! SIGN(a,b) with one real and one dual number as inputs, 
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    ELEMENTAL FUNCTION SIGN_RD_D(val1, val2) RESULT(res)
        REAL(DBL_AD),INTENT(IN)::val1
        TYPE (AD_D), INTENT(IN) :: val2
        TYPE (AD_D)::res
        allocate(res%xp_ad_(size(val2%xp_ad_)))
   
        IF (val2%x_ad_ < 0.D0) THEN
            res = -ABS(val1)
        ELSE
            res = ABS(val1)
        ENDIF
        
     END FUNCTION SIGN_RD_D
  
  
    !-----------------------------------------
    ! SIN of dual numbers
    ! SIN<u,up>=<SIN(u),up COS(u)>
    !----------------------------------------
    ELEMENTAL FUNCTION SIN_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         TYPE (AD_D)::res
         REAL (DBL_AD):: tmp
         allocate(res%xp_ad_(size(u%xp_ad_)))

         res%x_ad_ = SIN(u%x_ad_) 
         tmp=COS(u%x_ad_)
         res%xp_ad_= u%xp_ad_*tmp

    END FUNCTION SIN_D_D


    !-----------------------------------------
    ! SQRT of dual numbers
    ! SQRT<u,up>=<SQRT(u),up/2/sqrt(u) >
    !----------------------------------------
    !IMPURE ELEMENTAL FUNCTION SQRT_D_D(u) RESULT(res)
    ELEMENTAL FUNCTION SQRT_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u
         REAL(DBL_AD):: tmp
         TYPE (AD_D)::res
         allocate(res%xp_ad_(size(u%xp_ad_)))

         tmp=SQRT(u%x_ad_) 
  
         res%x_ad_ = tmp
         IF(tmp==0.0d0) THEN
            !print*, 'hi'
            res%xp_ad_ = Set_NaN_D()
            !res%xp_ad_ = 0.0d0
         ELSE
            tmp= 0.5D0/tmp
            res%xp_ad_ = u%xp_ad_*tmp
         ENDIF
    END FUNCTION SQRT_D_D


 
    !-----------------------------------------
    ! SUM OF A DUAL ARRAY
    !----------------------------------------
    FUNCTION SUM_D_D(u) RESULT(res)
         TYPE (AD_D), INTENT(IN)::u(:)
         TYPE (AD_D)::res
         INTEGER:: i,j
         REAL(DBL_AD):: tmp
!         allocate(res%xp_ad_(size(u(1)%xp_ad_)))

         res%x_ad_ = SUM(u%x_ad_)
         DO i=1,size(u(1)%xp_ad_)
         tmp = REAL(0.0,DBL_AD)
         
            DO j=1,size(u)
                tmp = tmp + u(j)%xp_ad_(i)
            ENDDO
            res%xp_ad_(i) = tmp
         ENDDO
    END FUNCTION SUM_D_D

    ELEMENTAL FUNCTION Set_NaN_D() RESULT(res)
         REAL(DBL_AD)::res

         res=SQRT(negative_one) 

    END FUNCTION Set_NaN_D
 







!    subroutine destructor(self)
!        type(AD_D), intent(inout) :: self
!
!
!    end subroutine








END MODULE  DNAD_D


