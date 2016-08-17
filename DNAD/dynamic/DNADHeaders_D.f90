!******************************************************************************
!* Dual Number Automatic Differentiation (DNAD) of Fortran Codes
!* -------------------------------------------------------------
!* COPYRIGHT (c) Wenbin Yu, All rights reserved, you are free to copy, 
!* modify or translate this code to other languages such as c/c++. If 
!* you find a bug please let me know through wenbinyu.heaven@gmail.com. If 
!* you added new functions and want to share with others, please let me know too.
!* You are welcome to share your successful stories with us through 
!* http://groups.google.com/group/hifi-comp. 
!******************************************************************************
!* Simple Instruction:
!*---------------------
!* This is the header file: define the interface needed for overloading intrinsic 
!* functions and operators.  This file should put in the same folder as dnad.f90. 
!* If a function or operation is not defined (unlikely), you need to create a 
!* corresponding interface in this file (DNADHeaders.f90) and a corresponding 
!* function/subroutine in DNAD.f90.
!*********************************************************************************
!* Acknowledgements
!----------------------
!* The development of DNAD is supported, in part, by the Chief Scientist 
!* Innovative Research Fund at AFRL/RB WPAFB, and by Department of Army 
!* SBIR (Topic A08-022) through Advanced Dynamics Inc. The views and 
!* conclusions contained herein are those of the authors and should not be 
!* interpreted as necessarily representing the official policies or endorsement,
!* either expressed or implied, of the funding agency.    
!*********************************************************************************

!******** Interfaces for operator overloading
    PUBLIC ASSIGNMENT (=)
    INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE ASSIGN_DI_D  ! dual=integer, ELEMENTAL
        MODULE PROCEDURE ASSIGN_DR_D  ! dual=real, ELEMENTAL
        MODULE PROCEDURE ASSIGN_DS_D  ! dual=integer, ELEMENTAL
        MODULE PROCEDURE ASSIGN_ID_D  ! integer=dual, ELEMENTAL

!It is found out that compilers will know how to assign a scalar to vectors and matrices 
!if the assignment is defined for assigning a scalar to a dual number. Hence it
!is unnecessary to define such assignment overloadings. 
                 
    END INTERFACE

    
    PUBLIC OPERATOR (+)
    INTERFACE OPERATOR (+)    
        MODULE PROCEDURE ADD_D_D   ! +dual number, ELEMENTAL
        MODULE PROCEDURE ADD_DD_D  ! dual+ dual, ELEMENTAL
        MODULE PROCEDURE ADD_DI_D  ! dual+ integer, ELEMENTAL
        MODULE PROCEDURE ADD_DR_D  ! dual+ real, ELEMENTAL
        MODULE PROCEDURE ADD_DS_D  ! dual+ single, ELEMENTAL
        MODULE PROCEDURE ADD_ID_D  ! integer+dual, ELEMENTAL
        MODULE PROCEDURE ADD_RD_D  ! real+ dual, ELEMENTAL
        MODULE PROCEDURE ADD_SD_D  ! single+dual, ELEMENTAL
!It is found out that these overloads also cover the cases when one of the operand is a matrix or vector.
!Of course, if both operands are vectors or matrices, they should be of the same shape 
    END INTERFACE

    PUBLIC OPERATOR (-)
    INTERFACE OPERATOR (-)
        MODULE PROCEDURE MINUS_D_D   ! negate a dual number,ELEMENTAL
        MODULE PROCEDURE MINUS_DD_D  ! dual -dual,ELEMENTAL
        MODULE PROCEDURE MINUS_DI_D  ! dual-integer,ELEMENTAL
        MODULE PROCEDURE MINUS_DR_D  ! dual-real,ELEMENTAL
        MODULE PROCEDURE MINUS_DS_D  ! dual-single,ELEMENTAL
        MODULE PROCEDURE MINUS_ID_D  ! integer-dual,ELEMENTAL
        MODULE PROCEDURE MINUS_RD_D  ! real-dual,ELEMENTAL
        MODULE PROCEDURE MINUS_SD_D  ! single-dual,ELEMENTAL
!It is found out that these overloads also cover the cases when one of the operand is a matrix or vector.
!Of course, if both operands are vectors or matrices, they should be of the same shape 
    END INTERFACE

    PUBLIC OPERATOR (*)
    INTERFACE OPERATOR (*)
        MODULE PROCEDURE MULT_DD_D    ! dual*dual, ELEMENTAL
        MODULE PROCEDURE MULT_DI_D    ! dual*integer,ELEMENTAL
        MODULE PROCEDURE MULT_DR_D    ! dual*real,ELEMENTAL
        MODULE PROCEDURE MULT_DS_D    ! dual*single,ELEMENTAL
        MODULE PROCEDURE MULT_ID_D    ! integer*dual,ELEMENTAL
        MODULE PROCEDURE MULT_RD_D    ! real*dual,ELEMENTAL
        MODULE PROCEDURE MULT_SD_D    ! single*dual,ELEMENTAL
!It is found out that these overloads also cover the cases when one of the operand is a matrix or vector.
!Of course, if both operands are vectors or matrices, they should be of the same shape 
    END INTERFACE

    PUBLIC OPERATOR (/)
    INTERFACE OPERATOR (/)
        MODULE PROCEDURE DIV_DD_D ! dual/dual,ELEMENTAL
        MODULE PROCEDURE DIV_DI_D ! dual/integer, ELEMENTAL
        MODULE PROCEDURE DIV_DR_D ! dual/real,EMENTAL
        MODULE PROCEDURE DIV_DS_D ! dual/single,EMENTAL
        MODULE PROCEDURE DIV_ID_D ! integer/dual, ELEMENTAL
        MODULE PROCEDURE DIV_RD_D ! real/dual, ELEMENTAL
        MODULE PROCEDURE DIV_SD_D ! single/dual, ELEMENTAL
    END INTERFACE

    PUBLIC OPERATOR (**)
    INTERFACE OPERATOR (**)
        MODULE PROCEDURE POW_I_D ! power a dual number to an integer power,ELEMENTAL
        MODULE PROCEDURE POW_R_D ! power a dual number to a real (double precision) power, ELEMENTAL
        MODULE PROCEDURE POW_S_D ! power a dual number to a real (single precision) power, ELEMENTAL
        MODULE PROCEDURE POW_D_D ! power a dual number to a dual power, ELEMENTAL
    END INTERFACE
   
    PUBLIC OPERATOR (==)
    INTERFACE OPERATOR (==)
        MODULE PROCEDURE EQ_DD_D ! compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE EQ_DI_D ! compare a dual and an integer, ELEMENTAL
        MODULE PROCEDURE EQ_DR_D ! compare a dual and a real, ELEMENTAL
        MODULE PROCEDURE EQ_DS_D ! compare a dual and a single, ELEMENTAL
        MODULE PROCEDURE EQ_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE EQ_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE EQ_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE
    
    PUBLIC OPERATOR (<=)
    INTERFACE OPERATOR (<=)
        MODULE PROCEDURE LE_DD_D  ! compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE LE_DI_D  ! compare a dual and an integer, ELEMENTAL
        MODULE PROCEDURE LE_DR_D  ! compare a dual and a real,ELEMENTAL
        MODULE PROCEDURE LE_DS_D  ! compare a dual and a single,ELEMENTAL
        MODULE PROCEDURE LE_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE LE_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE LE_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE
    
    PUBLIC OPERATOR (<)
    INTERFACE OPERATOR (<)
        MODULE PROCEDURE LT_DD_D  !compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE LT_DI_D  !compare a dual and an integer, ELEMENTAL
        MODULE PROCEDURE LT_DR_D  !compare dual with a real, ELEMENTAL
        MODULE PROCEDURE LT_DS_D ! compare <u,up> and a single
        MODULE PROCEDURE LT_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE LT_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE LT_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE

    PUBLIC OPERATOR (>=)
    INTERFACE OPERATOR (>=)
        MODULE PROCEDURE GE_DD_D ! compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE GE_DI_D ! compare dual with integer, ELEMENTAL
        MODULE PROCEDURE GE_DR_D ! compare a dual number with a real number, ELEMENTAL
        MODULE PROCEDURE GE_DS_D ! compare dual with a single, ELEMENTAL
        MODULE PROCEDURE GE_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE GE_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE GE_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE

    PUBLIC OPERATOR (>)
    INTERFACE OPERATOR (>)
        MODULE PROCEDURE GT_DD_D  !compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE GT_DI_D  !compare a dual and an integer, ELEMENTAL
        MODULE PROCEDURE GT_DR_D  !compare dual with a real, ELEMENTAL
        MODULE PROCEDURE GT_DS_D ! compare <u,up> and a single
        MODULE PROCEDURE GT_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE GT_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE GT_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE
  
    PUBLIC OPERATOR (/=)
    INTERFACE OPERATOR (/=)
        MODULE PROCEDURE NE_DD_D  !compare two dual numbers, ELEMENTAL
        MODULE PROCEDURE NE_DI_D  !compare a dual and an integer, ELEMENTAL
        MODULE PROCEDURE NE_DR_D  !compare dual with a real, ELEMENTAL
        MODULE PROCEDURE NE_DS_D ! compare <u,up> and a single
        MODULE PROCEDURE NE_ID_D ! compare an integer with a dual number, ELEMENTAL
        MODULE PROCEDURE NE_RD_D ! compare a real with a dual number, ELEMENTAL
        MODULE PROCEDURE NE_SD_D ! compare a single with a dual number, ELEMENTAL
    END INTERFACE
    
    
!------------------------------------------------
! Interfaces for intrinsic functions overloading
!------------------------------------------------
   PUBLIC ABS
   INTERFACE ABS
        MODULE PROCEDURE ABS_D_D  ! obtain the absolute value of a dual number, ELEMENTAL
   END INTERFACE
 
   PUBLIC DABS
   INTERFACE DABS
        MODULE PROCEDURE ABS_D_D ! the same as ABS, used for some old fortran commands
   END INTERFACE

   PUBLIC ACOS
   INTERFACE ACOS
        MODULE PROCEDURE ACOS_D_D ! obtain the arccosine of a dual number, ELEMENTAL
   END INTERFACE
 
   PUBLIC ASIN
   INTERFACE ASIN
        MODULE PROCEDURE ASIN_D_D ! obtain the arcsine of a dual number, ELEMENTAL
   END INTERFACE


!   ! Added by Nathan A. Wukie - 2/23/2016
!   PUBLIC ATAN
!   INTERFACE ATAN
!        MODULE PROCEDURE ATAN_D_D ! obtain the arcsine of a dual number, ELEMENTAL
!   END INTERFACE

 
   PUBLIC COS
   INTERFACE COS
        MODULE PROCEDURE COS_D_D ! obtain the cosine of a dual number, ELEMENTAL
   END INTERFACE
 
   PUBLIC DCOS
   INTERFACE DCOS
        MODULE PROCEDURE COS_D_D ! obtain the cosine of a dual number, ELEMENTAL
   END INTERFACE
  
   
!   PUBLIC DOT_PRODUCT
!   INTERFACE DOT_PRODUCT
!        MODULE PROCEDURE DOT_PRODUCT_DD_D ! dot product two dual number vectors
!   END INTERFACE

   PUBLIC EXP 
   INTERFACE EXP
        MODULE PROCEDURE EXP_D_D ! obtain the exponential of a dual number, ELEMENTAL
   END INTERFACE


   PUBLIC INT 
   INTERFACE INT
        MODULE PROCEDURE INT_D_D ! obtain the integer part of a dual number, ELEMENTAL
   END INTERFACE

   PUBLIC LOG 
   INTERFACE LOG
        MODULE PROCEDURE LOG_D_D ! obtain the log of a dual number, ELEMENTAL
   END INTERFACE
   
   PUBLIC LOG10 
   INTERFACE LOG10
        MODULE PROCEDURE LOG10_D_D ! obtain the log of a dual number, ELEMENTAL
   END INTERFACE
 
   PUBLIC MATMUL
   INTERFACE MATMUL
        MODULE PROCEDURE MATMUL_DD_D     ! matrix multiplies of two dual matrices
        MODULE PROCEDURE MATMUL_DV_D     ! matrix multiplies of a dual matrix with a dual vector
        MODULE PROCEDURE MATMUL_VD_D     ! matrix multiplies of a dual vector with a dual matrix
        MODULE PROCEDURE MATMUL_MV_RD_D  ! matrix-vector :: real-dual
   END INTERFACE
   

   PUBLIC MAX
   INTERFACE MAX
        MODULE PROCEDURE MAX_DD_D ! obtain the max of from two to four dual numbers, ELEMENTAL
       MODULE PROCEDURE MAX_DI_D ! obtain the max of from a dual number and an integer, ELEMENTAL
       MODULE PROCEDURE MAX_DR_D ! obtain the max of from a dual number and a real, ELEMENTAL
       MODULE PROCEDURE MAX_DS_D ! obtain the max of from a dual number and a real in single precision, ELEMENTAL
       MODULE PROCEDURE MAX_RD_D ! obtain the max of from a real,and a dual number,  ELEMENTAL
   END INTERFACE

   PUBLIC DMAX1
   INTERFACE DMAX1
        MODULE PROCEDURE MAX_DD_D ! obtain the max of from two to four dual numbers, ELEMENTAL
   END INTERFACE
   
   PUBLIC MAXVAL
   INTERFACE MAXVAL
        MODULE PROCEDURE MAXVAL_D_D ! obtain the maxval  of a dual number vectgor
   END INTERFACE
   
   PUBLIC MIN
   INTERFACE MIN
        MODULE PROCEDURE MIN_DD_D ! obtain the min of from two to four dual numbers, ELEMENTAL
         MODULE PROCEDURE MIN_DR_D ! obtain the min of a dual and a real, ELEMENTAL
         MODULE PROCEDURE MIN_DS_D ! obtain the min of a dual and a single, ELEMENTAL
   END INTERFACE

   PUBLIC DMIN1
   INTERFACE DMIN1
        MODULE PROCEDURE MIN_DD_D ! obtain the min of from two to four dual numbers, ELEMENTAL
   END INTERFACE
 
   PUBLIC MINVAL
   INTERFACE MINVAL
        MODULE PROCEDURE MINVAL_D_D ! obtain the maxval  of a dual number vectgor
   END INTERFACE
  
   PUBLIC NINT
   INTERFACE NINT
        MODULE PROCEDURE NINT_D_D ! Returns the nearest integer to the argument, ELEMENTAL
   END INTERFACE

   PUBLIC SIGN   
   INTERFACE  SIGN
     MODULE PROCEDURE  SIGN_DD_D ! SIGN(a,b) with two dual numbers as inputs, the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
     MODULE PROCEDURE  SIGN_RD_D ! SIGN(a,b) with a real and a dual, the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
   END INTERFACE

   PUBLIC SIN  
   INTERFACE SIN
        MODULE PROCEDURE SIN_D_D ! obtain sine of a dual number, ELEMENTAL
   END INTERFACE

   PUBLIC DSIN  
   INTERFACE DSIN
        MODULE PROCEDURE SIN_D_D ! obtain sine of a dual number, ELEMENTAL
   END INTERFACE

   PUBLIC SQRT  
   INTERFACE SQRT
        MODULE PROCEDURE SQRT_D_D ! obtain the sqrt of a dual number, ELEMENTAL
   END INTERFACE

   PUBLIC SUM  
   INTERFACE SUM
        MODULE PROCEDURE SUM_D_D ! sum a dual array
   END INTERFACE

