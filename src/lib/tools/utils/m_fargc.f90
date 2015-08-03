MODULE M_FArgC
  !=============================================================== 
  ! Fortran Command Arguments Functions
  ! Wrapper to System Functions : IARGC and GETARG
  !
  ! These Functions are not compatible with standard F90 
  ! Please Modify this file if you encounter difficulties
  ! concerning these functions with your compiler
  !
  !===============================================================
  IMPLICIT NONE
  PRIVATE
  !---
  PUBLIC :: F_IARGC
  PUBLIC :: F_GETARG

CONTAINS

  FUNCTION F_IARGC()
  
    INTEGER :: F_IARGC
    INTEGER IARGC
    
    F_IARGC = 0
    F_IARGC = IARGC() ! Comment this line if problems

  END FUNCTION F_IARGC

  !---
  
  SUBROUTINE F_GETARG(I, STR)
  
    INTEGER, INTENT(IN) :: I
    CHARACTER(LEN=*), INTENT(OUT):: STR
    
    STR = ""
    CALL GETARG(I, STR)       ! Comment this line if problems
    
  END SUBROUTINE F_GETARG


END MODULE M_FArgC
