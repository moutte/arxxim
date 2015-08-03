MODULE M_Safe_Functions

  !=========================================
  ! Safe Functions Exp, Log, Log10
  ! Truncated to avoid Out of Range Errors
  !=========================================

  IMPLICIT NONE
  PRIVATE

  REAL(kind=8), PARAMETER:: Ln10= 2.30258509299404D0
  REAL(kind=8), PARAMETER:: Zero = 0.D0
  REAL(kind=8), PARAMETER:: One = 1.D0

  REAL(kind=8) :: MinExpDP= -300.D0
  REAL(kind=8) :: MaxExpDP= +300.D0
  REAL(kind=8) :: TinyDP= 1.D-300

  !// Set Parameters
  PUBLIC :: Set_MinExpDP
  PUBLIC :: Set_MaxExpDP
  PUBLIC :: Set_TinyDP

  !// Scalar values
  PUBLIC :: Safe_Exp
  PUBLIC :: Safe_Log
  PUBLIC :: Safe_Log10

  !// Function Scalar values
  PUBLIC :: FSafe_Exp
  PUBLIC :: FSafe_Log
  PUBLIC :: FSafe_Log10
  
  !// Vector values
  PUBLIC :: Safe_vExp
  PUBLIC :: Safe_vLog
  PUBLIC :: Safe_vLog10

  !// Function Vector values
  PUBLIC :: FSafe_vExp
  PUBLIC :: FSafe_vLog
  PUBLIC :: FSafe_vLog10

CONTAINS

  !==================== PARAMETER SETTING =====================

  SUBROUTINE Set_MinExpDP(value)
    IMPLICIT NONE
    REAL(kind=8) :: value
    MinExpDP= value
  END SUBROUTINE Set_MinExpDP

  !---

  SUBROUTINE Set_MaxExpDP(value)
    IMPLICIT NONE
    REAL(kind=8) :: value
    MaxExpDP= value
  END SUBROUTINE Set_MaxExpDP

  !---
  
  SUBROUTINE Set_TinyDP(value)
    IMPLICIT NONE
    REAL(kind=8) :: value
    TinyDP= value
  END SUBROUTINE Set_TinyDP

  !====================>>  SCALAR VALUES <<=====================

  SUBROUTINE Safe_Exp(X,Expx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8), intent(out) :: Expx
    !--
    Expx = 0.; 
    IF( (X>MinExpDP) .AND. (X<MaxExpDP) ) THEN
       Expx= EXP(X)
    ELSE
       IF(X<=MinExpDP) Expx= EXP(MinExpDP)
       IF(X>=MaxExpDP) Expx= EXP(MaxExpDP)
    END IF

  END SUBROUTINE Safe_Exp

  !---
  
  FUNCTION FSafe_Exp(X) result (Expx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8) :: Expx
    !--
    Expx = 0.; 
    IF( (X>MinExpDP) .AND. (X<MaxExpDP) ) THEN
       Expx= EXP(X)
    ELSE
       IF(X<=MinExpDP) Expx= EXP(MinExpDP)
       IF(X>=MaxExpDP) Expx= EXP(MaxExpDP)
    END IF
  END FUNCTION FSafe_Exp

  !---

  SUBROUTINE Safe_Log(X,Logx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8), intent(out) :: Logx
    !--
    IF( X>TinyDP ) THEN
       Logx= log(X)
    ELSE
       IF ( X>= Zero ) Logx= LOG(TinyDP)
       IF ( X< Zero )  Logx= LOG(X) ! Error Log Undefined !    
    END IF

  END SUBROUTINE Safe_Log

  !---
  
  FUNCTION FSafe_Log(X) result(LogX) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8) :: Logx
    !--
    IF( X>TinyDP ) THEN
      Logx= log(X)
    ELSE
      IF ( X>= Zero ) Logx= LOG(TinyDP)
      IF ( X< Zero )  Logx= LOG(X) ! Error Log Undefined !    
    END IF
    
  END FUNCTION FSafe_Log

  !---

  SUBROUTINE Safe_Log10(X,Logx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8), intent(out) :: Logx
    !--
    CALL Safe_Log(X,Logx) 
    Logx = LogX + Ln10;

  END SUBROUTINE Safe_Log10

  !---
  
   FUNCTION FSafe_Log10(X) result(Logx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: X
    REAL(kind=8) :: Logx
    !--
    CALL Safe_Log(X,Logx) 
    Logx = LogX + Ln10;

  END FUNCTION FSafe_Log10

  !====================>>  VECTOR VALUES <<=====================


  SUBROUTINE Safe_vExp(vX,vExpx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8), intent(out) :: vExpx(:)
    !--
    vExpx = 0.
    WHERE((vX>MinExpDP) .AND. (vX<MaxExpDP))
      vExpx=EXP(vX)
    ELSEWHERE
      WHERE(vX<=MinExpDP) vExpx= EXP(MinExpDP)
      WHERE(vX>=MaxExpDP) vExpx= EXP(MaxExpDP)
    END WHERE

  END SUBROUTINE Safe_vExp

  !---

  FUNCTION FSafe_vExp(vX) result(vExpx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8) :: vExpx(size(vX))
    !--
    vExpx = 0.
    WHERE((vX>MinExpDP) .AND. (vX<MaxExpDP))
      vExpx=EXP(vX)
    ELSEWHERE
      WHERE(vX<=MinExpDP) vExpx= EXP(MinExpDP)
      WHERE(vX>=MaxExpDP) vExpx= EXP(MaxExpDP)
    END WHERE
  END FUNCTION FSafe_vExp
    
  !---
  
  SUBROUTINE Safe_vLog(vX,vLogx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8), intent(out) :: vLogx(:)
    !
    WHERE( vX> TinyDP )
       vLogx= LOG(vX)
    ELSEWHERE
       WHERE(vX>=0.) vLogx= LOG(TinyDP)
       WHERE(vX<0.)  vLogx= LOG(vX)     ! Error Log Undefined !    
    END WHERE

  END SUBROUTINE Safe_vLog

  !---

  FUNCTION FSafe_vLog(vX) result (vLogx)
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8) :: vLogx(size(vX))
    !
    WHERE( vX> TinyDP )
       vLogx= LOG(vX)
    ELSEWHERE
       WHERE(vX>=0.) vLogx= LOG(TinyDP)
       WHERE(vX<0.)  vLogx= LOG(vX)     ! Error Log Undefined !    
    END WHERE

  END FUNCTION FSafe_vLog

  !---

  SUBROUTINE Safe_vLog10(vX,vLogx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8), intent(out) :: vLogx(:)
    !--
    CALL Safe_vLog(vX,vLogx) 
    vLogx = vLogx + Ln10

  END SUBROUTINE Safe_vLog10

  !---

  FUNCTION FSafe_vLog10(vX) result (vLogx) 
    IMPLICIT NONE
    REAL(kind=8), intent(in)  :: vX(:)
    REAL(kind=8)  :: vLogx(size(vX))
    !--
    CALL Safe_vLog(vX,vLogx) 
    vLogx = vLogx + Ln10

  END FUNCTION FSafe_vLog10

  !---

END MODULE M_Safe_Functions
