MODULE M_Path_Vars !variables for path calculations
  USE M_Kinds
  USE M_T_Tpcond,ONLY: T_TPCond
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  LOGICAL, ALLOCATABLE:: vLPath(:),vFasFound(:)
  INTEGER, ALLOCATABLE:: vPhasBegin(:),vPhasFinal(:) !index of phases in vMixFas
  REAL(dp),ALLOCATABLE:: tPathData(:,:)
  REAL(dp),ALLOCATABLE:: vPathLogK(:)
  !
  INTEGER:: DimPath
  INTEGER:: iLogK
  !
  INTEGER,PARAMETER:: TotalMixStep= 100
  !
  REAL(dp),ALLOCATABLE:: tGrt(:,:)
  REAL(dp),ALLOCATABLE:: tPathResults(:,:)
  !
  !CHARACTER(LEN=15):: PathMode
  !
  TYPE(T_TPCond),ALLOCATABLE:: vTPpath(:)
  !
CONTAINS

SUBROUTINE Path_Vars_Clean
  IF(ALLOCATED(vLPath))      DEALLOCATE(vLPath)
  IF(ALLOCATED(vFasFound))   DEALLOCATE(vFasFound)
  IF(ALLOCATED(vPhasBegin))  DEALLOCATE(vPhasBegin)
  IF(ALLOCATED(vPhasFinal))  DEALLOCATE(vPhasFinal)
  IF(ALLOCATED(tPathData))   DEALLOCATE(tPathData)
  IF(ALLOCATED(vTPpath))     DEALLOCATE(vTPpath)
  IF(ALLOCATED(vPathLogK))   DEALLOCATE(vPathLogK)
ENDSUBROUTINE Path_Vars_Clean 

ENDMODULE M_Path_Vars

