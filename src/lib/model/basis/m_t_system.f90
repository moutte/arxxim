MODULE M_T_System

  USE M_Kinds
  USE M_T_Component,ONLY: T_Component
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: &
  & T_System,   &
  & System_New, &
  & System_Free
  
  TYPE:: T_System
    !PRIVATE
    CHARACTER(LEN=7):: SysName ! MASTER, INJECT, BOX, ...
    CHARACTER(LEN=7):: SysType ! AQUEOUS, ...
    !
    REAL(dp):: SysTdgK, SysPbar
    !
    INTEGER :: SysNCp
    TYPE(T_Component),ALLOCATABLE:: vCpn(:)
    !
  END TYPE T_System
  
CONTAINS

SUBROUTINE System_New(N,This)
  INTEGER,INTENT(IN):: N
  TYPE(T_System),INTENT(OUT):: This
  !
  This%SysName= "VOID"
  This%SysType= "VOID"
  
  This%SysTdgK= 0.D0
  This%SysPbar= 0.D0
  
  This%SysNCp=  N
  ALLOCATE(This%vCpn(N))
  !
END SUBROUTINE System_New

SUBROUTINE System_Free(This)
  TYPE(T_System),INTENT(INOUT):: This
  !
  This%SysName="VOID"
  This%SysType="VOID"
  This%SysNCp= 0
  DEALLOCATE(This%vCpn)
  !
END SUBROUTINE System_Free

END MODULE M_T_System

