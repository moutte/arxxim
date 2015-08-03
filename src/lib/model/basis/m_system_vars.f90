MODULE M_System_Vars
!--
!-- variables for system definition
!--
  USE M_Kinds
  !
  USE M_T_System
  USE M_T_Component,ONLY: T_Component
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: System_Zero
  PUBLIC:: System_Clean
  !
  REAL(dp),PUBLIC:: TdgK, Pbar
  !
  CHARACTER(LEN=7),PUBLIC:: System_Type="AQUEOUS" !AQUEOUS, SIMPLEX, MOMAS, GLOBAL
  !
  LOGICAL,PUBLIC:: BufferIsExtern= .FALSE.  !.TRUE.  !
  LOGICAL,PUBLIC:: CpnIsSpc= .TRUE.
  !
  TYPE(T_Component),ALLOCATABLE,PUBLIC:: vCpn   (:) !-> default system
  TYPE(T_Component),ALLOCATABLE,PUBLIC:: vCpnTot(:) !-> "global" system
  !
  ! stoikiometry table of the species in terms of the components
  REAL(dp),ALLOCATABLE,PUBLIC:: tStoikioCpn(:,:)
  !
  TYPE(T_System):: SysDefault,SysTotal,SysInject,SysBox
  !
  !-- vars for GEM computations
  !-- ? should be in a specific M_GEM_Vars module ?
  !~ TYPE(T_Component),ALLOCATABLE,PUBLIC:: vCpnGEM(:)
  !~ REAL(dp),         ALLOCATABLE,PUBLIC:: tStoikioGEM(:,:)
  
CONTAINS

SUBROUTINE System_Zero
  USE M_Dtb_Const,ONLY: Tref,Pref
  
  CALL System_Clean
  
  ALLOCATE(vCpn(0))
  ALLOCATE(vCpnTot(0))
  !ALLOCATE(vCpnInj(0))
  !ALLOCATE(vCpnBox(0))
  
  TdgK= Tref
  Pbar= Pref
  
ENDSUBROUTINE System_Zero

SUBROUTINE System_Clean
  IF(ALLOCATED(vCpn))     DEALLOCATE(vCpn)
  IF(ALLOCATED(vCpnTot))  DEALLOCATE(vCpnTot)
  IF(ALLOCATED(tStoikioCpn)) DEALLOCATE(tStoikioCpn)
ENDSUBROUTINE System_Clean

ENDMODULE M_System_Vars



