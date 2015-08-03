MODULE M_Equil_Specia
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Specia
  !
CONTAINS

SUBROUTINE Equil_Specia(iErr)
!--
!-- calculate equilibrium speciation of the aqueous phase
!-- (i.e. no consideration of saturation with minerals or gases)
!--
  USE M_Equil_Solve,ONLY: Equil_Solve
  !
  INTEGER,INTENT(OUT):: iErr
  !
  INTEGER :: nIts
  !~ REAL(dp):: CpuBegin,CpuEnd
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Specia >"
  !~ IF(iDebug>0) CALL CPU_TIME(CpuBegin)
  !
  CALL Equil_Solve(nIts,iErr)
  !
  !~ IF(iDebug>0) CALL CPU_TIME(CpuEnd)
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A,G15.6,/)') "CPU-Time",1.D3*(CpuEnd-CpuBegin)
  
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Specia >"
  
  RETURN
ENDSUBROUTINE Equil_Specia

ENDMODULE M_Equil_Specia


