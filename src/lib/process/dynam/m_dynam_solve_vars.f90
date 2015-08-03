MODULE M_Dynam_Solve_Vars
!--
!-- dynamic variables used only at "low level" (residual / jacobian)
!-- -> separated from other variable sets for dynam modules
!-- (M_Dynam_Vars,M_Dynam_OneStep_Vars)
!--
  USE M_Kinds
  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC:: Dynam_Solve_Vars_Alloc
  PUBLIC:: Dynam_Solve_Vars_Clean
  !
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: &
  !! & dVmAdLnXf, &   !1:nMk,1:nAq, @VmAct / @LnXf
  & dVmQdLnXf, &   !1:nMk,1:nCi, @VmQsK / @LnXf
  ! 
  & dSRdLnXm,  &   !1:nMk,1:nMk, derivative of vSurfK / Ln(nMol(Mineral))
  & dSRdXm,    &   !1:nMk,1:nMk, derivative of vSurfK / nMol(Mineral)
  !
  & dVMdLnXf,  &   !1:nMk,1:nCi, derivative of vRateM / Ln(nMol(Prim'Aqu.Spcecies))
  & dVMdLnXm,  &   !1:nMk,1:nMk, derivative of vRateM / Ln(nMol(Mineral))
  & dVMdXm         !1:nMk,1:nMk, derivative of vRateM / nMol(Mineral)
  !
CONTAINS

SUBROUTINE Dynam_Solve_Vars_Alloc(nMk,nCi)
  INTEGER, INTENT(IN):: nMk,nCi
  !
  ALLOCATE(dSRdLnXm(1:nMk,1:nMk))     ;  dSRdLnXm= Zero
  ALLOCATE(dSRdXm(1:nMk,1:nMk))       ;  dSRdXm=   Zero
  !
  ALLOCATE(dVmQdLnXf(1:nMk,1:nCi))    ;  dVmQdLnXf= Zero !VmAct depends on fluid species only
  !
  !! ALLOCATE(dVmAdLnXf(1:nMk,1:nCi)); dVmAdLnXf= Zero !VmAct depends on fluid species only
  !
  ALLOCATE(dVmdLnXf(1:nMk,1:nCi))     ;  dVmdLnXf= Zero !
  ALLOCATE(dVmdLnXm(1:nMk,1:nMk))     ;  dVmdLnXm= Zero !
  ALLOCATE(dVmdXm  (1:nMk,1:nMk))     ;  dVmdXm=   Zero !
  !
ENDSUBROUTINE Dynam_Solve_Vars_Alloc

SUBROUTINE Dynam_Solve_Vars_Clean

  IF(ALLOCATED(dSRdLnXm))   DEALLOCATE(dSRdLnXm)
  IF(ALLOCATED(dSRdXm))     DEALLOCATE(dSRdXm)
  !
  !! IF(ALLOCATED(dVmAdLnXf))  DEALLOCATE(dVmAdLnXf)      
  IF(ALLOCATED(dVmQdLnXf))  DEALLOCATE(dVmQdLnXf)
  !
  IF(ALLOCATED(dVmdXm))     DEALLOCATE(dVmdXm)
  IF(ALLOCATED(dVmdLnXf))   DEALLOCATE(dVmdLnXf) 
  IF(ALLOCATED(dVmdLnXm))   DEALLOCATE(dVmdLnXm)

ENDSUBROUTINE Dynam_Solve_Vars_Clean

ENDMODULE M_Dynam_Solve_Vars
