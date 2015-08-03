MODULE M_SolModel_Alloc
  USE M_Kinds
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: SolModel_Alloc
  PUBLIC:: SolPhase_Alloc
  !
CONTAINS

SUBROUTINE SolModel_Alloc
  USE M_T_SolModel
  USE M_Global_Vars,ONLY: vSolModel
  !
  INTEGER:: nSolModel
  INTEGER:: I
  
  nSolModel= 1
  
  IF(ALLOCATED(vSolModel)) DEALLOCATE(vSolModel)
  ALLOCATE(vSolModel(nSolModel))
  
  DO I=1,nSolModel
    !~ vSolModel(I)%Name= TRIM(vSolModelName(I))
    vSolModel(I)%Typ=  "LIQ"
    !~ vSolModel(I)%ActModel= "IDEAL"
    vSolModel(I)%iActModel= 1
  ENDDO
  
END SUBROUTINE SolModel_Alloc

SUBROUTINE SolPhase_Alloc
  USE M_T_SolPhase
  USE M_Global_Vars,ONLY: vSolFas
  !
  IF(ALLOCATED(vSolFas)) DEALLOCATE(vSolFas)
  ALLOCATE(vSolFas(1))
  
  vSolFas(1)%Name=   "WATER"
  vSolFas(1)%iModel= 1
  
END SUBROUTINE SolPhase_Alloc

ENDMODULE M_SolModel_Alloc
