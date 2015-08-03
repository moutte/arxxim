MODULE M_Global_Alloc
  USE M_Kinds
  USE M_Trace,  ONLY: fTrc,T_,Stop_,iDebug
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Elements_Alloc_forDtb
  PUBLIC:: SpeciesDtb_Alloc
  PUBLIC:: MixModels_Alloc
  PUBLIC:: MixPhases_Alloc
  PUBLIC:: Phases_Alloc
  PUBLIC:: Phases_Alloc_New
  PUBLIC:: DiscretParam_Alloc
  !
CONTAINS

SUBROUTINE Elements_Alloc_forDtb
!--
!-- read all elements (name,valency, weight) from database -> builds vEle
!--
  USE M_T_Element,   ONLY: Element_Index
  USE M_Global_Vars, ONLY: vEle
  USE M_Element_Read,ONLY: T_LnkEle,Elements_BuildLnk,Elements_LnkToVec
  !
  TYPE(T_LnkEle),POINTER::LnkEle
  INTEGER:: nEle
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< Elements_Alloc_forDtb"
  !
  CALL Elements_BuildLnk(LnkEle,nEle)
  !
  IF(nEle>0) THEN
    !
    IF(ALLOCATED(vEle)) DEALLOCATE(vEle)
    ALLOCATE(vEle(1:nEle))
    !
    CALL Elements_LnkToVec(LnkEle,vEle)
    !
    IF(Element_Index("O__",vEle)==0) CALL Stop_("Oxygen not Found ???") !________STOP 
    IF(Element_Index("H__",vEle)==0) CALL Stop_("Hydrogen not Found ???") !______STOP 
    !
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,"(A,/)") "</ Elements_Alloc_forDtb"
ENDSUBROUTINE Elements_Alloc_forDtb

SUBROUTINE SpeciesDtb_Alloc !-> build vSpcDtb
  USE M_SpeciesDtb_Read,ONLY: T_LnkSpc,SpeciesDtb_LnkToVec,SpeciesDtb_BuildLnk
  !
  USE M_Global_Vars, ONLY: vSpcDtb
  !
  TYPE(T_LnkSpc), POINTER:: LnkSpc
  INTEGER:: N
  !
  CALL SpeciesDtb_BuildLnk(LnkSpc,N)
  !-> must come after Dtb_Read_HSV, which builds all vDtbXxx from HSV DATAbases
  !-> from vDtbXxxXxx tables, build linked list LnkSpc
  !
  IF(N>0) THEN !transfer LnkSpc to vSpcDtb
    IF(ALLOCATED(vSpcDtb)) DEALLOCATE(vSpcDtb)
    ALLOCATE(vSpcDtb(1:N))
    CALL SpeciesDtb_LnkToVec(LnkSpc,vSpcDtb)
  ENDIF
  !
  !~ IF(iDebug>3) THEN 
    !~ ALLOCATE(tFormul(1:SIZE(vEle),1:SIZE(vSpc)))
    !~ CALL FormulaTable_Calc(vEle,vSpc,tFormul)
    !~ CALL FormulaTable_Sho(vEle,vSpc,tFormul)
    !~ DEALLOCATE(tFormul)
  !~ ENDIF
  !
ENDSUBROUTINE SpeciesDtb_Alloc

SUBROUTINE MixModels_Alloc(vSpc) !-> build vMixModel
  USE M_T_Species,    ONLY: T_Species
  USE M_Global_Vars,  ONLY: vMixModel
  USE M_MixModel_Read,ONLY: T_LnkMixModel,MixModel_BuildLnk,MixModel_LnkToVec
  !
  TYPE(T_Species),DIMENSION(:),INTENT(IN)::vSpc
  !
  TYPE(T_LnkMixModel),POINTER:: LnkSol
  INTEGER::N
  LOGICAL:: Ok
  CHARACTER(LEN=80):: MsgError
  !
  CALL MixModel_BuildLnk(vSpc,LnkSol,N,Ok,MsgError)
  !
  IF(.NOT.Ok) THEN
    CALL Stop_("Error reading MIXTURE.MODEL: "//TRIM(MsgError))
  ENDIF
  !
  IF(N>0) THEN
    IF(ALLOCATED(vMixModel)) DEALLOCATE(vMixModel)
    ALLOCATE(vMixModel(1:N))
    CALL MixModel_LnkToVec(LnkSol,vMixModel)
  ENDIF
  !
ENDSUBROUTINE MixModels_Alloc

SUBROUTINE MixPhases_Alloc(vSpc,vMixModel) !-> build vMixFas
  USE M_T_Species,    ONLY: T_Species
  USE M_T_MixModel,   ONLY: T_MixModel
  USE M_T_MixPhase,   ONLY: MixPhase_Zero
  USE M_MixPhase_Read,ONLY: T_LnkFas,MixPhase_BuildLnk,MixPhase_LnkToVec
  !
  USE M_Global_Vars,  ONLY: vMixFas
  !
  TYPE(T_Species), DIMENSION(:),INTENT(IN) :: vSpc
  TYPE(T_MixModel),DIMENSION(:),INTENT(IN) :: vMixModel
  !
  TYPE(T_LnkFas),POINTER:: Lnk
  INTEGER::I,nMixFas
  LOGICAL:: Ok
  CHARACTER(LEN=80):: MsgError
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixPhases_Alloc"
  !
  CALL MixPhase_BuildLnk(vSpc,vMixModel,nMixFas,Lnk,Ok,MsgError)
  !
  IF(.NOT.Ok) THEN
    CALL Stop_("Error reading MIXTURE block, "//TRIM(MsgError))
  ENDIF
  !
  IF(nMixFas>0) THEN
    !
    IF(ALLOCATED(vMixFas)) DEALLOCATE(vMixFas)
    ALLOCATE(vMixFas(1:nMixFas))
    !
    DO I=1,nMixFas; CALL MixPhase_Zero(vMixFas(I)); ENDDO
    !
    CALL MixPhase_LnkToVec(Lnk,vMixFas)
    !
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ MixPhases_Alloc"
ENDSUBROUTINE MixPhases_Alloc

SUBROUTINE Phases_Alloc(vSpc,vMixFas)
!--
!-- -> build vFas
!-- vFas- (H2O) U (non-aqueous pure species) U (solution phases)
!--
  USE M_T_Species,  ONLY: T_Species,Species_Index
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_Global_Vars,ONLY: vFas
  USE M_T_Phase,    ONLY: Phase_Zero
  !
  TYPE(T_Species), DIMENSION(:),INTENT(IN):: vSpc
  TYPE(T_MixPhase),DIMENSION(:),INTENT(IN):: vMixFas
  !
  INTEGER:: I,N
  !
  ! calculate number of non-aqueous species -> included as "PURE" phase
  N=0
  IF(Species_Index("H2O",vSpc)/=0) N=N+1
  ! "H2O" is "AQU" but it should be included here as pure
  N= N &
  & + COUNT(vSpc(:)%Typ/="AQU") !& !skip the solute species
  !& + SIZE(vMixFas)
  !
  IF(ALLOCATED(vFas)) DEALLOCATE(vFas)
  ALLOCATE(vFas(1:N)); vFas(:)=Phase_Zero
  !
  N=0
  !
  IF(Species_Index("H2O",vSpc)/=0) THEN
    ! include H2O species as pure phase in phase table
    N=N+1
    vFas(N)%NamFs= "H2O"
    !! vFas(N)%Typ=   "PURE"
    vFas(N)%iSpc=  Species_Index("H2O",vSpc)
    vFas(N)%iMix=  0
    vFas(N)%iSol=  0
  ENDIF
  !
  DO I=1,SIZE(vSpc)
    ! include all non'aqu thermodyn'species and discretized species
    ! as pure phases in phase table
    IF(vSpc(I)%Typ/="AQU") THEN
      N=N+1
      vFas(N)%NamFs= TRIM(vSpc(I)%NamSp)
      !! vFas(N)%Typ=  "PURE"
      !! IF(vSpc(I)%iDiscret>0) vFas(N)%Typ= "DISCRET"
      vFas(N)%iSpc= I
      vFas(N)%iSol= 0
      vFas(N)%iMix= 0
    ENDIF
  ENDDO
  !
  ! DO I=1,SIZE(vMixFas)
  !   ! include all non'aqu solution phases in phase table
  !   vFas(N+I)%NamFs= TRIM(vMixFas(I)%Name)
  !   vFas(N+I)%Typ=   "MIXT"
  !   vFas(N+I)%iSol=  0
  !   vFas(N+I)%iSpc=  0
  !   vFas(N+I)%iMix=  I
  ! ENDDO
  
  ! write(11,'(A)') "<==Phases_Alloc====="
  ! do i=1,size(vFas)
  !   write(11,'(I3,1X,A)') i, vFas(i)%NamFs
  ! end do
  ! write(11,'(A)') "===================>"
  
  RETURN
END SUBROUTINE Phases_Alloc

SUBROUTINE Phases_Alloc_New(vSpc,vMixFas,vSolFas)
!--
!-- -> build vFas
!-- vFas- (non-aqueous pure species) U (mixture phases) U (solution phases)
!--
  USE M_T_Species,  ONLY: T_Species,Species_Index
  USE M_T_SolPhase, ONLY: T_SolPhase
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_Global_Vars,ONLY: vFas
  USE M_T_Phase,    ONLY: Phase_Zero
  !
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  TYPE(T_MixPhase),INTENT(IN):: vMixFas(:)
  TYPE(T_SolPhase),INTENT(IN):: vSolFas(:)
  !
  INTEGER:: I,N
  !
  ! calculate number of non-aqueous species -> included as "PURE" phase
  N= COUNT(vSpc(:)%Typ/="AQU") &
  &  + SIZE(vMixFas) &
  &  + SIZE(vSolFas)
  !
  !~ ! calculate number of non-aqueous species -> included as "PURE" phase
  !~ N=0
  !~ IF(Species_Index("H2O",vSpc)/=0) N=N+1
  !~ ! "H2O" is "AQU" but it should be included here as pure
  !~ N= N &
  !~ & + COUNT(vSpc(:)%Typ/="AQU") & !skip the solute species
  !~ & + SIZE(vMixFas)
  !
  IF(ALLOCATED(vFas)) DEALLOCATE(vFas)
  ALLOCATE(vFas(1:N)); vFas(:)=Phase_Zero
  !
  N=0
  !
  !~ IF(Species_Index("H2O",vSpc)/=0) THEN
    !~ ! include H2O species as pure phase in phase table
    !~ N=N+1
    !~ vFas(N)%NamFs= "H2O"
    !~ vFas(N)%Typ=   "PURE"
    !~ vFas(N)%iSpc=  Species_Index("H2O",vSpc)
  !~ ENDIF
  !
  DO I=1,SIZE(vSpc)
    ! include all non'aqu thermodyn'species (and discretized species)
    ! as pure phases in phase table
    IF(vSpc(I)%Typ/="AQU") THEN
      N= N+1
      vFas(N)%NamFs= TRIM(vSpc(I)%NamSp)
      !! vFas(N)%Typ=  "PURE"
      ! nota:
      !   discretized phases are comprised within the "pure phases"
      vFas(N)%iSpc= I
      vFas(N)%iSol= 0
      vFas(N)%iMix= 0
    ENDIF
  ENDDO
  !
  DO I=1,SIZE(vSolFas)
    ! include all non'aqu solution phases in phase table
    N= N+1
    vFas(N)%NamFs= TRIM(vSolFas(I)%Name)
    !! vFas(N)%Typ=  "SOLU"
    vFas(N)%iSol= I
    vFas(N)%iSpc= 0
    vFas(N)%iMix= 0
  ENDDO
  !
  DO I=1,SIZE(vMixFas)
    ! include all non'aqu solution phases in phase table
    N= N+1
    vFas(N)%NamFs= TRIM(vMixFas(I)%Name)
    !! vFas(N)%Typ=  "MIXT"
    vFas(N)%iSol= 0
    vFas(N)%iSpc= 0
    vFas(N)%iMix= I
  ENDDO
  
  ! write(11,'(A)') "<==Phases_Alloc_New=="
  ! do i=1,size(vFas)
  !   write(11,'(I3,1X,A)') i, vFas(i)%NamFs
  ! end do
  ! write(11,'(A)') "====================>"
  
  RETURN
END SUBROUTINE Phases_Alloc_New

SUBROUTINE DiscretParam_Alloc(vDiscretModel)
!--
!-- build vDiscretParam, array of T_DiscretParam
!-- that will be used to build pseudo-species derived from discretized mixtures(s)
!--
  USE M_T_MixModel,ONLY: T_MixModel
  USE M_T_DiscretModel
  !
  USE M_Global_Vars, ONLY: vDiscretParam
  !
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  !
  INTEGER:: N,iM
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretParam_Alloc"
  !
  N= 0
  DO iM= 1, SIZE(vDiscretModel)
    N= N + vDiscretModel(iM)%DimTot
  ENDDO
  DEALLOCATE(vDiscretParam)
  ALLOCATE(vDiscretParam(N))
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretParam_Alloc"
  !
ENDSUBROUTINE DiscretParam_Alloc

ENDMODULE M_Global_Alloc

