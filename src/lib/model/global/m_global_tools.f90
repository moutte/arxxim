MODULE M_Global_Tools
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,fHtm,T_,Stop_,Pause_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Global_Zero
  PUBLIC:: Global_Clean
  PUBLIC:: Global_Species_Select
  PUBLIC:: Global_TP_Update
  PUBLIC:: Global_Show
  !
CONTAINS

SUBROUTINE Global_Show
!--
!--
  USE M_IoTools
  USE M_T_Phase,ONLY: T_Phase
  USE M_Global_Vars, ONLY: vSpc,vFas
  USE M_Global_Vars, ONLY: vMixFas,vMixModel
  USE M_Global_Vars, ONLY: vSolFas,vSolModel
  !~ & vEle,vSpc,vSpcDtb,vMixModel,vMixFas,vFas,vDiscretModel,vDiscretParam
  !
  INTEGER:: I,J,K
  INTEGER:: F
  TYPE(T_Phase):: fs
  
  CALL GetUnit(F)
  OPEN(F,FILE="global_show.txt")

  WRITE(F,'(A)') "NamFs,iSpc,iMix,iSol"
  DO I=1,SIZE(vFas)
  
    fs= vFas(I)
    
    WRITE(F,'(A,3I4,3X)',ADVANCE="NO") &
    & fs%NamFs,fs%iSpc,fs%iMix,fs%iSol
    
    IF(fs%iSpc >0) WRITE(F,'(A)') vSpc(fs%iSpc)%NamSp
    
    IF(fs%iMix >0) THEN
      DO J=1,vMixModel(vMixFas(fs%iMix)%iModel)%nPole
        K= vMixModel(vMixFas(fs%iMix)%iModel)%vIPole(J)
        WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vSpc(K)%NamSp)
      ENDDO
      WRITE(F,*)
    ENDIF
    
    IF(fs%iSol >0) THEN
      K= vSolModel(vSolFas(fs%iSol)%iModel)%iSolvent
      WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vSpc(K)%NamSp)
      DO J=1,vSolModel(vSolFas(fs%iSol)%iModel)%nSpecies
        K= vSolModel(vSolFas(fs%iSol)%iModel)%vISpecies(J)
        WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vSpc(K)%NamSp)
      ENDDO
      WRITE(F,*)
    ENDIF
    
  ENDDO
  
  CLOSE(F)
  
  RETURN
END SUBROUTINE Global_Show

SUBROUTINE Global_TP_Update( &
& TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
& vSpc,vMixModel,vMixFas,vFas) !inout
!--
!-- update (T,P) dependent properties of species, sol'models, phases,
!=  to be called everytime T,P changes
!--
  USE M_T_Species,  ONLY: T_Species,T_SpeciesDtb
  USE M_T_MixModel, ONLY: T_MixModel
  USE M_T_MixPhase, ONLY: T_MixPhase,MixPhase_CalcActivs
  USE M_T_MixModel, ONLY: MixModel_Margul_Wg_Calc
  USE M_T_Phase,    ONLY: T_Phase,Phase_Calc
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  !
  REAL(dp),            INTENT(IN):: TdgK,Pbar
  TYPE(T_SpeciesDtb),  INTENT(IN):: vSpcDtb(:)
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  TYPE(T_DiscretParam),INTENT(IN):: vDiscretParam(:)
  !
  TYPE(T_Species), INTENT(INOUT):: vSpc(:)
  TYPE(T_MixModel),INTENT(INOUT):: vMixModel(:)
  TYPE(T_MixPhase),INTENT(INOUT):: vMixFas(:)
  TYPE(T_Phase),   INTENT(INOUT):: vFas(:)
  !
  INTEGER :: I
  !
  CALL Species_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel)
  !
  !--------------- update activ's of end-members in non-aqu'solutions --
  DO I=1,SIZE(vMixFas) !
    CALL MixPhase_CalcActivs( &
    & TdgK,Pbar,                    & !IN
    & vMixModel(vMixFas(I)%iModel), & !IN
    & vMixFas(I))                     !INOUT
  ENDDO
  !
  !------------------------------------------------ update all phases --
  DO I=1,SIZE(vFas)
    CALL Phase_Calc( &
    & TdgK,Pbar,vSpc,vMixModel,vMixFas, & !IN
    & vFas(I))                            !OUT
  ENDDO
  !
ENDSUBROUTINE Global_TP_Update

SUBROUTINE Species_TP_Update( &
& TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
& vSpc,vMixModel) !inout
!--
!-- update (T,P) dependent properties of species, sol'models, phases,
!=  to be called everytime T,P changes
!--
  USE M_T_DtbH2OHkf,ONLY: DtbH2OHkf_Calc,T_DtbH2OHkf
  USE M_T_Species,  ONLY: T_Species,T_SpeciesDtb
  USE M_T_MixModel, ONLY: T_MixModel
  USE M_T_MixModel, ONLY: MixModel_Margul_Wg_Calc
  USE M_Dtb_Calc,   ONLY: Species_TP_Update_fromDtb
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  USE M_DiscretModel_Tools
  !
  REAL(dp),            INTENT(IN):: TdgK,Pbar
  TYPE(T_SpeciesDtb),  INTENT(IN):: vSpcDtb(:)
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  TYPE(T_DiscretParam),INTENT(IN):: vDiscretParam(:)
  !
  TYPE(T_Species), INTENT(INOUT):: vSpc(:)
  TYPE(T_MixModel),INTENT(INOUT):: vMixModel(:)
  !
  INTEGER :: I,J
  TYPE(T_DtbH2OHkf):: PropsH2O
  !
  ! compute solvent properties, needed for aqu'species
  IF(COUNT(vSpc(:)%Typ=='AQU')>0) &
  & CALL DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
  
  !-------------------- update T,P dependent param's for pure species --
  DO J=1,SIZE(vSpc)
    IF(vSpc(J)%iDtb>0) &
    CALL Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,vSpc(J))
    !-> update vSpc(J)%G0rt !-> =G/RT
  ENDDO
  
  !----------------- update T,P dependent parameters in mixing models --
  DO I=1,SIZE(vMixModel)
    IF(vMixModel(I)%NMarg>0) &
    & CALL MixModel_Margul_Wg_Calc( &
    & TdgK,Pbar, &  !in
    & vMixModel(I)) !out
  ENDDO
  
  !--------------------------- update species built by discretization --
  IF(SIZE(vDiscretModel)>0) THEN
    CALL DiscretSpecies_TP_Update( &
    & vMixModel,     & !IN
    & TdgK,Pbar,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
  ENDIF
  
  RETURN
END SUBROUTINE Species_TP_Update

SUBROUTINE Global_Zero
  USE M_Global_Vars
  
  nAq=0; nMn=0; nGs=0
  !
  CALL Global_Clean
  !
  ALLOCATE(vEle(0))
  !
  ALLOCATE(vSpc(0))
  ALLOCATE(vSpcDtb(0))
  !
  ALLOCATE(vFas(0))
  !
  ALLOCATE(vDiscretModel(0))
  ALLOCATE(vDiscretParam(0))
  !
  ALLOCATE(vKinModel(0))
  ALLOCATE(vKinFas(0))
  !
  ALLOCATE(vMixModel(0))
  ALLOCATE(vMixFas(0))
  !
  ALLOCATE(vSolModel(0))
  ALLOCATE(vSolFas(0))
  !
  ALLOCATE(tFormula(0,0))
  !
  !default solvent= water at 25°C/1bar (1 atm ???)
  !~ Solvent%Name=        "WATER"
  !~ Solvent%ActModel=    "IDEAL"
  !
ENDSUBROUTINE Global_Zero

SUBROUTINE Global_Clean
  USE M_Global_Vars
  USE M_Dtb_Vars
  !
  IF(ALLOCATED(vEle))    DEALLOCATE(vEle)
  !
  IF(ALLOCATED(vSpc))    DEALLOCATE(vSpc)
  IF(ALLOCATED(vSpcDtb)) DEALLOCATE(vSpcDtb)
  !
  IF(ALLOCATED(vFas))    DEALLOCATE(vFas)
  !
  IF(ALLOCATED(vMixModel)) DEALLOCATE(vMixModel)
  IF(ALLOCATED(vMixFas))   DEALLOCATE(vMixFas)
  !
  IF(ALLOCATED(vSolModel)) DEALLOCATE(vSolModel)
  IF(ALLOCATED(vSolFas))   DEALLOCATE(vSolFas)
  !
  IF(ALLOCATED(vDiscretModel)) DEALLOCATE(vDiscretModel)
  IF(ALLOCATED(vDiscretParam)) DEALLOCATE(vDiscretParam)
  !
  IF(ALLOCATED(vKinModel)) DEALLOCATE(vKinModel)
  IF(ALLOCATED(vKinFas))   DEALLOCATE(vKinFas)
  !
  IF(ALLOCATED(tFormula))  DEALLOCATE(tFormula)
  !
  CALL Dtb_Vars_Clean
  !
ENDSUBROUTINE Global_Clean

SUBROUTINE Global_Species_Select
  USE M_T_Species,  ONLY: T_Species
  USE M_Global_Vars,ONLY: vSpc
  !
  LOGICAL,ALLOCATABLE:: vExclude(:)
  LOGICAL,ALLOCATABLE:: vInclude(:)
  INTEGER:: iSp,nSp
  !
  TYPE(T_Species),ALLOCATABLE:: vSpcNew(:)
  !
  !------------------ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  ALLOCATE(vExclude(1:SIZE(vSpc)))  ;  vExclude(:)= .FALSE.
  ALLOCATE(vInclude(1:SIZE(vSpc)))  ;  vInclude(:)= .TRUE.
  CALL Species_Read_Excluded(vSpc,vExclude,vInclude)
  !-----------------/ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  !
  !------------------------------------------------ species selection --
  !-- from the vSpc built from database,
  !-- retrieve species consistent
  !-- with vEle & redox & vExclude
  IF(iDebug>0) WRITE(fTrc,'(A)') &
  & "< Restrict database to species consistent with element list and redox state"
  !
  ALLOCATE(vSpcNew(SIZE(vSpc)))
  !
  nSp=0
  DO iSp=1,SIZE(vSpc)
    IF( vExclude(iSp) .OR. (.NOT. vInclude(iSp)) ) CYCLE
    nSp=  nSp+1
    vSpcNew(nSp)= vSpc(iSp)
  ENDDO
  !
  DEALLOCATE(vExclude)
  DEALLOCATE(vInclude)
  !-----------------------------------------------/ species selection --
  !
  !-------------------------------------------- build new sorted vSpc --
  DEALLOCATE(vSpc)
  ALLOCATE(vSpc(1:nSp))
  vSpc(1:nSp)= vSpcNew(1:nSp)
  !
  RETURN
END SUBROUTINE Global_Species_Select

SUBROUTINE Species_Read_Excluded(vSpc,vExclude,vInclude)
!--
!-- process the SPECIES.EXCLUDE block
!-- and the SPECIES.INCLUDE block
!--
  USE M_Files,    ONLY: NamFInn
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_T_Species,ONLY: T_Species,Species_Index,Species_Rename
  !
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  LOGICAL,INTENT(OUT):: vExclude(:),vInclude(:)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios,I     
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Species_Read_Excluded"
  !
  vExclude= .FALSE.
  vInclude= .TRUE.
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO
  
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    
    CALL AppendToEnd(L,W,EoL)
    SELECT CASE(W)
    
      CASE("ENDINPUT"); EXIT DoFile
      
      CASE("SPECIES.EXCLUDE")
        !
        vExclude(:)= .FALSE.
        !
        DoLine1: DO
        
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoLine1 !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          
          SELECT CASE(W)
            CASE("ENDINPUT"); EXIT DoFile
            CASE("END","ENDSPECIES.EXCLUDE"); EXIT DoLine1
            !-> reading data from several BLOCK..END blocks 
            !CASE("END","ENDSPECIES.EXCLUDE"); EXIT DoFile
            !!-> reading only the first available BLOCK
          END SELECT
          !
          CALL Species_Rename(W)
          I= Species_Index(TRIM(W),vSpc)
          !
          !---------------------------------------------------- trace --
          IF(I<1) THEN
            WRITE(fTrc,'(3A)') "Species ",TRIM(W)," = is not in current vSpc"
          ELSE
            WRITE(fTrc,'(3A)') "Species ",TRIM(W)," = is excluded"
            vExclude(I)=.true.
          ENDIF
          !---------------------------------------------------/ trace --
          !
        ENDDO DoLine1
      !ENDCASE("SPECIES.EXCLUDE")
      
      CASE("SPECIES.INCLUDE")
        !
        vInclude(:)= .FALSE.
        !
        DoLine2: DO
        
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoLine2 !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          
          SELECT CASE(W)
            CASE("ENDINPUT"); EXIT DoFile
            CASE("END","ENDSPECIES.INCLUDE"); EXIT DoLine2
            !-> reading data from several BLOCK..END blocks 
            !CASE("END","ENDSPECIES.EXCLUDE"); EXIT DoFile
            !!-> reading only the first available BLOCK
          END SELECT
          
          CALL Species_Rename(W)
          I= Species_Index(TRIM(W),vSpc)
          !
          !---------------------------------------------------- trace --
          IF(I<1) THEN
            WRITE(fTrc,'(3A)') "Species ",TRIM(W)," = is not in current vSpc"
          ELSE
            WRITE(fTrc,'(3A)') "Species ",TRIM(W)," = is included"
            vInclude(I)=.true.
          ENDIF
          !---------------------------------------------------/ trace --
          !
        ENDDO DoLine2
      !ENDCASE("SPECIES.INCLUDE")
      
    END SELECT
  ENDDO DoFile
  !
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Species_Read_Excluded"
ENDSUBROUTINE Species_Read_Excluded

ENDMODULE M_Global_Tools
