MODULE M_SolModel_Calc
!--
!-- calculations of activity coefficients
!-- of solutes and solvent in aqueous solution
!-- to be applicable to any molality based solution model
!--
  USE M_Kinds
  USE M_Trace,ONLY: Stop_,fTrc,T_,iDebug
  USE M_T_DtbLogKTbl,ONLY: DimLogK_Max
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Solmodel_CalcGamma
  
CONTAINS

SUBROUTINE Solmodel_CalcGamma( &
& TdgK,Pbar,        & !IN
& SolModel,         & !IN
& isW,vSpcAq,vLAx,  & !IN
& vMolF,            & !IN
& vLnAct,vLnGam,    & !INOUT
& vTooLow_,OsmoSv) !OUT

  USE M_Solmodel_Calc_System
  
  USE M_T_Species,ONLY: T_Species
  USE M_T_SolModel,ONLY: T_SolModel
  !--
  TYPE(T_Species),INTENT(IN):: vSpcAq(:)
  TYPE(T_SolModel),INTENT(IN):: SolModel
  REAL(dp),INTENT(IN)   :: TdgK,Pbar
  INTEGER, INTENT(IN)   :: isW
  LOGICAL, INTENT(IN)   :: vLAx(:)  
  REAL(dp),INTENT(INOUT):: vMolF(:) 
  REAL(dp),INTENT(INOUT):: vLnAct(:) 
  REAL(dp),INTENT(INOUT):: vLnGam(:)
  LOGICAL, INTENT(OUT)  :: vTooLow_(:)
  REAL(dp),INTENT(OUT)  :: OsmoSv
  !---
  LOGICAL :: OkNew = .false. ! .true. => switch to new version
  !-----------

  IF (OkNew) THEN 
    CALL Solmodel_Calc_System(&
    & TdgK,Pbar,        & !IN
    & SolModel,         & !IN
    & isW,vSpcAq,vLAx,  & !IN
    & vMolF,            & !IN
    & vLnAct,vLnGam,    & !INOUT
    & vTooLow_,OsmoSv) !OUT
  ELSE
    CALL Solmodel_CalcGamma_Old(&
    & TdgK,Pbar,        & !IN
    & SolModel,         & !IN
    & isW,vSpcAq,vLAx,  & !IN
    & vMolF,            & !IN
    & vLnAct,vLnGam,    & !INOUT
    & vTooLow_,OsmoSv) !OUT
  END IF

END SUBROUTINE Solmodel_CalcGamma

SUBROUTINE Solmodel_CalcGamma_Old( &
!--
!-- WARNING
!-- when implementing a new activity model,
!-- add its name in the list of activity models
!--
& TdgK,Pbar,        & !IN
& SolModel,         & !IN
& isW,vSpcAq,vLAx,  & !IN
& vMolF,            & !IN
& vLnAct,vLnGam,    & !INOUT
& vTooLow_,OsmoSv) !OUT

!.IN = vXPi,vXAs=  MolNumbers of INERT aqu.species in system
!                  (i.e. for 1 kg water or 1 box)
!.OUT= vLnAct & vLnGam

  USE M_IOTools
  USE M_Numeric_Const,ONLY: TinyDP,Ln10,MinExpDP
  USE M_Dtb_Const, ONLY: T_CK
  USE M_T_Tpcond,  ONLY: T_TPCond
  USE M_T_Species, ONLY: T_Species
  USE M_T_SolModel,ONLY: T_SolModel
  
  USE M_Solmodel_Calc_Pitzer
  USE M_Solmodel_Calc_Debye_Hueckel
  USE M_Solmodel_Calc_Davies
  !~ USE M_Solmodel_Calc_SIT
  USE M_Solmodel_Calc_Dilute
  
  TYPE(T_Species), INTENT(IN):: vSpcAq(:)
  TYPE(T_SolModel),INTENT(IN):: SolModel
  
  REAL(dp),INTENT(IN)   :: TdgK,Pbar
  INTEGER, INTENT(IN)   :: isW
  LOGICAL, INTENT(IN)   :: vLAx(:)   !mobile aqu'species
  
  REAL(dp),INTENT(INOUT):: vMolF(:)  !mole numbers
  REAL(dp),INTENT(INOUT):: vLnAct(:) 
  REAL(dp),INTENT(INOUT):: vLnGam(:)
  
  LOGICAL, INTENT(OUT)  :: vTooLow_(:)
  REAL(dp),INTENT(OUT)  :: OsmoSv !SolModel's osmotic coeff'
  
  ! vMolF=  input for inert species, output for mobile species
  ! vLnAct= activities, input for mobile species, output for others
  ! activity coeff', input for mobile species, output for all (includ'g mobile)
  
  TYPE(T_Species),DIMENSION(:),ALLOCATABLE:: vSpcSolut
  REAL(dp),       DIMENSION(:),ALLOCATABLE:: vMolal,vLnGamSolut,vA0
  INTEGER,        DIMENSION(:),ALLOCATABLE:: vZSp
  
  REAL(dp):: LnActSv
  REAL(dp):: MolWeitSv,Rho,Eps,dhA,dhB,bDot !,Pbar
  INTEGER :: nAq,iAq,N,J !,nAi
  INTEGER :: iModel
  
  MolWeitSv= vSpcAq(isW)%WeitKg
  
  nAq= SIZE(vMolF)
  !nAi= nAq-nAx !nr of non mobile au'species
  vTooLow_=.FALSE.
  
  dhA= SolModel%Dat%dhA
  dhB= SolModel%Dat%dhB
  bDot=SolModel%Dat%bDot
  Rho= SolModel%Dat%Rho
  Eps= SolModel%Dat%Eps
  
  N= COUNT(vSpcAq(:)%Typ=="AQU") - 1 !mod 19/06/2008 09:38
  !-> solute species only
  ALLOCATE(vSpcSolut(N))
  ALLOCATE(vMolal(N))
  ALLOCATE(vLnGamSolut(N))
  ALLOCATE(vA0(N),vZsp(N))
  
  ! build arrays for solute species:
  ! vMolal (= molalities), vSpcSolut, vLnGamSolut
  ! Molality(iAq)=MoleNumber(iAq)/MassSolvent; vMolal(isW)=1/MolWeitSv=55.51
  ! vSpcSolut is used mainly as input for Pitzer calculations
  ! (separation solute / SolModel)
  J= 0
  DO iAq=1,SIZE(vSpcAq)
    IF(vSpcAq(iAq)%Typ=="AQU" .AND. vSpcAq(iAq)%NamSp/="H2O") THEN
      J= J+1
      vSpcSolut(J)= vSpcAq(iAq)
      vLnGamSolut(J)= vLnGam(iAq)
      !! vMolalSolut(J)= vMolal(iAq)
      IF(vLAx(iAq)) THEN 
        vMolal(J)= & !for mobile aqu'species, activity=IN
        & EXP(vLnAct(iAq) - vLnGam(iAq))
        vMolF(iAq)= vMolal(J) *vMolF(isW) *MolWeitSv
      ELSE !for inert aqu'species; MolWeitSv=MolWeitH2O
        vMolal(J)= vMolF(iAq) &
        & /vMolF(isW) &
        & /MolWeitSv
      ENDIF
    ENDIF
  ENDDO
  
  vZSp(:)= vSpcSolut(:)%Z
  vA0(:)=  vSpcSolut(:)%AquSize
  WHERE(vA0(:)==Zero .and. vZsp(:)/=0) vA0(:)= 3.72D0 ! provisional !!!!
  WHERE(vZsp(:)==0) vA0(:)= Zero
  
  vLnGam=Zero
  OsmoSv=One
  
  ! "IDEAL  ",              & ! 1
  ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
  ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
  ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
  ! "PITZER ", "SAMSON ",   & ! 8, 9
  ! "HKF81  ", "SIT    ",   & ! 10,11
  ! "name12 "               & ! 12
  
  !~ sModel= TRIM(SolModel%ActModel)
  iModel=SolModel%iActModel 
  
  SELECT CASE(iModel)

  CASE(1) !("IDEAL", actually should be "DILUTE" !!)
    CALL Solmodel_Calc_Dilute( &
    & vLnGamSolut, LnActSv, OsmoSv)
     
  CASE(2,3,4,5,10) !("DH1","DH2","DH1EQ3","DH2EQ3","HKF81")
    CALL Solmodel_Calc_Debye_Hueckel( &
    & TdgK, Pbar, vZSp, vMolal, &
    & dhA, dhB, vA0, bDot, iModel, MolWeitSv, &
    & vLnGamSolut, LnActSv, OsmoSv)
    
  CASE(8) !("PITZER")
    CALL Solmodel_Calc_Pitzer( &
    & vSpcSolut, &
    & MolWeitSv, &
    & Rho,Eps,TdgK,vMolal,&
    & vLnGamSolut,LnActSv,OsmoSv)
    
  CASE(6,7,9) !("DAV_1","DAV_2","SAMSON")
    CALL Solmodel_Calc_Davies( &
    & vZSp, vMolal, &
    & dhA, iModel, MolWeitSv, &
    & vLnGamSolut, LnActSv, OsmoSv)
  
  !~ CASE("SIT")
    !~ CALL Solmodel_Calc_SIT( &
    !~ & vZSp, vMolal, &
    !~ & dhA, sModel, MolWeitSv, &
    !~ & vLnGamSolut, LnActSv, OsmoSv)
    
  END SELECT
  
  !-------- back substitute solute species from vLnGamSolut in vLnGam --
  J= 0
  DO iAq=1,SIZE(vSpcAq)
    IF(vSpcAq(iAq)%Typ=="AQU" .AND. vSpcAq(iAq)%NamSp/="H2O") THEN
      J= J+1
      vLnGam(iAq)= vLnGamSolut(J)
      !
      IF(vLAx(iAq)) THEN 
        !-- MOBILE SPECIES --
        !-- -> apply act'coeff' to calculate #mole from activity
        IF(vLnAct(iAq) -vLnGam(iAq) <=MinExpDP) THEN
          vTooLow_(iAq)=.TRUE.
          vMolal(J)= EXP(MinExpDP)
        ELSE 
          vTooLow_(iAq)=.FALSE.
          vMolal(J)= EXP(vLnAct(iAq) -vLnGamSolut(J))
        ENDIF
        vMolF(iAq)= vMolal(J) *vMolF(isW) *MolWeitSv
      ELSE
        !-- INERT SPECIES --
        !-- -> compute activity
        IF(vMolal(J) <=TinyDP) THEN
          vTooLow_(iAq)=.TRUE.
          vLnAct(iAq)= LOG(TinyDP)
        ELSE
          vTooLow_(iAq)=.FALSE.
          vLnAct(iAq)= LOG(vMolal(J)) + vLnGamSolut(J)
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  !---------------------------------------------------/solute species--
  
  !--------------------------------------------------------- SolModel --
  vLnAct(isW)= LnActSv
  !IF(sModel /= "IDEAL") THEN
  IF(iModel /= 1) THEN
    vLnGam(isW)= LnActSv + LOG( One + MolWeitSv *SUM(vMolal(:)) )
  ELSE
    vLnGam(isW)= Zero
  ENDIF
  ! for SolModel, Gamma is the rational activity coeff', generally called Lambda
  !--------------------------------------------------------/ SolModel --
  
  DEALLOCATE(vSpcSolut)
  DEALLOCATE(vMolal)
  DEALLOCATE(vLnGamSolut)
  DEALLOCATE(vA0,vZsp)

ENDSUBROUTINE Solmodel_CalcGamma_Old

ENDMODULE M_SolModel_Calc

