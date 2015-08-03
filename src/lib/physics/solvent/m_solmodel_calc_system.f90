MODULE M_Solmodel_Calc_System

  !---------------------------------------------------------------
  ! Test for a New organization to compute Gamma
  ! Solvent Activity Calc Interface for all Models   
  ! Compute intermediate properties and CALL SUBROUTINEs 
  !---------------------------------------------------------------

  USE M_Kinds
  USE M_Trace

  IMPLICIT NONE
  PRIVATE 

  PUBLIC :: Solmodel_Calc_System
  PUBLIC :: Solmodel_Calc_Physics

CONTAINS
  
  !-- WARNING
  !-- when implementing a new activity model,
  !-- add its name to the list of activity model names
  
  SUBROUTINE Solmodel_Calc_System(&
  & TdgK,Pbar,        & !IN
  & SolModel,         & !IN  ! solvent phase model
  & isW,vSpcAq,vLAx,  & !IN
  & vMolF,            & !IN
  & vLnAct,vLnGam,    & !INOUT
  & vTooLow_,OsmoSv)    !OUT
  
    ! IN = vXPi,vXAs
    !    = MolNumbers of INERT aqu.species in system
    !      (i.e. for 1 kg water or 1 box)
    ! OUT= vLnAct & vLnGam
    
    USE M_IOTools
    USE M_Numeric_Const,ONLY: TinyDP,Ln10,MinExpDP
    USE M_Dtb_Const,ONLY: T_CK
    USE M_T_Tpcond, ONLY: T_TPCond
    USE M_T_Species,ONLY: T_Species
    USE M_T_SolModel
    USE M_Safe_Functions
    !--
    IMPLICIT NONE
    !
    TYPE(T_Species), INTENT(IN):: vSpcAq(:)
    TYPE(T_SolModel),INTENT(IN):: SolModel
    REAL(dp),INTENT(IN)   :: TdgK,Pbar
    INTEGER, INTENT(IN)   :: isW
    LOGICAL, INTENT(IN)   :: vLAx(:)  !mobile aqu'species
    REAL(dp),INTENT(INOUT):: vMolF(:)  !mole numbers
    REAL(dp),INTENT(INOUT):: vLnAct(:) 
    REAL(dp),INTENT(INOUT):: vLnGam(:)
    LOGICAL, INTENT(OUT)  :: vTooLow_(:)
    REAL(dp),INTENT(OUT)  :: OsmoSv !solvent's osmotic coeff'
    !vMolF=  input for inert species, output for mobile species
    !vLnAct= activities, input for mobile species, output for others
    !activity coeff', input for mobile species, output for all (includ'g mobile)
    !
    !! TYPE(T_Species),DIMENSION(:),ALLOCATABLE:: vSpcSolut
    !! REAL(dp),       DIMENSION(:),ALLOCATABLE:: vMolal,vLnGamSolut,vA0
    !! INTEGER,        DIMENSION(:),ALLOCATABLE:: vZSp
    !! REAL(dp):: TdgK
    !! REAL(dp):: LnActSv
    REAL(dp):: MolWeitSv
    INTEGER :: nAq,iAq,J
    !  
    !--
    ! compute molality for external species and test vTooLow_
    !--
    CALL Solmodel_Calc_Physics( &
    & TdgK,Pbar,        & !IN  ! temperature and pressure
    & SolModel,         & !IN  ! solvent phase model
    & vSpcAq,           & !IN  ! vector of Aqueous Species
    & isW,              & !IN  ! index of solvent
    & vMolF,            & !IN  ! mole numbers
    !--------------------------
    & vLnGam,           & !OUT ! gamma coefficient
    & vLnAct,           & !OUT ! activity
    & OsmoSv )            !OUT ! osmotic coefficient 
    !---
    nAq = SIZE(vMolF)
    MolWeitSv = vSpcAq(isW)%WeitKg
    !
    !------------ vMolF of External Species
    !
    J= 0 
    DO iAq=1,nAq
      IF (.not.(iAq==isW)) THEN
        J= J+1
        IF(vLAx(iAq)) THEN 
          vMolF(iAq)= FSafe_Exp(vLnAct(iAq)-vLnGam(iAq)) &
          &         * vMolF(isW)*MolWeitSv
        END IF
      END IF
    END DO
    !
    !------------ vTooLow Logical
    !
    vTooLow_ = .false.
    J= 0 
    DO iAq=1,nAq
      IF (.not.(iAq==isW)) THEN
        IF( ( vLnAct(iAq) - vLnGam(iAq) ) <= MinExpDP) THEN
          vTooLow_(iAq)=.TRUE.
        ELSE 
          vTooLow_(iAq)=.FALSE.
        ENDIF
      END IF
    ENDDO
    
  END SUBROUTINE Solmodel_Calc_System
  !---
  SUBROUTINE Solmodel_Calc_Physics( &
  & TdgK,Pbar,        & !IN  ! temperature and pressure
  & SolModel,         & !IN  ! sol'model activ'model
  & vSpcAq,           & !IN  ! vector of Aqueous Species
  & isW,              & !IN  ! index of solvent
  & vMolF,            & !IN  ! mole numbers
  !--------------------------
  & vLnGam,           & !OUT ! gamma coefficient
  & vLnAct,           & !OUT ! activity
  & OsmoSv )            !OUT ! osmotic coefficient 
  !
  USE M_IOTools
  USE M_Safe_Functions
  !--
  USE M_Solmodel_Calc_Pitzer
  USE M_Solmodel_Calc_Debye_Hueckel
  USE M_Solmodel_Calc_Davies
  !~ USE M_Solmodel_Calc_SIT
  USE M_Solmodel_Calc_Dilute
  !--
  USE M_T_Species
  USE M_T_SolModel
  !----
  IMPLICIT NONE
  !
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  TYPE(T_SolModel),INTENT(IN):: SolModel
  TYPE(T_Species), INTENT(IN):: vSpcAq(:)
  INTEGER,         INTENT(IN):: isW       
  REAL(dp),        INTENT(IN):: vMolF(:)  ! mole numbers
  REAL(dp),INTENT(OUT)  :: vLnGam(:) ! gamma coefficient
  REAL(dp),INTENT(OUT)  :: vLnAct(:) ! activity
  REAL(dp),INTENT(OUT)  :: OsmoSv    ! osmotic coefficient
  !---
  INTEGER :: nAq, nSolute, iAq, iSolute
  REAL(dp), DIMENSION(:),ALLOCATABLE:: vMolal
  REAL(dp), DIMENSION(:),ALLOCATABLE:: vLnGamSolute
  REAL(dp), DIMENSION(:),ALLOCATABLE:: vA0
  INTEGER,  DIMENSION(:),ALLOCATABLE:: IdxSolute
  INTEGER,  DIMENSION(:),ALLOCATABLE:: vZSp
  !---
  REAL(dp):: MolWeitSv, MassSv
  REAL(dp):: LnActSv
  REAL(dp):: dhA, dhB, bDot, Rho, Eps
    !
    !---
    !// Retrieve Debye-Hueckel Coeffcients
    dhA=  SolModel%Dat%dhA
    dhB=  SolModel%Dat%dhB
    bDot= SolModel%Dat%bDot
    
    Rho=  SolModel%Dat%Rho
    Eps=  SolModel%Dat%Eps
    !
    !// Build solutes indices   
    nAq = SIZE(vMolF)
    nSolute = nAq - 1 
    ALLOCATE(IdxSolute(nSolute))
    !
    iSolute = 0
    DO iAq = 1, nAq
      IF (.not.(iAq==isW)) THEN
        iSolute = iSolute+1
        IdxSolute(iSolute) = iAq
      END IF
    END DO
    
    !// Compute solvent properties
    MolWeitSv = vSpcAq(isW)%WeitKg
    MassSv = MolWeitSv * vMolF(isW)
    OsmoSv = One ! default value

    !// Compute solute properties
    ALLOCATE(vMolal(1:nSolute))
    vMolal(1:nSolute) = vMolF(IdxSolute(:))/MassSv
    !
    ALLOCATE(vZSp(1:nSolute))
    vZSp(1:nSolute)= vSpcAq(IdxSolute(:))%Z
    !
    ALLOCATE(vA0(1:nSolute))
    vA0(:)= vSpcAq(IdxSolute(:))%AquSize
    WHERE(vA0(:)==Zero .and. vZSp(:)/=0) vA0(:)= 3.72D0 ! provisional !!!
    !
    ALLOCATE(vLnGamSolute(1:nSolute))
    vLnGamSolute = Zero

    !// Compute activities
    
    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12
    
    ! SELECT CASE(SolModel%ActModel)
    SELECT CASE(SolModel%iActModel)

    CASE(1) ! ("IDEAL")
      CALL Solmodel_Calc_Dilute( &
      !! & vMolal,MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv) 

    CASE(8) !("PITZER")
      CALL Solmodel_Calc_Pitzer( &
      & vSpcAq(IdxSolute), &
      & MolWeitSv, &
      & Rho,Eps,TdgK,vMolal,&
      & vLnGamSolute, LnActSv, OsmoSv)

    CASE(2,3,4,5) !("DH1","DH2","DH1EQ3","DH2EQ3")
      CALL Solmodel_Calc_Debye_Hueckel( &
      & TdgK, Pbar, vZSp, vMolal, &
      & dhA, dhB, vA0, bDot, SolModel%iActModel, MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv)

    CASE(6,7,9) !("DAV_1","DAV_2","SAMSON")
      CALL Solmodel_Calc_Davies( &
      & vZSp, vMolal, &
      & dhA, SolModel%iActModel, MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv)

    !~ CASE(11) !("SIT")
      !~ CALL Solmodel_Calc_SIT( &
      !~ & vZSp, vMolal, &
      !~ & dhA, SolModel%ActModel, MolWeitSv, &
      !~ & vLnGamSolute, LnActSv, OsmoSv)
    
    END SELECT

    !// Store activity and gamma for solute species 
    vLnGam(IdxSolute(:)) = vLnGamSolute(1:nSolute)
    vLnAct(IdxSolute(:)) = FSafe_vLog(vMolal(1:nSolute)) + vLnGamSolute(1:nSolute) 

    !// Store activity and gamma for solvent 
    vLnAct(isW)= LnActSv
    vLnGam(isW)= LnActSv + LOG( One + MolWeitSv *SUM(vMolal(:)) )

    DEALLOCATE(vMolal)
    DEALLOCATE(vLnGamSolute)
    DEALLOCATE(vZSp)
    DEALLOCATE(vA0)

  END SUBROUTINE Solmodel_Calc_Physics

ENDMODULE M_Solmodel_Calc_System
