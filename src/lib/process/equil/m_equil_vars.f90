MODULE M_Equil_Vars
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  USE M_Kinds
  USE M_T_MixModel,ONLY: MaxPole
  IMPLICIT NONE
  !
  PUBLIC
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: &
  & vTotS,    & !1:nCp=  total amounts of COMPONENTS (incl. charge balance) in the SYSTEM at time T
  ! in SPC, vTotS = vTot_inFluid
  ! in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
  & vMolF,    & !1:nAq=  nr moles of aqu'species in the system, e.g. for 1 kg H2O
  & vMolal,   & !1:nAq=  molalities
  & vMolSec,  & !1:nAs= mole numbers second species, for DirectSub
  & vFasMole, & !1:nFs= mole nrs of phases (added for SpeciaEquil) -> SAVEd to vFas(:)%Mole
  & vLnAct,   & !
  & vLnGam,   & !
  & vDeltaG_As, & ! DeltaGibbsFreeEnergy of formation reaction, for Second'Aqu'Species
  & vDGapp_As,  & ! "apparent deltaG": includes delta(log(gammas solute species))
  & vDeltaG_Fs    ! DeltaGibbsFreeEnergy of formation reaction, for phases of fixed composition
  !
  LOGICAL, DIMENSION(:),ALLOCATABLE:: vYesList != the phase is taken into account
  !~ & vFasYes     !vFasYes= the phase is "active"
  !
  ! variables added for Specia_CalcEquil_N
  ! affinity of minerals or gases, scaling factor, rate (for cEquMode=3)
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vFasAff,vFasRat,vAffScale 
  !
  !--- added for M_Equil_3 / M_Equil_3_Mix
  REAL(dp),ALLOCATABLE:: &
  ! for equilibrium non'aqu'phases,
  & vDeltaG_Eq(:), &  ! DeltaGibbsFreeEnergy of formation reaction
  & tAlfEq(:,:),   &  !
  & tNuEq(:,:)        !
  ! Relation betweeen (vDeltaG_Eq,tAlfEq,tNuEq)
  !               and (vDeltaG_Fs,tAlfFs,tNuFas)
  !~ DO iFs=1,SIZE(vFasYes)
    !~ IF(vFasYes(iFs)) THEN
      !~ I=I+1
      !~ vDeltaG_Eq(I)= vDeltaG_Fs(iFs)
      !~ tAlfEq(:,I)= tAlfFs(:,iFs)
      !~ tNuEq(I,:)= tNuFas(iFs,:)
    !~ ENDIF
  !~ ENDDO
  !
  TYPE:: T_EquPhase
    CHARACTER(LEN=23):: NamEq
    INTEGER :: iPur
    INTEGER :: iMix
    INTEGER :: nPole
    INTEGER :: vIPole(1:MaxPole)
    REAL(dp):: vXPole(1:MaxPole)
    REAL(dp):: vLnAct(1:MaxPole)
    !REAL(dp):: Grt
    REAL(dp):: Mole
  END TYPE T_EquPhase
  !
  TYPE(T_EquPhase),ALLOCATABLE:: vEquFas(:)
  INTEGER :: nEquFas
  !---/
  !
  !~ REAL(dp):: OsmoSv
  REAL(dp):: dXi
  !
  LOGICAL,DIMENSION(:),ALLOCATABLE:: vTooLow !,vTooHih
  !
  INTEGER:: nEvalFunc
  !
  LOGICAL:: LogForAqu= .TRUE.
  !
  LOGICAL:: Complementarity= .FALSE.  !! .TRUE.  !!
  !
  LOGICAL:: DirectSub= .FALSE.
  ! DirectSub is "direct substitution":
  ! second'species mole numbers substituted in material consevation equations
  !
  CHARACTER(LEN=30):: cMethod= "NEWTONPRESS"
  !
  INTEGER:: & !
  & fTrAct=  0, & !"trace" file for activity calculations, OPENed IF DebActiv
  & fActiz=  0, & !log file for activities
  & fTrcEq=  0, & !
  & iNewt=   0, & !
  & iCountEq=0
  !
  LOGICAL:: & 
  & DebNewt=   .FALSE., &
  & DebJacob=  .FALSE., &
  & bFinDIF=   .FALSE., &
  & DebActiv=  .FALSE., &
  & DebFormula=.FALSE.
  !
  ! variables for Equil_EqN
  CHARACTER(LEN=3):: cEquMode= "EQ2"
  ! flag for the TYPE of "equilibrium speciation"
  !   EQ1= step method, titration-like
  !   EQ2= direct method, suitable for non-singular phase assemblage, given a priori
  !
  !~ LOGICAL:: UseOsmoSv= .FALSE. !.TRUE. !
  ! the system (residual=zero) is solved at constant activity coeff's of solute
  ! currently, concerning the solvent, the system can be solved
  !   (1) at constant ln(Activ(solvent))  (= was the only option available until 2008)
  !or (2) at constant osmotic coefficient (= new option, under test)
  !
  !~ LOGICAL:: bMolal=  .FALSE. !.TRUE. !
  !
  LOGICAL:: bDirect= .FALSE. !.TRUE. !
  !
  !--- parameters for convergence on gammas --
  INTEGER, PARAMETER:: GamMaxIts=100    !maximum number of iterations on gamma's
  REAL(dp),PARAMETER:: TolGam=   1.0D-6 !convergence condition on delta(Log(Gamma))
  REAL(dp),PARAMETER:: TolAct=   1.0D-6 !convergence condition on delta(Log(Activ))
  !---/
  !
  !--- parameters for Newton convergence --
  INTEGER, PARAMETER :: NewtMaxIts=200 !maximum number of iterations !!!!!!200
  REAL(dp),PARAMETER :: &
  !& NewtTolF=   1.0E-03_dp,& !convergence criterion on function values
  !& NewtTolMin= 1.0E-06_dp,& !criterion for spurious convergence
  & NewtTolF=   1.0E-09_dp, &  !convergence criterion on function values
  & NewtTolMin= 1.0E-09_dp, &  !criterion for spurious convergence
  !!& NewtTolX=   1.0E-12_dp
  & NewtTolX=   2.0E-16_dp !EPSILON(REAL(dp))
  ! convergence criterion on dx -> approx. 2.E-16
  ! -> NewtTolX is used to check whether Newton loops on a stationary point
  !---/
  !
CONTAINS

SUBROUTINE Equil_Vars_Alloc(nCp,nSp,nAq,nAs,nFs)
  INTEGER,INTENT(IN):: nCp,nSp,nAq,nAs,nFs
  !
  CALL Equil_Vars_Clean
  !
  ALLOCATE(vTotS(1:nCp))       ;  vTotS= Zero
  !
  ALLOCATE(vTooLow(1:nSp))     ;  vTooLow=.FALSE.
  ALLOCATE(vLnAct(1:nSp))      ;  vLnAct= Zero 
  ALLOCATE(vLnGam(1:nSp))      ;  vLnGam= Zero
  !
  ALLOCATE(vMolF(1:nAq))
  ALLOCATE(vMolal(1:nAq))
  !
  ALLOCATE(vEquFas(1:nCp))
  !
  ALLOCATE(vDeltaG_As(1:nAs))
  ALLOCATE(vDGapp_As(1:nAs))
  !
  !~ ALLOCATE(vFasYes(1:nFs))     ; vFasYes=  .FALSE.
  ALLOCATE(vFasMole(1:nFs))    ; vFasMole= Zero
  ALLOCATE(vYesList(1:nFs))    ; vYesList= .FALSE.
  !
  ALLOCATE(vDeltaG_Fs(1:nFs))  ; vDeltaG_Fs= Zero 
  ALLOCATE(vFasAff(1:nFs))     ; vFasAff=    Zero 
  ALLOCATE(vFasRat(1:nFs))     ; vFasRat=    Zero
  ALLOCATE(vAffScale(1:nFs))   ; vAffScale=  One
  !
ENDSUBROUTINE Equil_Vars_Alloc

SUBROUTINE Equil_Vars_Clean

  IF(ALLOCATED(vTotS))   DEALLOCATE(vTotS)
  IF(ALLOCATED(vMolF))   DEALLOCATE(vMolF)
  IF(ALLOCATED(vMolal))  DEALLOCATE(vMolal)
  IF(ALLOCATED(vLnAct))  DEALLOCATE(vLnAct)
  IF(ALLOCATED(vLnGam))  DEALLOCATE(vLnGam)
  IF(ALLOCATED(vTooLow)) DEALLOCATE(vTooLow)
  !
  IF(ALLOCATED(vEquFas)) DEALLOCATE(vEquFas)
  !
  IF(ALLOCATED(vDeltaG_As)) DEALLOCATE(vDeltaG_As)
  IF(ALLOCATED(vDeltaG_Fs)) DEALLOCATE(vDeltaG_Fs)
  IF(ALLOCATED(vDGapp_As))  DEALLOCATE(vDGapp_As)
  !
  IF(ALLOCATED(vYesList) ) DEALLOCATE(vYesList) 
  !~ IF(ALLOCATED(vFasYes)  ) DEALLOCATE(vFasYes)  
  IF(ALLOCATED(vFasAff)  ) DEALLOCATE(vFasAff)  
  IF(ALLOCATED(vFasRat)  ) DEALLOCATE(vFasRat)  
  IF(ALLOCATED(vFasMole) ) DEALLOCATE(vFasMole)  
  IF(ALLOCATED(vAffScale)) DEALLOCATE(vAffScale)

ENDSUBROUTINE Equil_Vars_Clean

ENDMODULE M_Equil_Vars

