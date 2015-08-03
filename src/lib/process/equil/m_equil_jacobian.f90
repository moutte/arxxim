MODULE M_Equil_Jacobian
  USE M_Trace,        ONLY: iDebug,fTrc,T_,Stop_,Pause_
  USE M_Numeric_Const,ONLY: MinExpDP,MaxExpDP
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Jacobian
  !
CONTAINS

SUBROUTINE Equil_Jacobian(vX,tJac)
!--
!-- EquJacobian of the system EquResidual; used by Newton
!-- NB: in vX and in tJac, solvent has index 1 !!!!
!--
  USE M_Safe_Functions
  !
  USE M_Basis_Vars,ONLY: isW,MWsv
  USE M_Basis_Vars,ONLY: tAlfPr,tAlfAs,tNuAs,tAlfFs,tNuFas
  USE M_Basis_Vars,ONLY: vOrdPr,vOrdAq,nAx,nAs,nCi
  !
  USE M_Equil_Vars,ONLY: vMolal,vLnAct,vLnGam,vMolSec
  !~ USE M_Equil_Vars,ONLY: vFasYes
  USE M_Equil_Vars,ONLY: cEquMode
  USE M_Equil_Vars,ONLY: vAffScale,vFasAff,dXi
  !~ USE M_Equil_Vars,ONLY: OsmoSv,UseOsmoSv
  USE M_Equil_Vars,ONLY: DebNewt,DebJacob
  USE M_Equil_Vars,ONLY: LogForAqu,DirectSub
  USE M_Equil_Vars,ONLY: tAlfEq,tNuEq,nEquFas
  !
  REAL(dp),DIMENSION(:),                INTENT(IN) :: vX
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)),INTENT(OUT):: tJac
  !
  REAL(dp),DIMENSION(1:nCi):: dRatdX
  INTEGER :: iCi,jCi,iAs,iAx,I
  REAL(dp):: X
  REAL(dp):: AffScale
  ! REAL(dp):: SumXi !sum of solute mole nrs, used when working with osmo'coeff
  REAL(dp),ALLOCATABLE:: vLnXf(:),vXf(:)
  REAL(dp),ALLOCATABLE:: vSumNuAs(:)
  !
  INTEGER:: nF !nr of aqu'species in array vX of residual
  INTEGER:: nM !nr of non-aqu'species in array vX of residual
  !
  tJac=  Zero
  !
  !IF(DirectSub) CALL Stop_("No Analytic Jacobian ...")
  !
  IF(DirectSub) THEN  ;   nF= nCi
  ELSE                ;   nF= nCi +nAs
  ENDIF
  !
  nM= nEquFas
  !
  ALLOCATE(vXf(nF))  ;  vXf= Zero
  !
  ALLOCATE(vLnXf(nF))
  IF(LogForAqu) THEN
    vLnXf(1:nF)= vX(1:nF)
    vXf(:)= FSafe_vExp(vLnXf(:)) ! avoid overflow with EXP
  ELSE
    vXf(1:nF)= vX(1:nF)
    vLnXf(1:nF)= FSafe_vLog(vX(1:nF))
  ENDIF
  !
  ALLOCATE(vSumNuAs(nAs))
  DO iAs=1,nAs
    X= Zero
    DO I=1,nCi
      IF(I/=isW) X= X+tNuAs(iAs,I)
    ENDDO
    vSumNuAs(iAs)= X
  ENDDO
  !
  !DO I=1,nAx
  !  vMolal(vOrdPr(nCi+I))= EXP(vLnAct(vOrdPr(nCi+I))-vLnGam(vOrdPr(nCi+I)))
  !ENDDO
  !IF(nAx>0) vMolal(vOrdPr(nCi+1:nCi+nAx))= &
  !& EXP(vLnAct(vOrdPr(nCi+1:nCi+nAx))-vLnGam(vOrdPr(nCi+1:nCi+nAx))) !-> i.e. molality
  !
  IF(LogForAqu) THEN
    !------------------------------------------ material conservation --
    !--------------------------- rows 1:nCi - "internal" prim'species --
    DO iCi=1,nCi
      !--- solvent
      tJac(iCi,1)= tAlfPr(iCi,1)*vXf(1)
      !-------- prim'mobile'aqu'species -> add corresponding moles of Solvent --
      !~ IF(nAx>0) &
      !~ & tJac(iCi,1)= tJac(iCi,1) &
      !~ & + vXf(isW) *MWSv &
      !~ &   *DOT_PRODUCT(tAlfPr(iCi,nCi+1:nCi+nAx),vMolal(vOrdPr(nCi+1:nCi+nAx)))
      DO iAx=1,nAx
        tJac(iCi,1)= tJac(iCi,1) &
        &          + vXf(isW) *MWSv *tAlfPr(iCi,nCi+iAx) *vMolal(vOrdPr(nCi+iAx))
      ENDDO
      !--- prim'intern'species
      DO jCi=2,nCi
        tJac(iCi,jCi)= tAlfPr(iCi,jCi)*vXf(jCi)
      ENDDO
      !
      !--- second'species
      IF(DirectSub) THEN
        DO iAs=1,nAs
          tJac(iCi,1)= tJac(iCi,1) &
          &          + vMolSec(iAs) *tAlfAs(iCi,iAs) *(One - vSumNuAs(iAs))
        ENDDO
        DO jCi=2,nCi
          DO iAs=1,nAs
            tJac(iCi,jCi)= tJac(iCi,jCi) &
            &            + vMolSec(iAs) *tAlfAs(iCi,iAs) *tNuAs(iAs,jCi)
          ENDDO
        ENDDO
      ELSE
        DO iAs=1,nAs
          tJac(iCi,nCi+iAs)= vXf(nCi+iAs) *tAlfAs(iCi,iAs)
        ENDDO
      ENDIF
      !
    ENDDO
    !---/"internal" prim'species --
    !
    !~ IF(bMolal) THEN
      !~ tJac(isW,:)= Zero
      !~ tJac(isW,1)= vXf(isW)
    !~ ENDIF
    !
    IF(DirectSub) THEN
      !
    ELSE
      !--- rows nCi+1:nAq- sec'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,1)= One-SUM(tNuAs(iAs,2:nCi)) !for Solvent
        tJac(nCi+iAs,2:nCi)=     tNuAs(iAs,2:nCi)  !for other aqu.species
      ENDDO
      DO iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= -One
      ENDDO
      !---/
    ENDIF
    !
    !~ IF(UseOsmoSv) THEN
      !~ !
      !~ SumXi= SUM(vXf(1:nAs))-vXf(isW)
      !~ !
      !~ DO iAs=1,nAs !add terms for osmotic coeff'
        !~ tJac(nCi+iAs,1  )= tJac(nCi+iAs,1) &
        !~ &                + tNuAs(iAs,isW) *OsmoSv *SumXi /vXf(isW) != d/dXw(LnActW)
        !~ DO I=2,nAs !for other aqu.species
          !~ tJac(nCi+iAs,I)= tJac(nCi+iAs,I) &
          !~ &              - tNuAs(iAs,isW) *OsmoSv *vXf(I)/vXf(isW) != d/dXi(LnActW)
        !~ ENDDO
      !~ ENDDO
      !~ !
    !~ ENDIF

  ELSE
    !-------------------------------------------------- NOT LogForAqu --
    !~ IF(DirectSub) CALL Stop_("No Analytic Jacobian ...")
    !
    !------------------------------------------ material conservation --
    DO iCi=1,nCi
      !
      !--- solvent --
      tJac(iCi,1)= tAlfPr(iCi,1)
      DO iAx=1,nAx
        tJac(iCi,1)= tJac(iCi,1) &
        &          + tAlfPr(iCi,nCi+iAx)*vMolal(vOrdPr(nCi+iAx))*MWsv
      ENDDO
      !---/
      !--- prim'intern'species
      tJac(iCi,2:nCi)= tAlfPr(iCi,2:nCi)
      !
      IF(DirectSub) THEN
        !--- solvent --
        DO iAs=1,nAs
          IF(tNuAs(iAs,iCi)/= 0.0D0) THEN
            ! RESIDUAL ! X= X + tAlfAs(iCi,iAs) *vMolal(iAs)*(vX(isW)*MWSv)
            tJac(iCi,1)= tJac(iCi,1) &
            &          + tAlfAs(iCi,iAs) *vMolal(iAs)*MWSv
          ENDIF
        END DO
      ELSE
        !--- second'species
        tJac(iCi,nCi+1:nCi+nAs)= tAlfAs(iCi,1:nAs)
      ENDIF
      !
    END DO
    !-----------------------------------------/ material conservation --
    !
    IF(.NOT. DirectSub) THEN
      !-------------------------------------- sec'species equilibrium --
      !------------------------ deriv' second'species vs prim'species --
      DO iAs=1,nAs
        !~ tJac(nCi+iAs,1)= vMolal(iAs)*(vX(isW)*MWSv) *(One -vSumNuAs(iAs)) /vX(isW)
        tJac(nCi+iAs,1)= vMolSec(iAs) *(One -vSumNuAs(iAs)) /vX(isW)
        DO iCi=1,nCi
          IF(iCi/=isW .AND. tNuAs(iAs,iCi)/= 0.0D0) &
          !~ & tJac(nCi+iAs,iCi)= vMolal(iAs)*(vX(isW)*MWSv) *tNuAs(iAs,iCi) /vX(iCi)
          & tJac(nCi+iAs,iCi)= vMolSec(iAs) *tNuAs(iAs,iCi) /vX(iCi)
        END DO
      END DO
      !---------------------- deriv' second'species vs second'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= -One
      ENDDO
      !-------------------------------------/ sec'species equilibrium --
    ENDIF
    !
    !-------------------------------------------------/ NOT LogForAqu --
  ENDIF
  !
  IF(nEquFas>0) THEN

    SELECT CASE(cEquMode)

    CASE("EQ2")
      DO I=1,nEquFas
        DO iCi=1,nCi
          tJac(iCi,nF+I)= tAlfEq(iCi,I)
        ENDDO
        IF(LogForAqu) THEN
          tJac(nF+I,isW)=   SUM(tNuEq(I,2:nCi))
          tJac(nF+I,2:nCi)= - tNuEq(I,2:nCi)
        ELSE
          tJac(nF+I,isW)=   SUM(tNuEq(I,2:nCi)) / vX(isW)
          tJac(nF+I,2:nCi)= - tNuEq(I,2:nCi) / vX(2:nCi)
        ENDIF
      ENDDO

    CASE("EQ1")
      DO I=1,nEquFas
        AffScale= SUM(ABS(tNuEq(I,:)))
        dRatdX(isW)= - SUM(tNuEq(I,2:nCi)) &
        &            *FSafe_Exp(- vFasAff(I)/AffScale) &
        &            /AffScale
        dRatdX(2:nCi)= tNuEq(I,2:nCi) &
        &            *FSafe_Exp(- vFasAff(I)/AffScale) &
        &            /AffScale
        DO iCi=1,nCi
          tJac(iCi,1:nCi)= tJac(iCi,1:nCi) &
          &              + tAlfEq(iCi,I) *dXI *dRatdX(1:nCi)
        ENDDO
        !
        tJac(nF+I,1:nCi)= - dXi *dRatdX(1:nCi)
        tJac(nF+I,nF+I)=    One
      ENDDO

    ENDSELECT

  ENDIF
  !
  IF(DebNewt .AND. DebJacob) THEN
    CALL ShoJacob(nCi,nF,nAs,tJac)
    DebJacob=.FALSE.
  ENDIF
  !
  DEALLOCATE(vSumNuAs)
  DEALLOCATE(vXf)
  DEALLOCATE(vLnXf)
  !
ENDSUBROUTINE Equil_Jacobian

SUBROUTINE ShoJacob(nCi,nF,nAs,tJac)
  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirLog
  USE M_Global_Vars, ONLY: vSpc
  !
  INTEGER, INTENT(IN):: nCi,nAs,nF
  REAL(dp),INTENT(IN):: tJac(:,:)
  !
  INTEGER::I,J,f
  !
  WRITE(fTrc,'(/,A)') "< ShoJacob"
  WRITE(fTrc,'(A)') "Isaac Newton_&_Carl_Gustav_Jacob_Jacobi -> FILE=JACOB.LOG"
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirLog)//"debug_jacob.log")
  DO I=1,nCi; WRITE(f,'(A7,A1)',ADVANCE='NO') vSpc(I)%NamSp,    T_; ENDDO
  DO I=1,nAs; WRITE(f,'(A7,A1)',ADVANCE='NO') vSpc(nCi+I)%NamSp,T_; ENDDO
  WRITE(f,*)
  DO I=1,nF
    DO J=1,nF
      IF(tJac(I,J)/=Zero) THEN; WRITE(f,'(E8.1,A1)',ADVANCE='NO') tJac(I,J), T_
      ELSE;                     WRITE(f,'(A1,A1)',  ADVANCE='NO') ".",       T_
      ENDIF
    ENDDO
    WRITE(f,'(I3)') I
  ENDDO
  WRITE(f,'(A7)') "_______"
  DO I=1,nF
    DO J=1,nF
      IF(tJac(I,J)/=Zero) THEN; WRITE(f,'(G15.4,A1)',ADVANCE='NO') tJac(I,J), T_
      ELSE                    ; WRITE(f,'(A1,A1)',   ADVANCE='NO') "0",       T_
      ENDIF
    ENDDO
    WRITE(f,'(I3)') I
  ENDDO
  CLOSE(f)
  WRITE(fTrc,'(A,/)') "</ ShoJacob"
ENDSUBROUTINE ShoJacob

ENDMODULE M_Equil_Jacobian


