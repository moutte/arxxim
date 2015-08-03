MODULE M_Dynam_Jacobian
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_, Stop_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dynam_Jacobian
  !
CONTAINS

SUBROUTINE BuildVecs( & !
& LogForAqu,          & !IN
& LogForMin,          & !IN
& nF,nEqA,nMkA,       & !IN
& vX,                 & !IN
& vXf,vLnXf,vXeq,vXmk)
  !
  USE M_Safe_Functions
  USE M_Numeric_Const,ONLY: MinExpDP,MaxExpDP,Ln10
  !
  LOGICAL, INTENT(IN) :: LogForAqu
  LOGICAL, INTENT(IN) :: LogForMin
  INTEGER, INTENT(IN) :: nF,nEqA,nMkA
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vXf(:),vLnXf(:)
  REAL(dp),INTENT(OUT):: vXeq(:),vXmk(:)
  !
  INTEGER:: N
  !
  IF(LogForAqu) THEN
    vXf(1:nF)= FSafe_vExp(vX(1:nF)) !-> avoid overflow
    vLnXf(:)= vX(1:nF)
  ELSE
    vXf(1:nF)= vX(1:nF)
    vLnXf(:)= FSafe_vLog(vX(1:nF))
  ENDIF
  !
  IF(nEqA>0) THEN
    vXeq=Zero
    N= nF
    IF(LogForMin) THEN
      != works with log(mineral amount)
      vXeq(1:nEqA)= FSafe_vExp(vX(N+1:N+nEqA))
    ELSE
      != works directly with mineral amounts
      vXeq(1:nEqA)=vX(N+1:N+nEqA)
    ENDIF
  ENDIF
  !
  IF(nMkA>0) THEN
    vXmk= Zero
    N= nF+nEqA
    IF(LogForMin) THEN
      != works with log(mineral amount)
      vXmk(1:nMkA)= FSafe_vExp(vX(N+1:N+nMkA))
    ELSE
      != works directly with mineral amounts
      vXmk(1:nMkA)=vX(N+1:N+nMkA)
    ENDIF
  ENDIF
  !
END SUBROUTINE BuildVecs

SUBROUTINE Dynam_Jacobian(vX,tJac)
!--
!-- tJac(1:n,1:n)- Jacobian of Dynam_Residual(1:n)
!--
  USE M_Global_Vars,ONLY: vKinFas
  USE M_Basis_Vars, ONLY: nCi,nAs,MWsv
  USE M_Basis_Vars, ONLY: tAlfAs,tNuAs,tAlfPr
  !
  USE M_Dynam_Solve_Vars
  USE M_Dynam_Vars, ONLY: &
  & tAlfKin,tNu_Kin, &
  & DirectSub,  &
  & vMolSec,      &
  & LogForMin,LogForAqu, &
  & vMolarVol,  &
  & vLKinActiv,vLEquActiv, &
  & vSurfK,vVmQsK,vVmAct, &
  & UDarcy,dTime,VBox,dX,PhiF,FOut, &
  & DebNewt,DebJacob, &
  & CoupledCoores
  !
  REAL(dp),DIMENSION(:),                INTENT(IN) :: vX
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)),INTENT(OUT):: tJac
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vXf,vLnXf,vXmk,vXeq,vSumNuAs
  !
  REAL(dp):: dum(SIZE(vX))
  REAL(dp):: Coeff_PhiF,TotFF
  INTEGER :: iPr,jPr,iAs,jAs,iMk,jMk,iMkA,jMkA,I,J
  INTEGER :: nF,nMk,nMkA,nEqA,N
  !
  !nAi= nCi +nAs
  !nF= nAi
  !
  nMk=  SIZE(vKinFas)
  nMkA= COUNT(vLKinActiv)
  nEqA= COUNT(vLEquActiv)
  !
  IF(nMkA +nEqA >0) Coeff_PhiF= UDarcy*dTime/dX /PhiF/PhiF /VBox
  !
  IF(DirectSub) THEN  ;  nF= nCi
  ELSE                ;  nF= nCi +nAs
  ENDIF
  !
  ALLOCATE(vXf(nF))
  ALLOCATE(vLnXf(nF))
  ALLOCATE(vXmk(nMkA))
  ALLOCATE(vXeq(nEqA))
  ALLOCATE(vSumNuAs(nAs))
  !
  DO iAs=1,nAs
    vSumNuAs(iAs)= SUM(tNuAs(iAs,2:nCi))
  ENDDO
  !
  IF(nMkA>0) THEN
    !iMkA=0
    DO iMk=1,nMk
      IF(vLKinActiv(iMk)) THEN
        !iMkA= iMkA + 1
        !-------------------------------------- derivative of VM-Surf*VmA*VmQ --
        dVMdLnXf(iMk,1:nCi)= &
        & vSurfK(iMk) *   vVmAct(iMk)*dVmQdLnXf(iMk,1:nCi) ! &
        ! &    + vVmQsK(iMk)*dVmAdLnXf(iMk,1:nF) )
        !
        dVMdLnXm(iMk,1:nMk)= &
        & vVmQsK(iMk)*vVmAct(iMk) * dSRdLnXm(iMk,1:nMk)
        !
        dVMdXm(iMk,1:nMk)= &
        & vVmQsK(iMk)*vVmAct(iMk) * dSRdXm(iMk,1:nMk)
        !-------------------------------------/ derivative of VM-Surf*VmA*VmQ --
      ELSE !not(vLKinActiv(iMk))
        dVMdLnXf(iMk,1:nCi)= Zero
        dVMdLnXm(iMk,1:nMk)= Zero
        dVMdXm  (iMk,1:nMk)= Zero
      ENDIF
    ENDDO
  ENDIF
  !
  !-----------------------------------------------------------------------------
  CALL BuildVecs( & !
  & LogForAqu,LogForMin,nF,nEqA,nMkA,vX, & !
  & vXf,vLnXf,vXeq,vXmk)
  !
  ! CALL BuildVecs(LogForMin,nF,nEqA,nMkA,vX,vXf,vXeq,vXmk)
  !-----------------------------------------------------------------------------
  !
  tJac=Zero
  !
  !-----------------------------------------------------------------------------
  !derive Dynam_Residual(iPr)= & !variable terms ONLY
  !& (One + UDarcy*dTime/dX/ PhiF) * &
  !& ( DOT_PRODUCT(tAlfPr(iPr,1:nCp),vY(1:nCp))   & !# moles in fluid, from prim'species
  !& + DOT_PRODUCT(tAlfAs(iPr,1:nAs),vY(nCp+1:nCp+nAs))) & !# moles in fluid
  !- vTotF(iPr) & !mole nrs in fluid at previous step
  !- vTotInj(iPr)*UDarcy*dTime/dX &
  !+ dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),vVm(1:nMk))
  !
  !! PhiF= One - DOT_PRODUCT(vY(nF+1:nF+nMk),vMolarVol(1:nMk)) /VBox  !porosity at time t
  !
  !----------------------------------------------------------------- Jacobian --
  !
  !------------------------------------------------- rows 1:nCi: Prim'Species --
  DO iPr=1,nCi
    !
    !----------------------------------------------- 1:nCi-- Aqu'Prim'Species --
    DO jPr=1,nCi
      !
      !! NOTE: Fout=  UDarcy *Vbox /dX /PhiF !! m3.s-1 flow related to porous volume
      !
      IF(LogForAqu) THEN
        tJac(iPr,jPr)= &
        & (One + dTime *FOut /VBox) *tAlfPr(iPr,jPr) *vXf(jPr)
        !! + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdLnXf(1:nMk,jPr))
      ELSE
        tJac(iPr,jPr)= &
        & (One + dTime *FOut /VBox) *tAlfPr(iPr,jPr)
      ENDIF
      !
      IF(nMkA>0) THEN
        DO iMk=1,nMk
          IF(vLKinActiv(iMk)) THEN
            tJac(iPr,jPr)= tJac(iPr,jPr) &
            &            + dTime *tAlfKin(iPr,iMk) *dVMdLnXf(iMk,jPr)
          ENDIF
        ENDDO
      ENDIF
      !
    ENDDO
    !----------------------------------------------/ 1:nCi-- Aqu'Prim'Species --
    !
    IF(DirectSub) THEN
      DO jAs=1,nAs
        tJac(iPr,1)= tJac(iPr,1) &
        &          + (One + dTime *FOut /VBox) &
        &            *vMolSec(jAs) *tAlfAs(iPr,jAs) *(One - vSumNuAs(jAs))
      ENDDO
      DO jPr=2,nCi
        DO jAs=1,nAs
          tJac(iPr,jPr)= tJac(iPr,jPr) &
          &            + (One + dTime *FOut /VBox) &
          &              *vMolSec(jAs) *tAlfAs(iPr,jAs) *tNuAs(jAs,jPr)
        ENDDO
      ENDDO
    ELSE
    !--------------------------------------------- 1:nCi-- Aqu'Second'species --
      DO jAs=1,nAs
        ! same as prim'sp,
        ! but we write again because tAlfa is split into tAlfPr U tAlfAs
        IF(LogForAqu) THEN
          tJac(iPr,nCi+jAs)= &
          & (One + dTime *FOut /VBox) *vXf(nCi+jAs) *tAlfAs(iPr,jAs)
        ELSE
          tJac(iPr,nCi+jAs)= &
          & (One + dTime *FOut /VBox) *tAlfAs(iPr,jAs)
        ENDIF
        !! + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdLnXf(1:nMk,jAs))
        !<removed: VM does not depend on second'species>
        !IF(nMkA>0) THEN
        !  DO iMk=1,nMk
        !    IF(vLKinActiv(iMk)) THEN
        !      tJac(iPr,jAs)= tJac(iPr,jAs) &
        !      &            + dTime *tAlfKin(iPr,iMk) *dVMdLnXf(iMk,jAs)
        !    ENDIF
        !  ENDDO
        !ENDIF
        !</removed>
      ENDDO
    ENDIF
    !----------------------------------------------------/ Aqu'Second'species --
    !
    IF(nMkA +nEqA >0) THEN
      IF(DirectSub) THEN
        TotFF= DOT_PRODUCT(tAlfPr(iPr,1:nCi),vXf(1    :nCi)) &
        &    + DOT_PRODUCT(tAlfAs(iPr,1:nAs),vMolSec(1:nAs))
      ELSE
        TotFF= DOT_PRODUCT(tAlfPr(iPr,1:nCi),vXf(1    :nCi)) &
        &    + DOT_PRODUCT(tAlfAs(iPr,1:nAs),vXf(nCi+1:nCi+nAs))
      ENDIF
    ENDIF
    !
    !-------------------------------------------------- 1:nCi-- equil' phases --
    IF(nEqA>0) THEN
      !
      N= nF
      !
      J= 0
      DO iMk=1,nMk
        IF(vLEquActiv(iMk)) THEN
          J=J+1
          !! vFunc(I)= vFunc(I) + tAlfKin(iPr,iMk)*vXeq(J)
          !! CAVEAT: option LogForMin not implemented !!
          IF(CoupledCoores) THEN
            tJac(iPr,N+J)= tAlfKin(iPr,iMk)
          ELSE
            !dPhi/dXeq= - vMolarVol(iMk) /VBox
            !dPsi/dXeq= TotFF*(UDarcy*dTime/dX) *(-1/PhiF/PhiF) *dPhi/dXeq
            tJac(iPr,N+J)= Coeff_PhiF *TotFF *vMolarVol(iMk) + tAlfKin(iPr,iMk)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    !---------------------------------------------------------/ equil' phases --
    !
    !------------------------------------------------ 1:nCi -- kinetic phases --
    IF(nMkA>0) THEN
      !
      N= nF+nEqA
      !
      J=0
      DO iMk=1,nMk
        !
        IF(vLKinActiv(iMk)) THEN
          !
          J=J+1
          IF(LogForMin) THEN
            !
            IF(CoupledCoores) THEN
              tJac(iPr,N+J)= &
              & dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
            ELSE !-> Phi IMPLICIT
              tJac(iPr,N+J)= &
              & Coeff_PhiF *TotFF * vXmk(J) *vMolarVol(iMk) &
              + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
              !**! + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
            ENDIF
            !
          ELSE
            !
            IF(CoupledCoores) THEN
              tJac(iPr,N+J)= &
              & dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
              !**! & dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
            ELSE !-> Phi IMPLICIT
              tJac(iPr,N+J)= &
              & Coeff_PhiF *TotFF *vMolarVol(iMk) &
              + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
              !**! + dTime * DOT_PRODUCT(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
            ENDIF
            !
          ENDIF
          !
        ENDIF !IF(vLKinActiv(iMk))
        !
      ENDDO
    ENDIF
    !--------------------------------------------------------/ kinetic phases --
  ENDDO
  !-----------------------------------------------/ rows 1:nCi (Prim'Species) --
  !
  !----------------------------------- rows nCi+1:nCi+nAs: Aqu'Second'Species --
  !
  !derive Psi(nCp+iAs)=
  ! DOT_PRODUCT(tNuAs(iAs,1:nCp),vX(1:nCp)) - vX(nCp+iAs) &
  !- (SUM(tNuAs(iAs,1:nCp)) - One) *(vX(isW) + LOG(MWSv)) + Constant Terms
  !
  IF(.NOT. DirectSub) THEN
    !
    IF(LogForAqu) THEN
      !------------------------------------------------------------ LogForAqu --
      !-------------------------------- deriv' second'species vs prim'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,1)=     One - vSumNuAs(iAs) !Solvent
        DO iPr=2,nCi
          tJac(nCi+iAs,iPr)= tNuAs(iAs,iPr)
        ENDDO
      ENDDO
      !---/
      !------------------------------ deriv' second'species vs second'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= - One !-> matrix block = -(Identity)
      ENDDO
      !---/
      !------------------------------------------------------------/LogForAqu --
      !
    ELSE
      !
      !-------------------------------------------------------- NOT LogForAqu --
      !-------------------------------- deriv' second'species vs prim'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,1)= vMolSec(iAs) *(One -vSumNuAs(iAs)) /vX(1)
        DO iPr=2,nCi
          IF(tNuAs(iAs,iPr)/= 0.0D0) &
          & tJac(nCi+iAs,iPr)= vMolSec(iAs) *tNuAs(iAs,iPr) /vX(iPr)
        END DO
      END DO
      !---/
      !------------------------------ deriv' second'species vs second'species --
      DO iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= - One !-> matrix block = -(Identity)
      ENDDO
      !---/
      !--------------------------------------------------------/NOT LogForAqu --
    ENDIF
  ENDIF
  !
  !------------------------------------------ rows nF+1:nF+nEqA: equil'phases --
  IF(nEqA>0) THEN
    !
    N= nF
    !
    J= 0
    DO iMk=1,nMk
      IF(vLEquActiv(iMk)) THEN
        J= J +1
        tJac(N +J,1)= - SUM(tNu_Kin(iMk,2:nCi))
        tJac(N +J,2:nCi)=   tNu_Kin(iMk,2:nCi)
      ENDIF
    ENDDO
  ENDIF
  !-----------------------------------------/ rows nF+1:nF+nEqA: equil'phases --
  !
  !----------------------------------- rows nF+nEqA+1:nF+nEqA+nMk: kin'phases --
  !------ Dynam_Residual(nF+i)- X(nF+i)) -Xi(nF+i) - dTime*vSurfK(i)*vVmF_(i) --
  !
  IF(nMkA>0) THEN
    !
    N= nF +nEqA
    !
    iMkA=0
    DO iMk=1,nMk !!!CHECK!!!
      !
      IF(vLKinActiv(iMk)) THEN
        !
        iMkA=iMkA+1
        !
        tJac(N+iMkA,1:nCi)= -dTime * dVMdLnXf(iMk,1:nCi)
        !
        jMkA=0
        DO jMk=1,nMk
          IF(vLKinActiv(jMk)) THEN
            !
            jMkA=jMkA+1
            !
            IF(jMkA==iMkA) THEN
              IF(LogForMin) THEN
                tJac(N+iMkA,N+jMkA)= vXmk(iMkA) - dTime * dVMdLnXm(iMk,jMk)
              ELSE
                tJac(N+iMkA,N+jMkA)= One - dTime * dVMdXm(iMk,jMk)
              ENDIF
            ELSE
              IF(LogForMin) THEN
                tJac(N+iMkA,N+jMkA)= - dTime * dVMdLnXm(iMk,jMk)
              ELSE
                tJac(N+iMkA,N+jMkA)= - dTime * dVMdXm(iMk,jMk)
              ENDIF
            ENDIF
            !
          ENDIF
        ENDDO
      ENDIF
      !
    ENDDO
  ENDIF !IF(nMkA>0)
  !----------------------------------/ rows nF+nEqA+1:nF+nEqA+nMk: kin'phases --
  !
  DEALLOCATE(vXf)
  DEALLOCATE(vLnXf)
  DEALLOCATE(vXmk)
  DEALLOCATE(vXeq)
  DEALLOCATE(vSumNuAs)
  !
  !-----------------------------------------------------------/ DynamJacobian --
  !
  !DebNewt= (iDebug>2)
  IF(DebNewt .AND. DebJacob) THEN !
    !
    !
    WRITE(31,'(A)') "Dynamic Jacobian"
    !
    !--- print max value of each column
    dum=MAXVAL(ABS(tJac),DIM=2)
    DO I=1,nF+nEqA+nMkA; WRITE(31,'(E8.1,A1)',ADVANCE='NO') dum(I), T_; ENDDO
    WRITE(31,'(A)') "MAXVAL(COL)"
    !---/
    !
    DO I=1,nF+nEqA+nMkA
      DO J=1,nF+nEqA+nMkA
        IF(tJac(I,J)/=Zero) THEN; WRITE(31,'(F7.1,A1)',ADVANCE='NO') tJac(I,J), T_
        ELSE;                     WRITE(31,'(A1,  A1)',ADVANCE='NO') ".",       T_
        ENDIF
      ENDDO
      WRITE(31,'(I3)') I
    ENDDO
    !
    DebJacob=.FALSE.
    !
  ENDIF
  !
ENDSUBROUTINE Dynam_Jacobian

ENDMODULE M_Dynam_Jacobian
