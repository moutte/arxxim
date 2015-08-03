MODULE M_Dynam_Residual
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dynam_Residual
  PUBLIC:: Dynam_Converge
  PUBLIC:: vTolCoef,vNewtTolF
  !
  REAL(dp),ALLOCATABLE:: vTolCoef(:)
  REAL(dp),ALLOCATABLE:: vNewtTolF(:)
  !
  LOGICAL:: UseTolCoef= .FALSE.  !! .TRUE. !!

CONTAINS

LOGICAL FUNCTION Dynam_Converge(vF,vTolF)
!--
!-- makes possible adaptation of convergence criteria to type of condition
!-- e.g. material conservation vs potential, kinetic, etc.
!--
  REAL(dp),INTENT(IN):: vF(:),vTolF(:)
  !
  IF(UseTolCoef) THEN
    Dynam_Converge= ALL(ABS(vF(:)) < vNewtTolF(:) *vTolCoef(:))
  ELSE
    Dynam_Converge= ALL(ABS(vF(:)) < vNewtTolF(:))
  ENDIF
  !
ENDFUNCTION Dynam_Converge

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

FUNCTION Dynam_Residual(vX) RESULT(vFunc)
!--
!-- system to be solved
!--
  USE M_Global_Vars,  ONLY: vKinFas,nAq,vFas
  USE M_T_KinFas,     ONLY: KinFas_CalcSurf
  USE M_KinFas_Surf,  ONLY: KinFas_Surf_Calc
  USE M_KinRate
  !
  USE M_Basis_Vars, ONLY: MWSv,isW,nCi,nAs,nAx,nMx
  USE M_Basis_Vars, ONLY: tAlfPr,tAlfAs,tNuAs
  USE M_Dynam_Vars, ONLY: LogForAqu, LogForMin, DirectSub, CoupledCoores
  USE M_Dynam_Vars, ONLY: PhiF,PhiInert,VBox,UDarcy,dX,dTime
  USE M_Dynam_Vars, ONLY: vQTotInj,vTotInj,vTotF,Fout
  USE M_Dynam_Vars, ONLY: vLnGam,vLnBuf,vLnAct,vMolSec
  USE M_Dynam_Vars, ONLY: vLKinActiv,vLEquActiv
  USE M_Dynam_Vars, ONLY: vKinMod,vMolarVol
  USE M_Dynam_Vars, ONLY: vMolK,vSurfK,vMolK0,vSurfK0
  USE M_Dynam_Vars, ONLY: vVmQsK,vVmAct
  USE M_Dynam_Vars, ONLY: tNu_Kin,tAlfKin,vDG_Kin,vDG_As,vDLGam_As
  USE M_Dynam_Vars, ONLY: nEvalFunc
  !
  USE M_Dynam_Solve_Vars
  !
  REAL(dp),DIMENSION(:),INTENT(IN):: vX
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vXf,vLnXf,vXmk,vXeq,vVm,vMolal,vSumNuAs
  !
  REAL(dp),DIMENSION(1:SIZE(vX)):: vFunc
  !
  REAL(dp):: QTotInjPr,X,TotFF
  INTEGER :: I,J,iPsi,iPr,iAs,iMk,iMkA,nF,nMk,nMkA,nEqA,nCx
  !
  !-----------------------------------------------under construction ...
  IF(nAx>0) THEN
    ALLOCATE(vMolal(nAx))
    !update molalities of aqueous "mobile" species (i.e. buffered)
    !-> used to calculate mole number= molality * solvent'mass
    !! vMolal(1:nAx)= &
    !! & EXP(vLnAct(vOrdPr(nCi+1:nCi+nAx))-vLnGam(vOrdPr(nCi+1:nCi+nAx)))
    !-> molality
  ENDIF
  !----------------------------------------------/under construction ...
  !
  IF(iDebug>2) nEvalFunc= nEvalFunc +1
  !
  nCx=  nAx+nMx
  nMk=  SIZE(vKinFas)
  nMkA= COUNT(vLKinActiv)
  nEqA= COUNT(vLEquActiv)
  !
  IF(DirectSub) THEN  ;  nF= nCi
  ELSE                ;  nF= nAq
  ENDIF
  !
  ALLOCATE(vXf(nF))
  ALLOCATE(vLnXf(nF))
  ALLOCATE(vXmk(nMkA))
  ALLOCATE(vXeq(nEqA))
  ALLOCATE(vVm(nMk))
  !
  CALL BuildVecs( & !
  & LogForAqu,LogForMin,nF,nEqA,nMkA,vX, & !
  & vXf,vLnXf,vXeq,vXmk)
  !
  CALL KinRate_Calc( & !
  & dTime,                 & !IN
  & vLKinActiv,            & !IN
  & vLnXf,vXmk,            & !IN
  & vSurfK,vVmQsK,vVmAct,  & !INOUT
  & dSRdXm,dSRdLnXm,       & !OUT
  & dVmQdLnXf,             & !OUT
  & dVMdLnXf,dVMdLnXm)       !OUT
  !
  vVm(:)= vSurfK(:) *vVmAct(:) *vVmQsK(:)
  !
  CALL PhiF_Calc( & !
  & PhiInert,              & !IN
  & vLKinActiv,vLEquActiv, & !IN
  & vXmk,vXeq,             & !IN
  & PhiF)                    !OUT
  !
  !~ IF(LogForAqu) THEN
  ALLOCATE(vSumNuAs(nAs))
  DO iAs=1,nAs
    X= Zero
    DO I=1,nCi
      IF(I/=isW) X= X+tNuAs(iAs,I)
    ENDDO
    vSumNuAs(iAs)= X
  ENDDO
  !~ ENDIF
  !
  !--------------------------------FINITE DIFFERENCE SYSTEM TO BE SOLVED
  !Psi(iCi)=
  !  N_iPi_(t) - N_iPi_(t-1) &
  !+ U.(dt/dX).Sum( Alfa(iCi,iAq).(n_iAq_(t) - n_iAq_(inj)) )
  !
  !Psi(iCi)=
  !  Sum(Alfa(iCi,iAq).n_iAq_t) + Sum(Alfa(iCi,iMn).n_iMn_t) - N_iPi_(t-1) &
  !+ (U.dt/dX).Sum(Alfa(iCi,iAq).n_iAq_t) - (U.dt/dX).Sum(Alfa(iCi,iAq).n_iAq_inj)
  !
  !Psi(iCi)=
  !  Sum(Alfa(iCi,iAq).n_iAq_t) . (1+(U.dt/dX)) &
  !+ Sum(Alfa(iCi,iMn).n_iMn_t) &
  !- N_iPi_(t-1) - U.dt/dX . n_iPi_inj
  !---------------------------------------------------------------------
  !
  vFunc=Zero
  !
  !-------------------------------Second'Aqu.Species (same as for Equil)
  IF(LogForAqu) THEN
    !----------------------------------------------------------LogForAqu
    DO iAs=1,nAs
      X=  vLnAct(1) *tNuAs(iAs,1) ! Solvent
      DO J=1,nCi
        IF(J/=1) X= X + vX(J) *tNuAs(iAs,J) ! Prim'Solute species
      ENDDO
      X= X - (vSumNuAs(iAs) -One)*(vX(1) +LOG(MWSv)) & ! mole nr -> molality
      &    - vDG_As(iAs)     & !
      &    - vDLGam_As(iAs)
      !
      !----------------"buffered" min'species -> add Sum(Nu_iMx*Act_iMx)
      DO J=1,nCx
        X=  X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      ENDDO
      !
      IF(DirectSub) THEN
        !--------------------------------------------direct substitution
        vMolSec(iAs)= EXP(X)
      ELSE
        ! Psi(nCi+1:nAq)=  equilibrium constraints (ChemicalAffinity=0)
        vFunc(nCi+iAs)= X - vX(nCi+iAs)
        vTolCoef(nCi+iAs)= MAX(One,ABS(vDG_As(iAs)+vDLGam_As(iAs)))
        !if(iDebug>2) write(73,'(I3,1X,G15.6)') nCi+iAs, vTolCoef(nCi+iAs)
      ENDIF
      !
    ENDDO
    !---------------------------------------------------------/LogForAqu
  ELSE
    !------------------------------------------------------NOT LogForAqu
    !~ print *,"size(vMolSec)=",size(vMolSec)  ;  pause
    DO iAs=1,nAs
      !
      X= vLnAct(1)*tNuAs(iAs,1) -vDG_As(iAs) -vDLGam_As(iAs)
      !
      !------------ "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
      DO J=1,nCx
        X= X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      ENDDO
      !
      X= EXP(X)
      !
      DO iPr=1,nCi
        IF(iPr /= 1) THEN
          IF(tNuAs(iAs,iPr)/=Zero) X= X *vX(iPr)**tNuAs(iAs,iPr)
          !& X= X *(vX(iPr)/vX(1)/MWSv)**tNuAs(iAs,iPr)
        ENDIF
      END DO
      vMolSec(iAs)= X *(vX(1)*MWSv)**(One -vSumNuAs(iAs))
      !~ print *,"iAs,vMolSec(iAs)=",iAs," _ ",vMolSec(iAs)  ;  pause
      !
      IF(.NOT. DirectSub) vFunc(nCi+iAs)= vMolSec(iAs) - vX(nCi+iAs)
      !
    END DO
    !-----------------------------------------------------/NOT LogForAqu
  ENDIF
  !--------------------------------------------------/Second'Aqu.Species
  !
  !-----------------------------------------------component conservation
  iPsi= 0
  !--------------------------------------------------------- iPsi- 1:nCi
  DO iPr=1,nCi
    !
    !---------------------------------------------mole nr injected /time
    IF (CoupledCoores) THEN
      QTotInjPr= vQTotInj(iPr) /Vbox
      ! mol.s-1 molar rate (explicit rates)
      ! Fout given
    ELSE
      ! nota: vTotInj= vTotInj /DotProd(vWeitCp,vTotInj) *RhoF *VBox
      ! -> mole/vol related to Box
      !
      ! UDarcy *Vbox /dX = fluid volume flushed through the box per time unit
      ! vTotInj(:) is in mole numbers for a volume VBox ->
      ! vQTotInj(:)= vTotInj(:) *UDarcy /dX
      !
      QTotInjPr= vTotInj(iPr) *UDarcy /dX   !! section S= (Vbox/dX)
      Fout=      UDarcy *Vbox /dX /PhiF     !! m3.s-1 flow related to porous volume
      !
    ENDIF
    !--------------------------------------------/mole nr injected /time
    !
    !! InFluid= DOT_PRODUCT(tAlfPr(iPr,1:nCi), vXf(1:nCi))          & ! # moles in fluid, from prim'species
    !! &      + DOT_PRODUCT(tAlfAs(iPr,1:nAs), vXf(nCi+1:nCi+nAs))    ! # moles in fluid
    !! InMin=   DOT_PRODUCT(tAlfKin(iPr,1:nMk),vVm(1:nMk)) * dTime    ! to/from minerals
    !
    !-------------------------------mole nrs in fluid, from prim'species
    TotFF= DOT_PRODUCT(tAlfPr(iPr,1:nCi),vXf(1:nCi))
    !
    !-------------------------------moles nrs in fluid, from sec'species
    DO iAs=1,nAs
      IF(DirectSub) THEN
        TotFF= TotFF + vMolSec(iAs) *tAlfAs(iPr,iAs)
      ELSE
        !TotFF= TotFF + EXP(vX(nCi+iAs)) *tAlfAs(iPr,iAs)
        IF(LogForAqu) THEN  ;  TotFF= TotFF + vXf(nCi+iAs) *tAlfAs(iPr,iAs)
        ELSE                ;  TotFF= TotFF + vX(nCi+iAs)  *tAlfAs(iPr,iAs)
        ENDIF
      ENDIF
    ENDDO
    !------------------------------/moles nrs in fluid, from sec'species
    !
    !--------------moles nrs in fluid, from the aqueous "mobile" species
    DO J=1,nAx
      TotFF= TotFF + MWSv *vXf(isW) *vMolal(J) *tAlfPr(iPr,nCi+J)
    ENDDO
    !---/
    !
    X=  TotFF *(One + dTime *Fout /Vbox) & ! Ci = (ni*Vbox )/ (Vbox*PhiW) = ni/Phiw
    & - vTotF(iPr)                       & ! # moles in fluid at previous step
    & - dTime *QTotInjPr                   ! # moles injection
    !
    !-------------------------------------------------------equil'phases
    IF(nEqA>0) THEN
      J=0
      DO iMk=1,nMk
        IF(vLEquActiv(iMk)) THEN
          J=J+1
          X= X + tAlfKin(iPr,iMk)*(vXeq(J) -vMolK(iMk))
        ENDIF
      ENDDO
    ENDIF
    !------------------------------------------------------/equil'phases
    !
    !---------------------------------------------------------kin'phases
    IF(nMkA>0) THEN
      DO iMk=1,nMk
        IF(vLKinActiv(iMk)) X= X + tAlfKin(iPr,iMk) *dTime *vVm(iMk)
      ENDDO
    ENDIF
    !--------------------------------------------------------/kin'phases
    !
    ! vFunc(iPr)= & !
    ! & ( One + dTime /dX *UDarcy /PhiF )           & ! Ci = (ni*Vbox )/ (Vbox*PhiW) = ni/Phiw
    ! * SUM(tAlfAq(iPr,1:nAq) * vXf(1:nAq)          & ! # moles in fluid
    ! - vTotF(iPr)                                  & ! # moles in fluid at previous step
    ! - dTime /dX *UDarcy *vTotInj(iPr)             & ! # moles injection
    ! + dTime *SUM(tAlfKin(iPr,1:nMk)*vVm(1:nMk))   & ! # moles from/to kin'phases
    ! + SUM(tAlfEqu(iPr,1:nEq)*(vVeq(:)-vMolEq(:))) & ! # moles from/to kin'phases
    !
    iPsi= iPsi+1
    !
    vTolCoef(iPsi)= MAX(One,ABS(TotFF))
    !if(iDebug>2) write(73,'(I3,1X,G15.6)') iPsi, vTolCoef(iPsi)
    !
    vFunc(iPsi)= X
    !
  ENDDO
  !----------------------------------------------/component conservation
  !
  if(iDebug>2 .AND. DirectSub) then
    write(71,'(A)') "==RESIDUAL============"
    do I=1,nAs
      write(71,'(A,G15.6)') "secc",vMolSec(I)
    enddo
    do I=1,nCi
      write(71,'(A,G15.6)') "prim",vXf(I)
    enddo
    do I=1,nCi
      write(71,'(A,G15.6)') "vFunc",vFunc(I)
    end do
  endif
  !
  IF(.NOT. DirectSub) iPsi= iPsi +nAs
  !
  !---------------------------------------------------equilibrium phases
  !---------------------------------------------------iPsi- nF+1:nF+nEqA
  IF(nEqA>0) THEN
    DO iMk=1,nMk
      IF(vLEquActiv(iMk)) THEN
        !
        !--- add one equation for each equ'phase
        X=   vLnAct(isW) *tNu_Kin(iMk,isW)                     & ! Solvent activity
        &  - SUM(tNu_Kin(iMk,2:nCi))*(vLnXf(isW) +LOG(MWSv))   & ! Solvent term (mole nr -> molality
        &  + DOT_PRODUCT(tNu_Kin(iMk,2:nCi),vLnGam(2:nCi))     & ! gamma prim'solute species
        &  - vDG_Kin(iMk)
        DO J=1,nCi
          IF(J/=isW) X= X + vLnXf(J) *tNu_Kin(iMk,J)
        ENDDO
        !---------- "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
        DO J=1,nCx
          X= X + vLnBuf(J) *tNu_Kin(iMk,nCi+J)
        ENDDO
        !---/
        !
        iPsi= iPsi+1
        vFunc(iPsi)= X
        !
        ! vFunc(iPsi)= &
        ! &             tNu_Kin(iMk,  isW) * vLnAct(isW)          & ! Solvent
        ! + DOT_PRODUCT(tNu_Kin(iMk,2:nCi), vX(2:nCi))            & ! prim' solute species
        ! + DOT_PRODUCT(tNu_Kin(iMk,2:nCi), vLnGam(2:nCi))        & ! gamma prim'solute species
        ! -         SUM(tNu_Kin(iMk,2:nCi))*(vX(isW) +LOG(MWSv))  & ! Solvent term (mole nr -> molality
        ! - vDG_Kin(iMk)
        ! !
        ! !---------- "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
        ! IF(nCx>0) &
        ! & vFunc(iPsi)= vFunc(iPsi) &
        ! + DOT_PRODUCT(tNu_Kin(iMk,nCi+1:nCi+nCx),vLnBuf(1:nCx))
      ENDIF
    ENDDO
  ENDIF
  !--------------------------------------------------/equilibrium phases
  !
  !-------------------------------------------------------kinetic phases
  !---------------------------- iPsi- nCi+nAs+nEqA+1 : nCi+nAs+nEqA+nMkA
  IF(nMkA>0) THEN
    iMkA= 0
    DO iMk=1,nMk
      IF(vLKinActiv(iMk)) THEN
        iMkA= iMkA + 1
        !--- add one equation for each kin'phase
        iPsi= iPsi+1
        vFunc(iPsi)= vXmk(iMkA) -vMolK(iMk) -dTime *vVm(iMk)
        vTolCoef(iPsi)= MAX(One,ABS(vMolK(iMk)))
        if(iDebug>2) write(73,'(I3,1X,G15.6)') iPsi, vTolCoef(iPsi)
        !---/
      ENDIF
    ENDDO
  ENDIF
  !------------------------------------------------------/kinetic phases
  !
  IF(nAx>0) DEALLOCATE(vMolal)
  !
  DEALLOCATE(vXf)
  DEALLOCATE(vLnXf)
  DEALLOCATE(vVm)
  DEALLOCATE(vXmk)
  DEALLOCATE(vXeq)
  DEALLOCATE(vSumNuAs)
  !
  !do i=1,size(vFunc)
  !  write(81,'(G12.3,1X)',advance="no") vFunc(I)
  !enddo
  !write(81,*)
  !Dynam_Residual= vFunc
  !
  RETURN
ENDFUNCTION Dynam_Residual

SUBROUTINE KinRate_Calc( & !
& dTime,                 & !IN
& vLKinActiv,            & !IN
& vX,vXmk,               & !IN
!
& vSurfK,vVmQsK,vVmAct,  & !INOUT
& dSRdXm,dSRdLnXm,       & !OUT
& dVmQdLnXf,             & !OUT
& dVMdLnXf,dVMdLnXm)       !OUT
!--
!-- compute mineral rates --
!--
  !
  USE M_T_KinFas,   ONLY: KinFas_CalcSurf
  USE M_KinFas_Surf,ONLY: KinFas_Surf_Calc
  USE M_KinRate
  !
  USE M_Global_Vars,ONLY: vKinFas,nAq,vFas
  USE M_Basis_Vars, ONLY: MWSv,isW,nCi,nAx,nMx
  USE M_Basis_Vars, ONLY: tAlfPr,tAlfAs,tNuAs
  !
  USE M_Dynam_Vars, ONLY:  & !
  & Implicit_Surface,      & !
  & Implicit_ActivFactor,  & !
  & QsK_Iota, VBox,        & !
  & vLnGam,vLnBuf,vLnAct,  & !
  & vKinMod,vMolarVol,vKinMinim, & !
  & vMolK0,vSurfK0,        & !
  & tNu_Kin,tAlfKin,       & !
  & vDG_Kin,vDG_As
  !
  REAL(dp),INTENT(IN) :: dTime
  LOGICAL, INTENT(IN) :: vLKinActiv(:)
  REAL(dp),INTENT(IN) :: vX(:),vXmk(:)
  REAL(dp),INTENT(INOUT):: vSurfK(:),vVmQsK(:),vVmAct(:)
  REAL(dp),INTENT(OUT):: dSRdXm(:,:),dSRdLnXm(:,:)
  REAL(dp),INTENT(OUT):: dVmQdLnXf(:,:)
  REAL(dp),INTENT(OUT):: dVMdLnXf(:,:),dVMdLnXm(:,:)
  !
  REAL(dp),DIMENSION(1:nCi):: dQsKdLnXi
  REAL(dp):: QsK,QskIota,varMaxM,X
  INTEGER :: nMkA,nMk,nCx
  INTEGER :: iPr,iMk,iMkA
  CHARACTER(LEN=7):: cSatur
  !
  nCx=  nAx+nMx
  nMk= SIZE(vKinFas)
  nMkA= COUNT(vLKinActiv)
  !
  vVmQsK=   Zero
  dSRdLnXm= Zero
  dSRdXm=   Zero
  dVmQdLnXf=Zero
  dVMdLnXf(:,1:nCi)= Zero
  dVMdLnXm(:,1:nMk)= Zero
  !
  !IF(Implicit_ActivFactor) THEN !update activities
  !  !Nota: as Activity Coeffs are not Implicited, Rate_Act is better not Implicited
  !  dVmAdLnXf=Zero
  !  vLnAct(2:nAq)=vX(2:nAq)+vLnGam(2:nAq)-vX(isW)-LOG(MWSv)
  !ENDIF
  !
  QskIota= QsK_Iota !Zero !mod 10/06/2008 17:02
  !
  ! IF(nMkA + nEqA >0) THEN
  IF(nMkA > 0) THEN
    !
    varMaxM= Zero
    iMkA= 0
    iMk=  0
    !
    DO
      !
      iMk= iMk+1
      IF(iMk>SIZE(vKinFas)) EXIT
      !
      !--------------------------------------------iMk is kinetic active
      IF(vLKinActiv(iMk)) THEN
        !
        iMkA= iMkA + 1
        !
        !--------------------------------------------------calc. surface
        IF(Implicit_Surface) THEN
          !X= vXmk(iMkA)
          X= MAX(vXmk(iMkA),vKinMinim(iMkA))
          CALL KinFas_Surf_Calc( & !
          & .TRUE.,             & !IN: bImplicit
          & vKinFas(iMk)%cMode, & !IN: kinetic model
          & X,                  & !IN: Nr Moles of mineral
          & vMolK0(iMk),        & !IN: Nr Moles of phase at given state
          & vSurfK0(iMk),       & !IN: Surface of phase at that state
          & vSurfK(iMk),        & !OUT: Surface
          & dSRdLnXm(iMk,iMk),  & !OUT:
          & dSRdXm(iMk,iMk))      !OUT:
        ENDIF
        !-------------------------------------------------/calc. surface
        !
        !----------------------------------calc mineral rate per surface
        !--------------------------------------(- only "chemical" terms)
        CALL KinRate_CalcQsK( & !
        & nCi,nCx,        & !
        & vDG_Kin(iMk),   & !
        & tNu_Kin(iMk,:), & !
        & vX,             & !IN:  Ln(aqu.species mole nrs) -> used for CalcQsK
        & vLnGam,         & !IN:  Ln(ActCoeff,aqu.species)
        & vLnBuf,         & !IN:  Ln(ActivitySpecies) of buffered species
        & vLnAct(isW),    & !IN:  Ln(ActivitySolvent)
        & QsK,            & !OUT, mineral saturation
        & dQsKdLnXi)        !OUT
        !---------------------------------/calc mineral rate per surface
        !
        !-----------------------------------------from QsK deduce cSatur
        !---------------------------------for branching to dissol/precip
        CALL KinRate_SatState( &
        & vKinFas(iMk),   & !IN:  mineral: for cSatur, QsKseuil
        & vXmk(iMkA),     & !IN:  Nr Moles of mineral, used to check nMol<=MolMinim
        !& vKinMinim(iMk), & !IN:
        & Zero,           & !IN: -> cSatur will not be "MINIMAL"
        & QsK,            & !IN:
        & QskIota,        & !IN
        & cSatur)           !OUT
        !----------------------------------------/from QsK deduce cSatur
        !
        !----------------------------------------compute affinity factor
        !! IF(cSatur/="MINIMAL") & !
        CALL KinRate_CalcQsKFactor( &
        & cSatur,      & !IN
        & vKinMod(vKinFas(iMk)%iKin),&  !IN
        & QsK,         & !IN
        & dQsKdLnXi,   & !IN
        & vVmQsK(iMk), & !OUT: the QsK factor in the rate law
        & dVmQdLnXf(iMk,1:nCi))
        !---------------------------------------/compute affinity factor
        !
        !---------------------------------------------------calc. vVmAct
        !provisionally, dVmAdLnX_M is not computed in this new version  !!!!!!!!!!!!!
        !-> we assume that this factor will not be implicited           !!!!!!!!!!!!!
        !!!  IF(Implicit_ActivFactor) & !else will use current values obtained in Dynam_OneStep_KinFasUpd
        !!!  & CALL KinRate_CalcActivFactor(&
        !!!  & vKinFas(iMk)%cSatur, & !IN
        !!!  & vKinMod(vKinFas(iMk)%iKin), & !IN: kinetic PARAMETERs
        !!!  & vLnAct,           & !LnActH_,LnActOH,LnActCO2, & !activators
        !!!  & vVmAct(iMk),      & !OUT
        !!!  & dVmAdLnXf(iMk,:))   !OUT
        !
        !! IF(TestMax .AND. vMolK(iMk)>1.0D-6) &
        !! & varMaxM=MAX(varMaxM, ABS( dTime*vVm(iMk) / vMolK(iMk) ))
        !
      ENDIF
      !-------------------------------------------/iMk is kinetic active
      !
    ENDDO
    !
    !! IF(TestMax .AND. varMaxM>1.0D0) THEN
    !!   dTime=dTime * 0.5_dp
    !!   IF(iDebug>0) WRITE(fTrc,'(A)') "dtime adjusted in vFunc"
    !! ENDIF
  ENDIF
  !
  !----------------------------------------------------/CalcMineralRates
END SUBROUTINE KinRate_Calc

SUBROUTINE PhiF_Calc( & !
!--
!-- compute porosity
!--
& PhiInert,              & !IN
& vLKinActiv,vLEquActiv, & !IN
& vXmk,vXeq,             & !IN
& PhiF)                    !OUT
  !
  USE M_T_KinFas,   ONLY: KinFas_CalcSurf
  USE M_KinFas_Surf,ONLY: KinFas_Surf_Calc
  USE M_KinRate
  !
  USE M_Global_Vars,ONLY: vKinFas
  USE M_Dynam_Vars, ONLY: VBox, vMolarVol
  !
  REAL(dp),INTENT(IN) :: PhiInert
  LOGICAL, INTENT(IN) :: vLKinActiv(:),vLEquActiv(:)
  REAL(dp),INTENT(IN) :: vXmk(:),vXeq(:)
  REAL(dp),INTENT(OUT):: PhiF
  !
  INTEGER :: nMk
  INTEGER :: iMk,iMkA,iEqA
  !
  nMk= SIZE(vKinFas)
  !
  PhiF= One - PhiInert
  !-> take account of minerals that are "inactive" during current time step
  !
  IF(COUNT(vLKinActiv) + COUNT(vLEquActiv) >0) THEN
    !
    iMkA= 0
    iEqA= 0
    iMk=  0
    !
    DO
      iMk= iMk+1
      IF(iMk>nMk) EXIT
      !
      !---------------------------------------- iMk is kinetic active --
      IF(vLKinActiv(iMk)) THEN
        iMkA= iMkA + 1
        PhiF= PhiF - vXmk(iMkA) * vMolarVol(iMk) /VBox
      ENDIF
      !---------------------------------------/ iMk is kinetic active --
      !
      !------------------------------------ iMk is equilibrium active --
      IF(vLEquActiv(iMk)) THEN
        iEqA= iEqA +1
        PhiF= PhiF - vXeq(iEqA) * vMolarVol(iMk) /VBox
      ENDIF
      !-----------------------------------/ iMk is equilibrium active --
      !
    ENDDO
    !
  ENDIF
  !
END SUBROUTINE PhiF_Calc

ENDMODULE M_Dynam_Residual
