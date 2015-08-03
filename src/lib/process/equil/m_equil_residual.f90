MODULE M_Equil_Residual
!--
!-- implements the equation system
!-- for speciation and equilibrium calculations
!--
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  USE M_Kinds
  USE M_Numeric_Const, ONLY: MinExpDP,MaxExpDP
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Residual
  PUBLIC:: Equil_Converge
  PUBLIC:: vTolCoef
  !
  REAL(dp),ALLOCATABLE:: vTolCoef(:)
  !
CONTAINS

LOGICAL FUNCTION Equil_Converge(vF,vTolF)
!--
!-- makes possible adapatation of convergence criteria to type of condition
!-- e.g. material conservation vs potential, etc.
!--
  REAL(dp),INTENT(IN):: vF(:),vTolF(:)
  !
  Equil_Converge= ALL(ABS(vF(:)*vTolCoef(:))<vTolF(:))
  !
ENDFUNCTION Equil_Converge

FUNCTION Equil_Residual(vX)
!--
!-- equation system to be solved for speciation and equilibrium calculations
!-- find vX(1:nF+nM) that makes Equil_Residual(1:nF+nM)-0(1:nF+nM)
!--
  USE M_Safe_Functions
  USE M_IOTools
  !
  USE M_Global_Vars,ONLY: vSpc
  USE M_Basis_Vars,ONLY: &
  & iOx,isW,MWSv, &
  & tAlfPr,tAlfAs,tNuAs,tAlfFs,tNuFas, &
  & vOrdPr,vOrdAq,vOrdAs, &
  & nCi,nAx,nAs,nMx
  !
  !~ USE M_Equil_Vars,ONLY: OsmoSv,UseOsmoSv
  USE M_Equil_Vars,ONLY: vMolal,vMolSec,vFasMole,vLnAct,vLnGam,vTotS
  USE M_Equil_Vars,ONLY: vFasAff,vFasRat,vAffScale
  USE M_Equil_Vars,ONLY: vDeltaG_Fs,vDGapp_As
  USE M_Equil_Vars,ONLY: cEquMode,nEvalFunc,dXi
  USE M_Equil_Vars,ONLY: LogForAqu,DirectSub,Complementarity
  USE M_Equil_Vars,ONLY: vDeltaG_Eq,tAlfEq,tNuEq,nEquFas
  !
  REAL(dp),DIMENSION(:),INTENT(IN) :: vX
  REAL(dp),DIMENSION(SIZE(vX))     :: Equil_Residual
  !
  REAL(dp),DIMENSION(SIZE(vX)):: vFunc
  !
  REAL(dp),ALLOCATABLE:: vLnXf(:),vXf(:),vXm(:)
  REAL(dp),ALLOCATABLE:: vSumNuAs(:)
  !
  INTEGER:: nF !nr of aqu'species in array vX of residual
  INTEGER:: nM !nr of non-aqu'species in array vX of residual
  !
  INTEGER :: iAs,iAx,iCi,I,J
  REAL(dp):: LnActW,X !, xAff !ln(Solmodel_activity)
  REAL(dp):: AffScale

  IF(iDebug>2) nEvalFunc= nEvalFunc +1

  vFunc=Zero

  IF(DirectSub) THEN  ;   nF= nCi
  ELSE                ;   nF= nCi +nAs
  ENDIF

  ALLOCATE(vXf(nF))  ;  vXf(:)= Zero

  nM= nEquFas

  IF(nM>0) THEN
    ALLOCATE(vXm(nM))
    vXm(1:nM)= vX(nF+1:nF+nM)
  ENDIF
  !
  ALLOCATE(vLnXf(nF))
  IF(LogForAqu) THEN
    vLnXf(1:nF)= vX(1:nF)
    vXf(:)= FSafe_vExp(vLnXf(:)) ! avoid overflow with EXP
  ELSE
    vXf(1:nF)= vX(1:nF)
    vLnXf(1:nF)= FSafe_vLog(vX(1:nF))
  ENDIF
  !IF(iOx/=0) THEN
  !  vTotS(iOx)= SUM(tStoikio(iOx,1:nAq)*vSpc(1:nAq)%Mole)
  !ENDIF
  !
  ! alkalinity
  ! Alc = HCO3- + 2*CO3-2 + OH- - H+
  ! -> condition on vMolF(H+)
  !
  !-------------------- update molalities of aqueous "mobile" species --
  !----------- used to calculate mole number- molality * solvent'mass --
  DO I=1,nAx
    vMolal(vOrdPr(nCi+I))= EXP(vLnAct(vOrdPr(nCi+I))-vLnGam(vOrdPr(nCi+I)))
  ENDDO
  !-------------------/ update molalities of aqueous "mobile" species --
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
  ! caveat!
  ! Solvent's gamma refers to
  ! deviation from standard state (= pure solvent) along a mole fraction scale
  !
  ! -> for Solvent,
  !   vLnAct(isW) = vLnGam(isW) &
  !               + vX(isW)/( SUM(mole numbers all species, including Solvent) )
  !
  ! -> using vLnGam(isW) would need SUM(EXP(vX(1:nF)))+SUM(mobileSp)
  ! -> for this reason, considering vLnAct(isW) does not change much within Newton,
  !    it is taken constant from its last evaluation by Solmodel_CalcGamma in outer loop
  !
  !!??? preferable to work with fixed OsmoSv, instead of fixed vLnAct(isW) ???
  !!vLnAct(isW)= - OsmoSv * SUM(#mole_solute_species) / #mole_solvent
  !!vLnAct(isW)= - OsmoSv * ( SUM(vXf(2:nAs)) /vX(isW) + MWSv*SUM(vMolal(vOrdPr(nCi+1:nCi+nAx))) )
  !!                          "internal" solute'sp     + "external" solute'sp
  !
  !~ IF(UseOsmoSv) THEN
    !~ LnActW= - OsmoSv &
    !~ &       * ( SUM(vXf(2:nAs)) /vXf(1) + MWSv*SUM(vMolal(vOrdPr(nCi+1:nCi+nAx))) )
  !~ ELSE
    !~ LnActW= vLnAct(isW)
  !~ ENDIF
  !
  LnActW= vLnAct(isW)
  !
  !---------------------- chemical equilibrium for Second'Aqu.Species --
  IF(LogForAqu) THEN
    !--------------------------- sec'species, USING LOG(MOLE NUMBERS) --
    DO iAs=1,nAs
      !
      ! vDGapp_As(iAs)=
      ! Nu(iAs,isW)*LnActW + sum( Nu(iAs,iPr) *[vX(iPr)-vX(isW)-ln(MWsv)] )
      ! -[ vX(iAs) -vX(isW) -ln(MWsv) ]
      !
      X  = LnActW *tNuAs(iAs,isW) &
      &  + (One - vSumNuAs(iAs))*(vX(isW) +LOG(MWSv)) &
      &  - vDGapp_As(iAs)
      !
      !-- "buffered" aqu'species -> add Sum(Nu_iAx*ln(Act_iAx))
      !-- "buffered" min'species -> add Sum(Nu_iMx*ln(Act_iMx))
      !here in "mass action", we take account of all mobile species,
      !                       including non-aqueous, nF+1:nF+nAx+nMx,
      !in contrast to "mass balance", which takes ONLY aqu'species: nF+1:nF+nAx
      DO iAx=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+iAx)) *tNuAs(iAs,nCi+iAx)
      ENDDO
      !
      DO iCi=1,nCi
        IF(iCi/=isW) X= X +tNuAs(iAs,iCi)*vX(iCi)
      ENDDO
      !
      IF(DirectSub) THEN
        ! X is log(mole nr sec'species)
        !~ vMolal(iAs)= EXP(X) /vX(isW) /MWSv  != molality nrs sec'species
        vMolSec(iAs)= EXP(X) != molality nrs sec'species
      ELSE
        !---------------------- Psi(nCi+iAs)-  equilibrium constraint --
        vFunc(nCi+iAs)= X  - vX(nCi+iAs)
      ENDIF
      !
    ENDDO
    !---------------------------/sec'species, USING LOG(MOLE NUMBERS) --
  ELSE
    !-------------------------------- sec'species, USING MOLE NUMBERS --
    DO iAs=1,nAs
      !
      X = LnActW *tNuAs(iAs,isW) &
      & - vDGapp_As(iAs)
      DO iAx=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+iAx)) *tNuAs(iAs,nCi+iAx)
      ENDDO
      !
      X= EXP(X)
      !
      DO iCi=1,nCi
        IF(iCi/=isW .AND. tNuAs(iAs,iCi)/= 0.0D0) X= X *vX(iCi)**tNuAs(iAs,iCi)
        ! & X= X *(vX(iCi)/vX(isW)/MWSv)**tNuAs(iAs,iCi)
      END DO
      !~ X= X /(vX(isW)*MWSv)**vSumNuAs(iAs)
      !~ vMolal(iAs)= X
      vMolSec(iAs)= X *(vX(isW)*MWSv)**(One -vSumNuAs(iAs))
      !
      !------------------------ Psi(nCi+iAs)-  equilibrium constraint --
      !~ IF(.NOT. DirectSub) vFunc(nCi+iAs)= X *(vX(isW) *MWSv) - vX(nCi+iAs)
      IF(.NOT. DirectSub) vFunc(nCi+iAs)= vMolSec(iAs) - vX(nCi+iAs)
      !
    END DO
    !--------------------------------/sec'species, USING MOLE NUMBERS --
  ENDIF
  !----------------------/chemical equilibrium for Second'Aqu.Species --
  !
  !-------------------------------------------- material conservation --
  IF(LogForAqu) THEN
    !------------------------------- material conservation, LogForAqu --
    DO iCi=1,nCi
      !
      !-------------------------------- contrib' from primary species --
      X= SUM(tAlfPr(iCi,1:nCi)*vXf(1:nCi))
      !------------------------------ contrib' from secondary species --
      DO iAs=1,nAs
        IF(tAlfAs(iCi,iAs)/= 0.0D0) THEN
          IF(DirectSub) THEN  ;  X= X + tAlfAs(iCi,iAs)*vMolSec(iAs)
          ELSE                ;  X= X + tAlfAs(iCi,iAs)*vXf(nCi+iAs)
          ENDIF
        ENDIF
      ENDDO
      !------------------- contrib' from the aqueous "mobile" species --
      DO I=1,nAx
        X= X + tAlfPr(iCi,nCi+I) *vMolal(vOrdPr(nCi+I)) *vXf(isW) *MWSv
      ENDDO
      !---/
      !
      vFunc(iCi)= X - vTotS(iCi)
      !
    ENDDO
    !
    !~ if(iDebug>2) then
      !~ write(71,'(A)') "==RESIDUAL============"
      !~ do I=1,nAs
        !~ write(71,'(A,G15.6)') "secc",vMolSec(I)
      !~ enddo
      !~ do I=1,nCi
        !~ write(71,'(A,G15.6)') "prim",vXf(I)
      !~ enddo
      !~ do I=1,nCi
        !~ write(71,'(A,G15.6)') "vFunc",vFunc(I)
      !~ end do
    !~ endif
    !
    !~ IF(bMolal) vFunc(isW)= vXf(isW) - One/MWSv
    !------------------------------/ material conservation, LogForAqu --
  !
  ELSE
    !--------------------------- material conservation, NOT LogForAqu --
    DO iCi=1,nCi
      !
      !-------------------------------- contrib' from primary species --
      X= SUM(tAlfPr(iCi,1:nCi)*vX(1:nCi))
      !
      !------------------------------ contrib' from secondary species --
      DO iAs=1,nAs
        IF(tAlfAs(iCi,iAs)/= 0.0D0) THEN
          IF(DirectSub) THEN  ;  X= X + tAlfAs(iCi,iAs) *vMolSec(iAs)
          ELSE                ;  X= X + tAlfAs(iCi,iAs) *vX(nCi+iAs)
          ENDIF
        ENDIF
      END DO
      !------------------- contrib' from the aqueous "mobile" species --
      DO iAx=1,nAx
        X= X + tAlfPr(iCi,nCi+iAx)*vMolal(vOrdPr(nCi+iAx)) *vX(isW)*MWsv
      ENDDO
      !
      vFunc(iCi)= X -vTotS(iCi)
      ! in SPC, vTotS = vTot_inFluid
      ! in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
      !
    END DO
    !
    if( DirectSub .AND. iDebug>2) then
      write(71,'(A)') "==RESIDUAL============"
      do iAs=1,nAs
        !~ write(71,'(A,G15.6)') "secc",vMolal(iAs)*(vX(isW)*MWSv)
        write(71,'(A,G15.6)') "secc",vMolSec(iAs)
      enddo
      do iCi=1,nCi
        write(71,'(A,G15.6)') "prim",vX(iCi)
      enddo
      do iCi=1,nCi
        write(71,'(A,G15.6)') "vFunc",vFunc(iCi)
      end do
    endif
    !--------------------------/ material conservation, NOT LogForAqu --
  ENDIF
  !-------------------------------------------/ material conservation --
  !
  !----------------------in case of "global equilibrium" calculation, --
  !-------------------------------- contribution from non aqu'species --
  IF(nM>0) THEN
    SELECT CASE(cEquMode)
    !
    ! cEquMode is flag for the type of "equilibrium speciation", i.e.
    ! EQ4= reaction progress method, titration-like
    ! EQ2= direct method
    !
    !--------------------------------------- reaction progress method --
    CASE("EQ1") !
      DO I=1,nEquFas !nFs
        AffScale= SUM(ABS(tNuEq(I,:)))
        !----------------------------------- compute phase affinities --
        X=   vDeltaG_Eq(I)    & !
        &  - tNuEq(I,1) *vLnAct(isW)  &                    ! solvent
        &  + SUM(tNuEq(I,2:nCi)) *(vLnXf(isW) +LOG(MWSv))  ! mole nr to molality
        DO J=2,nCi
          X= X - ( vLnXf(J) +vLnGam(vOrdAq(J)) )*tNuEq(I,J)
        ENDDO
        !----------------------------- add terms for mobile species --
        DO J=1,nAx+nMx
          X= X -tNuEq(I,nCi+J) *vLnAct(vOrdPr(nCi+J))
        ENDDO
        !
        vFasAff(I)= X
        !--------------------------------/ compute phase affinities --
        !
        !------------ mineral production rate (precip' - vFasRat>0) --
        !----------------------- rate - EXP(Affinity) - 1 - Q/K - 1 --
        vFasRat(I)= FSafe_Exp(-vFasAff(I) /AffScale) - One
        !
        !------- add mineral production in mass balance constraints --
        DO iCi=1,nCi !solvent involved in mass balance
        !DO iCi=2,nCi  !solvent not involved in mass balance
          vFunc(iCi)= vFunc(iCi) &
          &         + tAlfEq(iCi,I) &
          &           *(vFasMole(I) + vFasRat(I) *dXI) !XI= progress param'
        ENDDO
        !
        !--- add material balance equations --
        !--- Rate(xi) - [vFasMole(xi+dxi) - vFasMole(xi)] / dxi
        vFunc(nF+I)= vX(nF+I) - vFasMole(I) - dXi *vFasRat(I)
        !
      ENDDO ! DO iFs
      !---------------------------------------/ loop on active phases --
    !--------------------------------------/ reaction progress method --
    !
    CASE("EQ2") !
      !
      IF(nEquFas>0) THEN
        DO I=1,nEquFas
          !--- POTENTIAL BALANCE--
          !-- add potential balance equations imposed by non-aqu'phase
          !-- i.e. impose vFasAff(I)-Zero
          X=  vDeltaG_Eq(I) &
          -             tNuEq(I,isW)    *vLnAct(isW)            & ! solvent
          - DOT_PRODUCT(tNuEq(I,2:nCi),  vLnXf(2:nCi)        )  & !- vX(nCi+I)
          +         SUM(tNuEq(I,2:nCi))*(vLnXf(isW) +LOG(MWSv)) & !mole nr -> molality
          - DOT_PRODUCT(tNuEq(I,2:nCi),  vLnGam(vOrdAq(2:nCi)))
          !
          DO J=1,nAx ! aqu'mobile prim'species
            X= X - tNuEq(I,nCi+J) *vLnAct(vOrdPr(nCi+J))
          ENDDO
          DO J=1,nMx ! min'mobile prim'species
            X= X - tNuEq(I,nCi+nAx+J) *vLnAct(vOrdPr(nCi+nAx+J))
          ENDDO
          !
          IF(Complementarity) THEN
            ! vFunc(nF+I)= MIN(X,vX(nF+I))
            !~ IF(X < vX(nF+I)) THEN  ;  vFunc(nF+I)= X
            !~ ELSE                   ;  vFunc(nF+I)= vX(nF+I)
            !~ ENDIF
            ! vFunc(nF+I)= Fischer - Burmeister
            vFunc(nF+I)= X + vX(nF+I) - SQRT(X**2 + vX(nF+I)**2)
          ELSE
            vFunc(nF+I)= X
          ENDIF
          !---/
          !
          !--- MATERIAL BALANCE--
          !-- add contribution of non aqu'phase I on material balance
          !-- equations
          DO iCi=1,nCi
            vFunc(iCi)= vFunc(iCi) + tAlfEq(iCi,I) *vX(nF+I)
          ENDDO
          !---/
        ENDDO
      ENDIF
    !-------------------------------------------------/ direct method --
    !
    ENDSELECT
    !
  ENDIF ! IF(nM>0)
  !-------------------------------/ contribution from non aqu'species --
  !
  DEALLOCATE(vXf)
  DEALLOCATE(vLnXf)
  DEALLOCATE(vSumNuAs)
  IF(nM>0) DEALLOCATE(vXm)
  !
  Equil_Residual= vFunc
  !
  RETURN
ENDFUNCTION Equil_Residual

ENDMODULE M_Equil_Residual


