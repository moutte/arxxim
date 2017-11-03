module M_Equil_Residual
!--
!-- implements the equation system
!-- for speciation and equilibrium calculations
!--
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  use M_Kinds
  use M_Numeric_Const, only: MinExpDP,MaxExpDP
  implicit none
  !
  private
  !
  public:: Equil_Residual
  public:: Equil_Converge
  public:: vTolCoef
  !
  real(dp),allocatable:: vTolCoef(:)
  !
contains

logical function Equil_Converge(vF,vTolF)
!--
!-- makes possible adapatation of convergence criteria to type of condition
!-- e.g. material conservation vs potential, etc.
!--
  real(dp),intent(in):: vF(:),vTolF(:)
  !
  Equil_Converge= all(abs(vF(:)*vTolCoef(:))<vTolF(:))
  !
end function Equil_Converge

function Equil_Residual(vX)
!--
!-- equation system to be solved for speciation and equilibrium calculations
!-- find vX(1:nF+nM) that makes Equil_Residual(1:nF+nM)-0(1:nF+nM)
!--
  use M_Safe_Functions
  use M_IOTools
  !
  use M_Global_Vars,only: vSpc,SolModel
  !
  use M_Basis_Vars,only: tAlfPr,tAlfAs,tNuAs,tAlfFs,tNuFas
  use M_Basis_Vars,only: vOrdPr,vOrdAq,vOrdAs
  use M_Basis_Vars,only: nCi,nAx,nAs,nMx
  !
  !! use M_Equil_Vars,only: OsmoSv,UseOsmoSv
  use M_Equil_Vars,only: vMolal,vMolSec,vFasMole,vLnAct,vLnGam,vTotS
  use M_Equil_Vars,only: vFasAff,vFasRat,vAffScale
  use M_Equil_Vars,only: vDeltaG_Fs,vDGapp_As
  use M_Equil_Vars,only: cEquMode,nEvalFunc,dXi
  use M_Equil_Vars,only: LogForAqu,DirectSub,Complementarity
  use M_Equil_Vars,only: vDeltaG_Eq,tAlfEq,tNuEq,nEquFas
  !---------------------------------------------------------------------
  real(dp),dimension(:),intent(in) :: vX
  real(dp),dimension(size(vX))     :: Equil_Residual
  !---------------------------------------------------------------------
  real(dp),dimension(size(vX)):: vFunc
  !
  real(dp),allocatable:: vLnXf(:),vXf(:),vXm(:)
  real(dp),allocatable:: vSumNuAs(:)
  !
  integer:: nF !nr of aqu'species in array vX of residual
  integer:: nM !nr of non-aqu'species in array vX of residual
  !
  integer :: iAs,iAx,iCi,I,J,isW
  real(dp):: LnActW,X,MWSv !, xAff !ln(Solmodel_activity)
  real(dp):: AffScale
  !---------------------------------------------------------------------
  if(iDebug>2) nEvalFunc= nEvalFunc +1

  isW= SolModel%iSolvent
  MWSv = SolModel%MolWeitSv
  
  vFunc=Zero

  if(DirectSub) then  ;   nF= nCi
  else                ;   nF= nCi +nAs
  end if

  allocate(vXf(nF))  ;  vXf(:)= Zero

  nM= nEquFas

  if(nM>0) then
    allocate(vXm(nM))
    vXm(1:nM)= vX(nF+1:nF+nM)
  end if
  !
  allocate(vLnXf(nF))
  if(LogForAqu) then
    vLnXf(1:nF)= vX(1:nF)
    vXf(:)= FSafe_vExp(vLnXf(:)) ! avoid overflow with EXP
  else
    vXf(1:nF)= vX(1:nF)
    vLnXf(1:nF)= FSafe_vLog(vX(1:nF))
  end if
  !
  ! alkalinity
  ! Alc = HCO3- + 2*CO3-2 + OH- - H+
  ! -> condition on vMolF(H+)
  !
  !----------------------update molalities of aqueous "mobile" species--
  !-------------used to calculate mole number- molality * solvent'mass--
  do I=1,nAx
    vMolal(vOrdPr(nCi+I))= exp(vLnAct(vOrdPr(nCi+I))-vLnGam(vOrdPr(nCi+I)))
  end do
  !--/
  !
  allocate(vSumNuAs(nAs))
  do iAs=1,nAs
    X= Zero
    do I=1,nCi
      if(I/=isW) X= X+tNuAs(iAs,I)
    end do
    vSumNuAs(iAs)= X
  end do
  !
  ! caveat!
  ! Solvent's gamma refers to
  ! deviation from standard state (= pure solvent) along a mole fraction scale
  !
  ! -> for Solvent,
  !   vLnAct(isW) = vLnGam(isW) &
  !               + vX(isW)/( sum(mole numbers all species, including Solvent) )
  !
  ! -> using vLnGam(isW) would need sum(exp(vX(1:nF)))+sum(mobileSp)
  ! -> for this reason, considering vLnAct(isW) does not change much within Newton,
  !    it is taken constant from its last evaluation by Solmodel_CalcGamma in outer loop
  !
  !!??? preferable to work with fixed OsmoSv, instead of fixed vLnAct(isW) ???
  !!vLnAct(isW)= - OsmoSv * sum(#mole_solute_species) / #mole_solvent
  !!vLnAct(isW)= - OsmoSv * ( sum(vXf(2:nAs)) /vX(isW) + MWSv*sum(vMolal(vOrdPr(nCi+1:nCi+nAx))) )
  !!                          "internal" solute'sp     + "external" solute'sp
  !
  !! if(UseOsmoSv) then
  !!   LnActW= - OsmoSv &
  !!   &       * ( sum(vXf(2:nAs)) /vXf(1) + MWSv*sum(vMolal(vOrdPr(nCi+1:nCi+nAx))) )
  !! else
  !!   LnActW= vLnAct(isW)
  !! end if
  !
  LnActW= vLnAct(isW)
  ! print *,"vLnAct(isW)=",  LnActW   ;  call pause_             !!DEBUG!!
  !
  !--------------------------chemical equilibrium for Second'Aqu.Species
  if(LogForAqu) then
    !-------------------------------sec'species, USING log(MOLE NUMBERS)
    do iAs=1,nAs
      !
      ! vDGapp_As(iAs)=
      ! Nu(iAs,isW)*LnActW + sum( Nu(iAs,iPr) *[vX(iPr)-vX(isW)-ln(MWsv)] )
      ! -[ vX(iAs) -vX(isW) -ln(MWsv) ]
      !
      X  = LnActW *tNuAs(iAs,isW) &
      &  + (One - vSumNuAs(iAs))*(vX(isW) +log(MWSv)) &
      &  - vDGapp_As(iAs)
      !
      !-- "buffered" aqu'species -> add Sum(Nu_iAx*ln(Act_iAx))
      !-- "buffered" min'species -> add Sum(Nu_iMx*ln(Act_iMx))
      !here in "mass action", we take account of all mobile species,
      !                       including non-aqueous, nF+1:nF+nAx+nMx,
      !in contrast to "mass balance", which takes only aqu'species: nF+1:nF+nAx
      do iAx=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+iAx)) *tNuAs(iAs,nCi+iAx)
      end do
      !
      do iCi=1,nCi
        if(iCi/=isW) X= X +tNuAs(iAs,iCi)*vX(iCi)
      end do
      !
      if(DirectSub) then
        ! X is log(mole nr sec'species)
        !! vMolal(iAs)= exp(X) /vX(isW) /MWSv  != molality nrs sec'species
        vMolSec(iAs)= exp(X) != molality nrs sec'species
      else
        !---------------------- Psi(nCi+iAs)-  equilibrium constraint --
        vFunc(nCi+iAs)= X  - vX(nCi+iAs)
      end if
      !
    end do
    !------------------------------/sec'species, USING log(MOLE NUMBERS)
  else
    !------------------------------------sec'species, USING MOLE NUMBERS
    do iAs=1,nAs
      !
      X = LnActW *tNuAs(iAs,isW) &
      & - vDGapp_As(iAs)
      do iAx=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+iAx)) *tNuAs(iAs,nCi+iAx)
      end do
      !
      X= exp(X)
      !
      do iCi=1,nCi
        if(iCi/=isW .and. tNuAs(iAs,iCi)/= 0.0D0) X= X *vX(iCi)**tNuAs(iAs,iCi)
        ! & X= X *(vX(iCi)/vX(isW)/MWSv)**tNuAs(iAs,iCi)
      end do
      !! X= X /(vX(isW)*MWSv)**vSumNuAs(iAs)
      !! vMolal(iAs)= X
      vMolSec(iAs)= X *(vX(isW)*MWSv)**(One -vSumNuAs(iAs))
      !
      !------------------------ Psi(nCi+iAs)-  equilibrium constraint --
      !! if(.not. DirectSub) vFunc(nCi+iAs)= X *(vX(isW) *MWSv) - vX(nCi+iAs)
      if(.not. DirectSub) vFunc(nCi+iAs)= vMolSec(iAs) - vX(nCi+iAs)
      !
    end do
    !-----------------------------------/sec'species, USING MOLE NUMBERS
  end if
  !-------------------------/chemical equilibrium for Second'Aqu.Species
  !
  !------------------------------------------------material conservation
  if(LogForAqu) then
    !-----------------------------------material conservation, LogForAqu
    do iCi=1,nCi
      !
      !----------------------------------contrib' from primary species--
      X= sum(tAlfPr(iCi,1:nCi)*vXf(1:nCi))
      !--/
      do iAs=1,nAs
        if(tAlfAs(iCi,iAs)/= 0.0D0) then
          if(DirectSub) then  ;  X= X + tAlfAs(iCi,iAs)*vMolSec(iAs)
          else                ;  X= X + tAlfAs(iCi,iAs)*vXf(nCi+iAs)
          end if
        end if
      end do
      !---------------------contrib' from the aqueous "mobile" species--
      do I=1,nAx
        X= X + tAlfPr(iCi,nCi+I) *vMolal(vOrdPr(nCi+I)) *vXf(isW) *MWSv
      end do
      !
      vFunc(iCi)= X - vTotS(iCi)
      !
    end do
    !
    !! if(iDebug>2) then
    !!   write(71,'(A)') "==RESIDUAL============"
    !!   do I=1,nAs
    !!     write(71,'(A,G15.6)') "secc",vMolSec(I)
    !!   end do
    !!   do I=1,nCi
    !!     write(71,'(A,G15.6)') "prim",vXf(I)
    !!   end do
    !!   do I=1,nCi
    !!     write(71,'(A,G15.6)') "vFunc",vFunc(I)
    !!   end do
    !! end if
    !! !
    !! if(bMolal) vFunc(isW)= vXf(isW) - One/MWSv
    !----------------------------------/material conservation, LogForAqu
  !
  else
    !-------------------------------material conservation, NOT LogForAqu
    do iCi=1,nCi
      !
      !------------------------------------contrib' from primary species
      X= sum(tAlfPr(iCi,1:nCi)*vX(1:nCi))
      !
      !----------------------------------contrib' from secondary species
      do iAs=1,nAs
        if(tAlfAs(iCi,iAs)/= 0.0D0) then
          if(DirectSub) then  ;  X= X + tAlfAs(iCi,iAs) *vMolSec(iAs)
          else                ;  X= X + tAlfAs(iCi,iAs) *vX(nCi+iAs)
          end if
        end if
      end do
      !----------------------contrib' from the aqueous "mobile" species
      do iAx=1,nAx
        X= X + tAlfPr(iCi,nCi+iAx)*vMolal(vOrdPr(nCi+iAx)) *vX(isW)*MWsv
      end do
      !
      vFunc(iCi)= X -vTotS(iCi)
      ! in SPC, vTotS = vTot_inFluid
      ! in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
      !
    end do
    !
    if( DirectSub .and. iDebug>2) then
      write(71,'(A)') "==RESIDUAL============"
      do iAs=1,nAs
        !! write(71,'(A,G15.6)') "secc",vMolal(iAs)*(vX(isW)*MWSv)
        write(71,'(A,G15.6)') "secc",vMolSec(iAs)
      end do
      do iCi=1,nCi
        write(71,'(A,G15.6)') "prim",vX(iCi)
      end do
      do iCi=1,nCi
        write(71,'(A,G15.6)') "vFunc",vFunc(iCi)
      end do
    end if
    !------------------------------/material conservation, NOT LogForAqu
  end if
  !----------------------------------------------/material conservation
  !
  !----------------------in case of "global equilibrium" calculation, --
  !------------------------------------contribution from non aqu'species
  if(nM>0) then
    select case(cEquMode)
    !
    ! cEquMode is flag for the type of "equilibrium speciation", i.e.
    ! EQ1= reaction progress method, titration-like
    ! EQ2= direct method
    !
    !-------------------------------------------reaction progress method
    case("EQ1") !
      do I=1,nEquFas !nFs
        AffScale= sum(abs(tNuEq(I,:)))
        !--------------------------------------compute phase affinities
        X=   vDeltaG_Eq(I)    & !
        &  - tNuEq(I,1) *vLnAct(isW)  &                    ! solvent
        &  + sum(tNuEq(I,2:nCi)) *(vLnXf(isW) +log(MWSv))  ! mole nr to molality
        do J=2,nCi
          X= X - ( vLnXf(J) +vLnGam(vOrdAq(J)) )*tNuEq(I,J)
        end do
        !----------------------------- add terms for mobile species --
        do J=1,nAx+nMx
          X= X -tNuEq(I,nCi+J) *vLnAct(vOrdPr(nCi+J))
        end do
        !
        vFasAff(I)= X
        !--------------------------------------/compute phase affinities
        !
        !------------ mineral production rate (precip' - vFasRat>0) --
        !----------------------- rate - exp(Affinity) - 1 - Q/K - 1 --
        vFasRat(I)= FSafe_Exp(-vFasAff(I) /AffScale) - One
        !
        !------- add mineral production in mass balance constraints --
        do iCi=1,nCi !solvent involved in mass balance
        !do iCi=2,nCi  !solvent not involved in mass balance
          vFunc(iCi)= vFunc(iCi) &
          &         + tAlfEq(iCi,I) &
          &           *(vFasMole(I) + vFasRat(I) *dXI) !XI= progress param'
        end do
        !
        !--- add material balance equations --
        !--- Rate(xi) - [vFasMole(xi+dxi) - vFasMole(xi)] / dxi
        vFunc(nF+I)= vX(nF+I) - vFasMole(I) - dXi *vFasRat(I)
        !
      end do ! do iFs
      !-------------------------------------------/loop on active phases
    !------------------------------------------/reaction progress method
    !
    case("EQ2") !
      !
      if(nEquFas>0) then
        do I=1,nEquFas
          !--------------------------------------------POTENTIAL BALANCE
          !-- add potential balance equations imposed by non-aqu'phase
          !-- i.e. impose vFasAff(I)-Zero
          X=  vDeltaG_Eq(I) &
          -             tNuEq(I,isW)    *vLnAct(isW)            & ! solvent
          - dot_product(tNuEq(I,2:nCi),  vLnXf(2:nCi)        )  & !- vX(nCi+I)
          +         sum(tNuEq(I,2:nCi))*(vLnXf(isW) +log(MWSv)) & !mole nr -> molality
          - dot_product(tNuEq(I,2:nCi),  vLnGam(vOrdAq(2:nCi)))
          !
          do J=1,nAx ! aqu' mobile prim' species
            X= X - tNuEq(I,nCi+J) *vLnAct(vOrdPr(nCi+J))
          end do
          do J=1,nMx ! min' mobile prim' species
            X= X - tNuEq(I,nCi+nAx+J) *vLnAct(vOrdPr(nCi+nAx+J))
          end do
          !
          if(Complementarity) then
            ! vFunc(nF+I)= min(X,vX(nF+I))
            !! if(X < vX(nF+I)) then  ;  vFunc(nF+I)= X
            !! else                   ;  vFunc(nF+I)= vX(nF+I)
            !! end if
            ! vFunc(nF+I)= Fischer - Burmeister
            vFunc(nF+I)= X + vX(nF+I) - sqrt(X**2 + vX(nF+I)**2)
          else
            vFunc(nF+I)= X
          end if
          !---/
          !
          !---------------------------------------------MATERIAL BALANCE
          !-- add contribution of non aqu'phase I on material balance
          !-- equations
          do iCi=1,nCi
            vFunc(iCi)= vFunc(iCi) + tAlfEq(iCi,I) *vX(nF+I)
          end do
          !---/
        end do
      end if
    !-----------------------------------------------------/direct method
    !
    end select
    !
  end if ! if(nM>0)
  !-----------------------------------/contribution from non aqu'species
  !
  deallocate(vXf)
  deallocate(vLnXf)
  deallocate(vSumNuAs)
  if(nM>0) deallocate(vXm)
  !
  Equil_Residual= vFunc
  !
  return
end function Equil_Residual

end module M_Equil_Residual


