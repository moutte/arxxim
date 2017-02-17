module M_Dynam_Residual
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_
  implicit none
  !
  private
  !
  public:: Dynam_Residual
  public:: Dynam_Converge
  public:: vTolCoef,vNewtTolF
  !
  real(dp),allocatable:: vTolCoef(:)
  real(dp),allocatable:: vNewtTolF(:)
  !
  logical:: UseTolCoef= .false.  !! .true. !!

contains

logical function Dynam_Converge(vF,vTolF)
!--
!-- makes possible adaptation of convergence criteria to type of condition
!-- e.g. material conservation vs potential, kinetic, etc.
!--
  real(dp),intent(in):: vF(:),vTolF(:)
  !
  if(UseTolCoef) then
    Dynam_Converge= all(ABS(vF(:)) < vNewtTolF(:) *vTolCoef(:))
  else
    Dynam_Converge= all(ABS(vF(:)) < vNewtTolF(:))
  end if
  !
end function Dynam_Converge

subroutine BuildVecs( & !
& LogForAqu,          & !IN
& LogForMin,          & !IN
& nF,nEqA,nMkA,       & !IN
& vX,                 & !IN
& vXf,vLnXf,vXeq,vXmk)
  !
  use M_Safe_Functions
  use M_Numeric_Const,only: MinExpDP,MaxExpDP,Ln10
  !
  logical, intent(in) :: LogForAqu
  logical, intent(in) :: LogForMin
  integer, intent(in) :: nF,nEqA,nMkA
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vXf(:),vLnXf(:)
  real(dp),intent(out):: vXeq(:),vXmk(:)
  !
  integer:: N
  !
  if(LogForAqu) then
    vXf(1:nF)= FSafe_vExp(vX(1:nF)) !-> avoid overflow
    vLnXf(:)= vX(1:nF)
  else
    vXf(1:nF)= vX(1:nF)
    vLnXf(:)= FSafe_vLog(vX(1:nF))
  end if
  !
  if(nEqA>0) then
    vXeq=Zero
    N= nF
    if(LogForMin) then
      != works with log(mineral amount)
      vXeq(1:nEqA)= FSafe_vExp(vX(N+1:N+nEqA))
    else
      != works directly with mineral amounts
      vXeq(1:nEqA)=vX(N+1:N+nEqA)
    end if
  end if
  !
  if(nMkA>0) then
    vXmk= Zero
    N= nF+nEqA
    if(LogForMin) then
      != works with log(mineral amount)
      vXmk(1:nMkA)= FSafe_vExp(vX(N+1:N+nMkA))
    else
      != works directly with mineral amounts
      vXmk(1:nMkA)=vX(N+1:N+nMkA)
    end if
  end if
  !
end subroutine BuildVecs

function Dynam_Residual(vX) result(vFunc)
!--
!-- system to be solved
!--
  use M_Global_Vars,  only: vKinFas,nAq,vFas
  use M_T_KinFas,     only: KinFas_CalcSurf
  use M_KinFas_Surf,  only: KinFas_Surf_Calc
  use M_KinRate
  !
  use M_Basis_Vars, only: MWSv,isW,nCi,nAs,nAx,nMx
  use M_Basis_Vars, only: tAlfPr,tAlfAs,tNuAs
  use M_Dynam_Vars, only: LogForAqu, LogForMin, DirectSub, CoupledCoores
  use M_Dynam_Vars, only: PhiF,PhiInert,VBox,UDarcy,dX,dTime
  use M_Dynam_Vars, only: vQTotInj,vTotInj,vTotF,Fout
  use M_Dynam_Vars, only: vLnGam,vLnBuf,vLnAct,vMolSec
  use M_Dynam_Vars, only: vLKinActiv,vLEquActiv
  use M_Dynam_Vars, only: vKinMod,vMolarVol
  use M_Dynam_Vars, only: vMolK,vSurfK,vMolK0,vSurfK0
  use M_Dynam_Vars, only: vVmQsK,vVmAct
  use M_Dynam_Vars, only: tNu_Kin,tAlfKin,vDG_Kin,vDG_As,vDLGam_As
  use M_Dynam_Vars, only: nEvalFunc
  !
  use M_Dynam_Solve_Vars
  !
  real(dp),dimension(:),intent(in):: vX
  !
  real(dp),dimension(:),allocatable:: vXf,vLnXf,vXmk,vXeq,vVm,vMolal,vSumNuAs
  !
  real(dp),dimension(1:size(vX)):: vFunc
  !
  real(dp):: QTotInjPr,X,TotFF
  integer :: I,J,iPsi,iPr,iAs,iMk,iMkA,nF,nMk,nMkA,nEqA,nCx
  !
  !-----------------------------------------------under construction ...
  if(nAx>0) then
    allocate(vMolal(nAx))
    !update molalities of aqueous "mobile" species (i.e. buffered)
    !-> used to calculate mole number= molality * solvent'mass
    !! vMolal(1:nAx)= &
    !! & exp(vLnAct(vOrdPr(nCi+1:nCi+nAx))-vLnGam(vOrdPr(nCi+1:nCi+nAx)))
    !-> molality
  end if
  !----------------------------------------------/under construction ...
  !
  if(iDebug>2) nEvalFunc= nEvalFunc +1
  !
  nCx=  nAx+nMx
  nMk=  size(vKinFas)
  nMkA= count(vLKinActiv)
  nEqA= count(vLEquActiv)
  !
  if(DirectSub) then  ;  nF= nCi
  else                ;  nF= nAq
  end if
  !
  allocate(vXf(nF))
  allocate(vLnXf(nF))
  allocate(vXmk(nMkA))
  allocate(vXeq(nEqA))
  allocate(vVm(nMk))
  !
  call BuildVecs( & !
  & LogForAqu,LogForMin,nF,nEqA,nMkA,vX, & !
  & vXf,vLnXf,vXeq,vXmk)
  !
  call KinRate_Calc( & !
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
  call PhiF_Calc( & !
  & PhiInert,              & !IN
  & vLKinActiv,vLEquActiv, & !IN
  & vXmk,vXeq,             & !IN
  & PhiF)                    !OUT
  !
  !~ if(LogForAqu) then
  allocate(vSumNuAs(nAs))
  do iAs=1,nAs
    X= Zero
    do I=1,nCi
      if(I/=isW) X= X+tNuAs(iAs,I)
    end do
    vSumNuAs(iAs)= X
  end do
  !~ end if
  !
  !--------------------------------FINITE DifFERENCE SYSTEM TO BE SOLVED
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
  if(LogForAqu) then
    !----------------------------------------------------------LogForAqu
    do iAs=1,nAs
      X=  vLnAct(1) *tNuAs(iAs,1) ! Solvent
      do J=1,nCi
        if(J/=1) X= X + vX(J) *tNuAs(iAs,J) ! Prim'Solute species
      end do
      X= X - (vSumNuAs(iAs) -One)*(vX(1) +log(MWSv)) & ! mole nr -> molality
      &    - vDG_As(iAs)     & !
      &    - vDLGam_As(iAs)
      !
      !----------------"buffered" min'species -> add Sum(Nu_iMx*Act_iMx)
      do J=1,nCx
        X=  X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      end do
      !
      if(DirectSub) then
        !--------------------------------------------direct substitution
        vMolSec(iAs)= exp(X)
      else
        ! Psi(nCi+1:nAq)=  equilibrium constraints (ChemicalAffinity=0)
        vFunc(nCi+iAs)= X - vX(nCi+iAs)
        vTolCoef(nCi+iAs)= MAX(One,ABS(vDG_As(iAs)+vDLGam_As(iAs)))
        !if(iDebug>2) write(73,'(I3,1X,G15.6)') nCi+iAs, vTolCoef(nCi+iAs)
      end if
      !
    end do
    !---------------------------------------------------------/LogForAqu
  else
    !------------------------------------------------------NOT LogForAqu
    !~ print *,"SIZE(vMolSec)=",size(vMolSec)  ;  pause
    do iAs=1,nAs
      !
      X= vLnAct(1)*tNuAs(iAs,1) -vDG_As(iAs) -vDLGam_As(iAs)
      !
      !------------ "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
      do J=1,nCx
        X= X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      end do
      !
      X= exp(X)
      !
      do iPr=1,nCi
        if(iPr /= 1) then
          if(tNuAs(iAs,iPr)/=Zero) X= X *vX(iPr)**tNuAs(iAs,iPr)
          !& X= X *(vX(iPr)/vX(1)/MWSv)**tNuAs(iAs,iPr)
        end if
      end do
      vMolSec(iAs)= X *(vX(1)*MWSv)**(One -vSumNuAs(iAs))
      !~ print *,"iAs,vMolSec(iAs)=",iAs," _ ",vMolSec(iAs)  ;  pause
      !
      if(.not. DirectSub) vFunc(nCi+iAs)= vMolSec(iAs) - vX(nCi+iAs)
      !
    end do
    !-----------------------------------------------------/NOT LogForAqu
  end if
  !--------------------------------------------------/Second'Aqu.Species
  !
  !-----------------------------------------------component conservation
  iPsi= 0
  !--------------------------------------------------------- iPsi- 1:nCi
  do iPr=1,nCi
    !
    !---------------------------------------------mole nr injected /time
    if (CoupledCoores) then
      QTotInjPr= vQTotInj(iPr) /Vbox
      ! mol.s-1 molar rate (explicit rates)
      ! Fout given
    else
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
    end if
    !--------------------------------------------/mole nr injected /time
    !
    !! InFluid= dot_product(tAlfPr(iPr,1:nCi), vXf(1:nCi))          & ! # moles in fluid, from prim'species
    !! &      + dot_product(tAlfAs(iPr,1:nAs), vXf(nCi+1:nCi+nAs))    ! # moles in fluid
    !! InMin=   dot_product(tAlfKin(iPr,1:nMk),vVm(1:nMk)) * dTime    ! to/from minerals
    !
    !-------------------------------mole nrs in fluid, from prim'species
    TotFF= dot_product(tAlfPr(iPr,1:nCi),vXf(1:nCi))
    !
    !-------------------------------moles nrs in fluid, from sec'species
    do iAs=1,nAs
      if(DirectSub) then
        TotFF= TotFF + vMolSec(iAs) *tAlfAs(iPr,iAs)
      else
        !TotFF= TotFF + exp(vX(nCi+iAs)) *tAlfAs(iPr,iAs)
        if(LogForAqu) then  ;  TotFF= TotFF + vXf(nCi+iAs) *tAlfAs(iPr,iAs)
        else                ;  TotFF= TotFF + vX(nCi+iAs)  *tAlfAs(iPr,iAs)
        end if
      end if
    end do
    !------------------------------/moles nrs in fluid, from sec'species
    !
    !--------------moles nrs in fluid, from the aqueous "mobile" species
    do J=1,nAx
      TotFF= TotFF + MWSv *vXf(isW) *vMolal(J) *tAlfPr(iPr,nCi+J)
    end do
    !---/
    !
    X=  TotFF *(One + dTime *Fout /Vbox) & ! Ci = (ni*Vbox )/ (Vbox*PhiW) = ni/Phiw
    & - vTotF(iPr)                       & ! # moles in fluid at previous step
    & - dTime *QTotInjPr                   ! # moles injection
    !
    !-------------------------------------------------------equil'phases
    if(nEqA>0) then
      J=0
      do iMk=1,nMk
        if(vLEquActiv(iMk)) then
          J=J+1
          X= X + tAlfKin(iPr,iMk)*(vXeq(J) -vMolK(iMk))
        end if
      end do
    end if
    !------------------------------------------------------/equil'phases
    !
    !---------------------------------------------------------kin'phases
    if(nMkA>0) then
      do iMk=1,nMk
        if(vLKinActiv(iMk)) X= X + tAlfKin(iPr,iMk) *dTime *vVm(iMk)
      end do
    end if
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
  end do
  !----------------------------------------------/component conservation
  !
  if(iDebug>2 .and. DirectSub) then
    write(71,'(A)') "==RESIDUAL============"
    do I=1,nAs
      write(71,'(A,G15.6)') "secc",vMolSec(I)
    end do
    do I=1,nCi
      write(71,'(A,G15.6)') "prim",vXf(I)
    end do
    do I=1,nCi
      write(71,'(A,G15.6)') "vFunc",vFunc(I)
    end do
  end if
  !
  if(.not. DirectSub) iPsi= iPsi +nAs
  !
  !---------------------------------------------------equilibrium phases
  !---------------------------------------------------iPsi- nF+1:nF+nEqA
  if(nEqA>0) then
    do iMk=1,nMk
      if(vLEquActiv(iMk)) then
        !
        !--- add one equation for each equ'phase
        X=   vLnAct(isW) *tNu_Kin(iMk,isW)                     & ! Solvent activity
        &  - SUM(tNu_Kin(iMk,2:nCi))*(vLnXf(isW) +log(MWSv))   & ! Solvent term (mole nr -> molality
        &  + dot_product(tNu_Kin(iMk,2:nCi),vLnGam(2:nCi))     & ! gamma prim'solute species
        &  - vDG_Kin(iMk)
        do J=1,nCi
          if(J/=isW) X= X + vLnXf(J) *tNu_Kin(iMk,J)
        end do
        !---------- "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
        do J=1,nCx
          X= X + vLnBuf(J) *tNu_Kin(iMk,nCi+J)
        end do
        !---/
        !
        iPsi= iPsi+1
        vFunc(iPsi)= X
        !
        ! vFunc(iPsi)= &
        ! &             tNu_Kin(iMk,  isW) * vLnAct(isW)          & ! Solvent
        ! + dot_product(tNu_Kin(iMk,2:nCi), vX(2:nCi))            & ! prim' solute species
        ! + dot_product(tNu_Kin(iMk,2:nCi), vLnGam(2:nCi))        & ! gamma prim'solute species
        ! -         SUM(tNu_Kin(iMk,2:nCi))*(vX(isW) +log(MWSv))  & ! Solvent term (mole nr -> molality
        ! - vDG_Kin(iMk)
        ! !
        ! !---------- "buffered" min'species -> add Sum(Nu_iMx*Act_iMx) --
        ! if(nCx>0) &
        ! & vFunc(iPsi)= vFunc(iPsi) &
        ! + dot_product(tNu_Kin(iMk,nCi+1:nCi+nCx),vLnBuf(1:nCx))
      end if
    end do
  end if
  !--------------------------------------------------/equilibrium phases
  !
  !-------------------------------------------------------kinetic phases
  !---------------------------- iPsi- nCi+nAs+nEqA+1 : nCi+nAs+nEqA+nMkA
  if(nMkA>0) then
    iMkA= 0
    do iMk=1,nMk
      if(vLKinActiv(iMk)) then
        iMkA= iMkA + 1
        !--- add one equation for each kin'phase
        iPsi= iPsi+1
        vFunc(iPsi)= vXmk(iMkA) -vMolK(iMk) -dTime *vVm(iMk)
        vTolCoef(iPsi)= MAX(One,ABS(vMolK(iMk)))
        if(iDebug>2) write(73,'(I3,1X,G15.6)') iPsi, vTolCoef(iPsi)
        !---/
      end if
    end do
  end if
  !------------------------------------------------------/kinetic phases
  !
  if(nAx>0) deallocate(vMolal)
  !
  deallocate(vXf)
  deallocate(vLnXf)
  deallocate(vVm)
  deallocate(vXmk)
  deallocate(vXeq)
  deallocate(vSumNuAs)
  !
  !do i=1,size(vFunc)
  !  write(81,'(G12.3,1X)',advance="no") vFunc(I)
  !end do
  !write(81,*)
  !Dynam_Residual= vFunc
  !
  return
end function Dynam_Residual

subroutine KinRate_Calc( & !
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
  use M_T_KinFas,   only: KinFas_CalcSurf
  use M_KinFas_Surf,only: KinFas_Surf_Calc
  use M_KinRate
  !
  use M_Global_Vars,only: vKinFas,nAq,vFas
  use M_Basis_Vars, only: MWSv,isW,nCi,nAx,nMx
  use M_Basis_Vars, only: tAlfPr,tAlfAs,tNuAs
  !
  use M_Dynam_Vars, only:  & !
  & Implicit_Surface,      & !
  & Implicit_ActivFactor,  & !
  & QsK_Iota, VBox,        & !
  & vLnGam,vLnBuf,vLnAct,  & !
  & vKinMod,vMolarVol,vKinMinim, & !
  & vMolK0,vSurfK0,        & !
  & tNu_Kin,tAlfKin,       & !
  & vDG_Kin,vDG_As
  !
  real(dp),intent(in) :: dTime
  logical, intent(in) :: vLKinActiv(:)
  real(dp),intent(in) :: vX(:),vXmk(:)
  real(dp),intent(inout):: vSurfK(:),vVmQsK(:),vVmAct(:)
  real(dp),intent(out):: dSRdXm(:,:),dSRdLnXm(:,:)
  real(dp),intent(out):: dVmQdLnXf(:,:)
  real(dp),intent(out):: dVMdLnXf(:,:),dVMdLnXm(:,:)
  !
  real(dp),dimension(1:nCi):: dQsKdLnXi
  real(dp):: QsK,QskIota,varMaxM,X
  integer :: nMkA,nMk,nCx
  integer :: iPr,iMk,iMkA
  character(len=7):: cSatur
  !
  nCx=  nAx+nMx
  nMk= size(vKinFas)
  nMkA= count(vLKinActiv)
  !
  vVmQsK=   Zero
  dSRdLnXm= Zero
  dSRdXm=   Zero
  dVmQdLnXf=Zero
  dVMdLnXf(:,1:nCi)= Zero
  dVMdLnXm(:,1:nMk)= Zero
  !
  !if(Implicit_ActivFactor) then !update activities
  !  !Nota: as Activity Coeffs are not Implicited, Rate_Act is better not Implicited
  !  dVmAdLnXf=Zero
  !  vLnAct(2:nAq)=vX(2:nAq)+vLnGam(2:nAq)-vX(isW)-log(MWSv)
  !end if
  !
  QskIota= QsK_Iota !Zero !mod 10/06/2008 17:02
  !
  ! if(nMkA + nEqA >0) then
  if(nMkA > 0) then
    !
    varMaxM= Zero
    iMkA= 0
    iMk=  0
    !
    do
      !
      iMk= iMk+1
      if(iMk>size(vKinFas)) exit
      !
      !--------------------------------------------iMk is kinetic active
      if(vLKinActiv(iMk)) then
        !
        iMkA= iMkA + 1
        !
        !--------------------------------------------------calc. surface
        if(Implicit_Surface) then
          !X= vXmk(iMkA)
          X= MAX(vXmk(iMkA),vKinMinim(iMkA))
          call KinFas_Surf_Calc( & !
          & .true.,             & !IN: bImplicit
          & vKinFas(iMk)%cMode, & !IN: kinetic model
          & X,                  & !IN: Nr Moles of mineral
          & vMolK0(iMk),        & !IN: Nr Moles of phase at given state
          & vSurfK0(iMk),       & !IN: Surface of phase at that state
          & vSurfK(iMk),        & !OUT: Surface
          & dSRdLnXm(iMk,iMk),  & !OUT:
          & dSRdXm(iMk,iMk))      !OUT:
        end if
        !-------------------------------------------------/calc. surface
        !
        !----------------------------------calc mineral rate per surface
        !--------------------------------------(- only "chemical" terms)
        call KinRate_CalcQsK( & !
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
        call KinRate_SatState( &
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
        !! if(cSatur/="MINIMAL") & !
        call KinRate_CalcQsKFactor( &
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
        !!!  if(Implicit_ActivFactor) & !else will use current values obtained in Dynam_OneStep_KinFasUpd
        !!!  & call KinRate_CalcActivFactor(&
        !!!  & vKinFas(iMk)%cSatur, & !IN
        !!!  & vKinMod(vKinFas(iMk)%iKin), & !IN: kinetic parameters
        !!!  & vLnAct,           & !LnActH_,LnActOH,LnActCO2, & !activators
        !!!  & vVmAct(iMk),      & !OUT
        !!!  & dVmAdLnXf(iMk,:))   !OUT
        !
        !! if(TestMax .and. vMolK(iMk)>1.0D-6) &
        !! & varMaxM=MAX(varMaxM, ABS( dTime*vVm(iMk) / vMolK(iMk) ))
        !
      end if
      !-------------------------------------------/iMk is kinetic active
      !
    end do
    !
    !! if(TestMax .and. varMaxM>1.0D0) then
    !!   dTime=dTime * 0.5_dp
    !!   if(iDebug>0) write(fTrc,'(A)') "dtime adjusted in vFunc"
    !! end if
  end if
  !
  !----------------------------------------------------/CalcMineralRates
end subroutine KinRate_Calc

subroutine PhiF_Calc( & !
!--
!-- compute porosity
!--
& PhiInert,              & !IN
& vLKinActiv,vLEquActiv, & !IN
& vXmk,vXeq,             & !IN
& PhiF)                    !OUT
  !
  use M_T_KinFas,   only: KinFas_CalcSurf
  use M_KinFas_Surf,only: KinFas_Surf_Calc
  use M_KinRate
  !
  use M_Global_Vars,only: vKinFas
  use M_Dynam_Vars, only: VBox, vMolarVol
  !
  real(dp),intent(in) :: PhiInert
  logical, intent(in) :: vLKinActiv(:),vLEquActiv(:)
  real(dp),intent(in) :: vXmk(:),vXeq(:)
  real(dp),intent(out):: PhiF
  !
  integer :: nMk
  integer :: iMk,iMkA,iEqA
  !
  nMk= size(vKinFas)
  !
  PhiF= One - PhiInert
  !-> take account of minerals that are "inactive" during current time step
  !
  if(count(vLKinActiv) + count(vLEquActiv) >0) then
    !
    iMkA= 0
    iEqA= 0
    iMk=  0
    !
    do
      iMk= iMk+1
      if(iMk>nMk) exit
      !
      !---------------------------------------- iMk is kinetic active --
      if(vLKinActiv(iMk)) then
        iMkA= iMkA + 1
        PhiF= PhiF - vXmk(iMkA) * vMolarVol(iMk) /VBox
      end if
      !---------------------------------------/ iMk is kinetic active --
      !
      !------------------------------------ iMk is equilibrium active --
      if(vLEquActiv(iMk)) then
        iEqA= iEqA +1
        PhiF= PhiF - vXeq(iEqA) * vMolarVol(iMk) /VBox
      end if
      !-----------------------------------/ iMk is equilibrium active --
      !
    end do
    !
  end if
  !
end subroutine PhiF_Calc

end module M_Dynam_Residual
