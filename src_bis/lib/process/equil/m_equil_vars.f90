module M_Equil_Vars
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  use M_Kinds
  use M_T_MixModel,only: MaxPole
  implicit none
  !
  public
  !
  real(dp),dimension(:),allocatable:: &
  & vTotS,    & !1:nCp=  total amounts of COMPONENTS (incl. charge balance) in the SYSTEM at time T
  ! in SPC, vTotS = vTot_inFluid
  ! in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
  & vMolF,    & !1:nAq=  nr moles of aqu'species in the system, e.g. for 1 kg H2O
  & vMolal,   & !1:nAq=  molalities
  & vMolSec,  & !1:nAs= mole numbers second species, for DirectSub
  & vFasMole, & !1:nFs= mole nrs of phases (added for SpeciaEquil) -> saved to vFas(:)%Mole
  & vLnAct,   & !
  & vLnGam,   & !
  & vDeltaG_As, & ! DeltaGibbsFreeEnergy of formation reaction, for Second'Aqu'Species
  & vDGapp_As,  & ! "apparent deltaG": includes delta(log(gammas solute species))
  & vDeltaG_Fs    ! DeltaGibbsFreeEnergy of formation reaction, for phases of fixed composition
  !
  logical, dimension(:),allocatable:: vYesList != the phase is taken into account
  !~ & vFasYes     !vFasYes= the phase is "active"
  !
  ! variables added for Specia_CalcEquil_N
  ! affinity of minerals or gases, scaling factor, rate (for cEquMode=3)
  real(dp),dimension(:),allocatable:: vFasAff,vFasRat,vAffScale 
  !
  !--- added for M_Equil_3 / M_Equil_3_Mix
  real(dp),allocatable:: &
  ! for equilibrium non'aqu'phases,
  & vDeltaG_Eq(:), &  ! DeltaGibbsFreeEnergy of formation reaction
  & tAlfEq(:,:),   &  !
  & tNuEq(:,:)        !
  ! Relation betweeen (vDeltaG_Eq,tAlfEq,tNuEq)
  !               and (vDeltaG_Fs,tAlfFs,tNuFas)
  !~ do iFs=1,size(vFasYes)
    !~ if(vFasYes(iFs)) then
      !~ I=I+1
      !~ vDeltaG_Eq(I)= vDeltaG_Fs(iFs)
      !~ tAlfEq(:,I)= tAlfFs(:,iFs)
      !~ tNuEq(I,:)= tNuFas(iFs,:)
    !~ end if
  !~ enddo
  !
  type:: T_EquPhase
    character(len=23):: NamEq
    integer :: iPur
    integer :: iMix
    integer :: nPole
    integer :: vIPole(1:MaxPole)
    real(dp):: vXPole(1:MaxPole)
    real(dp):: vLnAct(1:MaxPole)
    !real(dp):: Grt
    real(dp):: Mole
  end type T_EquPhase
  !
  type(T_EquPhase),allocatable:: vEquFas(:)
  integer :: nEquFas
  !---/
  !
  !~ real(dp):: OsmoSv
  real(dp):: dXi
  !
  logical,dimension(:),allocatable:: vTooLow !,vTooHih
  !
  integer:: nEvalFunc
  !
  logical:: LogForAqu= .true.
  !
  logical:: Complementarity= .false.  !! .true.  !!
  !
  logical:: DirectSub= .false.
  ! DirectSub is "direct substitution":
  ! second'species mole numbers substituted in material consevation equations
  !
  character(len=30):: cMethod= "NEWTONPRESS"
  !
  integer:: & !
  & fTrAct=  0, & !"trace" file for activity calculations, opened if DebActiv
  & fActiz=  0, & !log file for activities
  & fTrcEq=  0, & !
  & iNewt=   0, & !
  & iCountEq=0
  !
  logical:: & 
  & DebNewt=   .false., &
  & DebJacob=  .false., &
  & bFinDif=   .false., &
  & DebActiv=  .false., &
  & DebFormula=.false.
  !
  ! variables for Equil_EqN
  character(len=3):: cEquMode= "EQ2"
  ! flag for the type of "equilibrium speciation"
  !   EQ1= step method, titration-like
  !   EQ2= direct method, suitable for non-singular phase assemblage, given a priori
  !
  !~ logical:: UseOsmoSv= .false. !.true. !
  ! the system (residual=zero) is solved at constant activity coeff's of solute
  ! currently, concerning the solvent, the system can be solved
  !   (1) at constant ln(Activ(solvent))  (= was the only option available until 2008)
  !or (2) at constant osmotic coefficient (= new option, under test)
  !
  !~ logical:: bMolal=  .false. !.true. !
  !
  logical:: bDirect= .false. !.true. !
  !
  !--- parameters for convergence on gammas --
  integer, parameter:: GamMaxIts=100    !maximum number of iterations on gamma's
  real(dp),parameter:: TolGam=   1.0D-6 !convergence condition on delta(Log(Gamma))
  real(dp),parameter:: TolAct=   1.0D-6 !convergence condition on delta(Log(Activ))
  !---/
  !
  !--- parameters for Newton convergence --
  integer, parameter :: NewtMaxIts=200 !maximum number of iterations !!!!!!200
  real(dp),parameter :: &
  !& NewtTolF=   1.0E-03_dp,& !convergence criterion on function values
  !& NewtTolMin= 1.0E-06_dp,& !criterion for spurious convergence
  & NewtTolF=   1.0E-09_dp, &  !convergence criterion on function values
  & NewtTolMin= 1.0E-09_dp, &  !criterion for spurious convergence
  !!& NewtTolX=   1.0E-12_dp
  & NewtTolX=   2.0E-16_dp !EPSILON(real(dp))
  ! convergence criterion on dx -> approx. 2.E-16
  ! -> NewtTolX is used to check whether Newton loops on a stationary point
  !---/
  !
contains

subroutine Equil_Vars_Alloc(nCp,nSp,nAq,nAs,nFs)
  integer,intent(in):: nCp,nSp,nAq,nAs,nFs
  !
  call Equil_Vars_Clean
  !
  allocate(vTotS(1:nCp))       ;  vTotS= Zero
  !
  allocate(vTooLow(1:nSp))     ;  vTooLow=.false.
  allocate(vLnAct(1:nSp))      ;  vLnAct= Zero 
  allocate(vLnGam(1:nSp))      ;  vLnGam= Zero
  !
  allocate(vMolF(1:nAq))
  allocate(vMolal(1:nAq))
  !
  allocate(vEquFas(1:nCp))
  !
  allocate(vDeltaG_As(1:nAs))
  allocate(vDGapp_As(1:nAs))
  !
  !~ allocate(vFasYes(1:nFs))     ; vFasYes=  .false.
  allocate(vFasMole(1:nFs))    ; vFasMole= Zero
  allocate(vYesList(1:nFs))    ; vYesList= .false.
  !
  allocate(vDeltaG_Fs(1:nFs))  ; vDeltaG_Fs= Zero 
  allocate(vFasAff(1:nFs))     ; vFasAff=    Zero 
  allocate(vFasRat(1:nFs))     ; vFasRat=    Zero
  allocate(vAffScale(1:nFs))   ; vAffScale=  One
  !
end subroutine Equil_Vars_Alloc

subroutine Equil_Vars_Clean

  if(allocated(vTotS))   deallocate(vTotS)
  if(allocated(vMolF))   deallocate(vMolF)
  if(allocated(vMolal))  deallocate(vMolal)
  if(allocated(vLnAct))  deallocate(vLnAct)
  if(allocated(vLnGam))  deallocate(vLnGam)
  if(allocated(vTooLow)) deallocate(vTooLow)
  !
  if(allocated(vEquFas)) deallocate(vEquFas)
  !
  if(allocated(vDeltaG_As)) deallocate(vDeltaG_As)
  if(allocated(vDeltaG_Fs)) deallocate(vDeltaG_Fs)
  if(allocated(vDGapp_As))  deallocate(vDGapp_As)
  !
  if(allocated(vYesList) ) deallocate(vYesList) 
  !~ if(allocated(vFasYes)  ) deallocate(vFasYes)  
  if(allocated(vFasAff)  ) deallocate(vFasAff)  
  if(allocated(vFasRat)  ) deallocate(vFasRat)  
  if(allocated(vFasMole) ) deallocate(vFasMole)  
  if(allocated(vAffScale)) deallocate(vAffScale)

end subroutine Equil_Vars_Clean

end module M_Equil_Vars

