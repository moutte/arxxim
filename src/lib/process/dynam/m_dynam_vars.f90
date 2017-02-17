module M_Dynam_Vars
  use M_Kinds
  use M_T_Component,only: T_Component
  use M_T_KinModel, only: T_KinModel
  use M_T_Phase,    only: T_Phase
  
  implicit none
  
  public
  
  logical:: CoupledCoores= .false. !-> coores
  !
  ! if(AdjustVolBox) then all calculations done for VolBox= VBox0
  ! else calcualtions done with user's VolBox
  logical:: AdjustVolBox=  .true.
  !
  !!<new 200912>
  logical,parameter:: TimeIsSeconds= .false. ! .true.  !
  real(dp):: TimeScale
  !!</new 200912>
  !
  logical,parameter:: OnlyActiv= .false.
  !
  type(T_Component),allocatable:: vCpnInj(:)
  type(T_Component),allocatable:: vCpnBox(:)
  !
  real(dp),dimension(:),allocatable:: & !
  & vQTotInj !1:nCp, flow rates of components by time unit (mol.s-1)

  real(dp),dimension(:),allocatable:: & !
  & vTotF,   &  !1:nCp, total amounts of components (incl. redox balance) in the FLUID at time T
  & vTotInj, &  !1:nCp, nr moles of components in given amount of injected fluid
  & vMolSec, &  !1:nAs, nr moles of aqu'second'species
  & vMolF       !1:nAq, nr moles of aqu'species in the box
  !
  real(dp),dimension(:),allocatable:: & !
  ! used in DynamCalc, Residual
  & vLnAct, & !Ln(Activity) aqu'species in box
  & vLnGam, & !Ln(Activ'Coeff) aqu'species in box 
  & vLnBuf    !Ln(Activity) buffering species
  
  real(dp),dimension(:),allocatable:: & !
  ! used in DynamCalc, Residual
  & vDG_As, & !DeltaGibbsFreeEnergy of reaction, for Second'Aqu'Species
  & vDG_Sp, & !DeltaGibbsFreeEnergy of reaction, for all species
  & vDLGam_As
  !
  !-----------------------------------------kinetic phases, allocatables
  character,dimension(:),allocatable:: & !
  & vStatusK       !1:nMk, saturation state
  !
  real(dp),dimension(:),allocatable:: & !
  & vMolK,     &   !1:nMk, nr moles of kin'species in the box
  & vSurfK,    &   !1:nMk, surface of kin'species in the box
  & vKinQsk        !1:nMk, saturation index
  !
  real(dp),dimension(:),allocatable:: & !
  & vKinMinim, &   !1:nMk, minimal nr moles of kin'species in the box
  & vMolarVol, &   !1:nMk, Molar volume (=vFas(vKinFas(1:nMk)%iFas)%V)
  & vDG_Kin        !1:nMk,
  !
  real(dp),dimension(:),allocatable:: & !
  & vVmAct,    &   !1:nMk, kinetic constant * (activators/inhibitors)
  & vVmQsK         !1:nMk, affinity factor (e.g. QsK-1)
  !!& vVm0         !1:nMk, value of vVm at begin time step (end former time step)
  !                !-> vVm(iMk)= vSurfK(iMk) *vVmAct(iMk) *vVmQsK(iMk)
  real(dp),dimension(:),allocatable:: & !
  & vSurfK0,   &   !1:nMk, surface of phase at ref'time
  & vMolK0         !1:nMk, mole number of phase at ref'time
  !
  logical,allocatable:: vLKinActiv(:)
  logical,allocatable:: vLEquActiv(:)
  integer,allocatable:: vKinPrm(:)
  !
  real(dp),allocatable:: tStoikioAqu(:,:) ! stoikio aqu'sp/elements
  real(dp),allocatable:: tStoikioKin(:,:) ! stoikio kin'sp/elements
  !
  real(dp),dimension(:),allocatable:: & !
  & vWeitSp, &  !1:nAq, vSpc(vOrdAq(1:nAq))%WeitKg
  & vWeitCp     !1:nCp, vEle(vOrdCp(1:nCp))%WeitKg
  !
  real(dp),dimension(:,:),allocatable:: & !
  & tNu_Kin,  &   !1:nMd,1:nCp
  & tAlfKin       !1:nCp,1:nMd
  !
  type(T_KinModel),allocatable:: vKinMod(:)
  !----------------------------------------/kinetic phases, allocatables
  !
  !---------------------------------------------------------------------
  real(dp):: QsK_Iota = Zero !1.D-9 !  1.D-6 !
  !
  character(len=7):: sModelSurf= "SPHERE" !"CRUNCH" !"SPHERE",
  !
  real(dp):: PhiInert !volume fraction occupied by inactive minerals
  !
  !------------------default texture parameters for secondary minerals--
  real(dp):: RadiusMinim= 1.0D-9  !radius (meter), 1.0D-9= 1 nm
  real(dp):: DensitMinim= 1.0D+18  !number of grains / m^3 fluid volume
  !---------------------------------------------------------------------
  ! vKinMinim(:)= lower limit for mole number in the box
  ! -> to insure that the mineral is always present (potentially)
  ! its mole number cannot decrease below vKinMinim(:)
  !!!! vKinMinim currently not used as constant !!!
  !!!! we use RadiusMinim & DensitMinim to compute VolMinim : 
  !!!! -> cf Dynam_Init_KinFas_Texture
  !
  !=< used in TP_Changed, to check whether TP-dependent param's must be updated
  real(dp):: TdgK0,Pbar0
  real(dp):: TolTdgK= One
  real(dp):: TolPbar= One
  !---/
  !
  integer:: fSavTime= 0 !file for saving results at (nearly) fixed time intervals (dTSav)
  integer:: fSavRate= 0 !file for details on the rate parameters of minerals
  !
  logical:: SteadyState_Stop= .false. !.true.
  !
  !-------------------------------------- box volume is FIXED or FREE --
  !------- FIXED- default, FREE- for closed system simulations mainly --
  logical:: UpdateMassFluid
  !---/
  !
  !---------------------------------------------------------------------
  type:: T_DynBox  ! lumping all box related parameters
    real(dp)::   & !
    & VBox,      & ! volume of the "box", m^3
    & dX,        & ! length of the "box", m
    & Fout,      & ! flow rate output, m^3.s-1
    & PhiF,      & ! fluid fraction, i.e. porosity
    & RhoF,      & ! fluid density
    & UDarcy       ! Darcy velocity of fluid
    integer:: nCell
    logical:: VFixed
  end type T_DynBox
  !
  type(T_DynBox):: DynBox,DynBoxUser
  !
  !---------------------------------------------------------------------
  real(dp)::  &
  & VBox,     & ! volume of the "box", m^3
  & dX,       & ! length of the "box", m
  & Fout        ! flow rate output, m^3.s-1
  !
  ! volume of the "box", m^3, default value
  real(dp):: VBox0= 1.D-3
  !
  real(dp)::  &  !
  & UDarcy,   &  !Darcy velocity of fluid
  & PhiF,     &  !fluid vol'fraction, i.e. porosity (and the gas ?)
  & PhiF0,    &  !idem at ref'time
  & RhoF,     &  !fluid density, kg.m^-3
  & pH_
  !
  !---------------------------------------------------------------------
  type:: T_DynColumn
    integer::  nCell
    integer::  Method
    real(dp):: UDarcy, Disp
    real(dp):: dt, dx, duration
    real(dp):: time_save
  end type T_DynColumn
  !
  type(T_DynColumn):: DynColumn
  !---------------------------------------------------------------------
  type:: T_DynTime
    character(len=6):: TUnit !"SECOND","MINUTE","HOUR","DAY","YEAR"
    real(dp):: TimeFactor
    real(dp)::  & !
    & Time,     & ! current time
    & TFinal,   & ! end of simulation
    & dTime,    & ! current time step
    & dTmin,    & ! minimal time step
    & dTMax,    & ! maximal time step
    & dTSav       ! time laps between two records on x_time.tab
  end type T_DynTime
  !
  type(T_DynTime):: DynTime
  !
  character(len=6):: TUnit !"DAY","YEAR","HOUR","SECOND","MINUTE"
  real(dp):: TimeFactor != the length, in seconds, of current time unit
  !
  real(dp)::  & !
  & Time,     & !current time
  & dTime,    & !current time step
  & TFinal      !end of simulation
  !
  real(dp)::  &
  & dTmin,    & !minimal time step
  & dTMax,    & !maximal time step
  & dTSav       !time laps between two records on x_time.tab
  !
  !---------------------------------------------------------------------
  logical:: Extrapole= .true. ! .false. !
  !
  logical:: LogForAqu=.true. ! .false. !
  !
  logical:: LogForMin=.false.
  !if .true. then works on log(MoleNr) instead of MoleNr for kin'phases
  !
  logical:: DirectSub=.false.
  !
  ! implicitation of rate parameters
  ! CAVEAT: as Activity Coeffs are not Implicited, Rate_Act is not Implicited !!!
  logical:: Implicit_ActivFactor !.true.  !used only in Residual
  logical:: Implicit_Surface     !.false. !used only in residual
  !
  !------------------------------------------- parameters for SOLVER ---
  !
  !integer :: iMethod
  character(len=30):: cMethod
  ! default value given in Dynam_Zero_Numeric (=7 -> Newton_Walker)
  ! can be changed at run time (keyword METHOD in DYNAMIC.NUMERIC)
  ! (read by Dynam_ReadNumeric)
  integer :: iCtrlTime
  real(dp):: &
  & Time_Decrease= 0.5_dp, & !timestep decrease factor
  & Time_Increase= 2.0_dp    !timestep increase factor
  !& Time_Increase= 1.5_dp    !timestep increase factor
  !
  integer:: nEvalFunc
  !
  integer:: &
  & NewtMaxIts,  &
  & NewtIterMax, & !if(iter>NewtIterMax) dT0=dT0*Time_Decrease
  & NewtIterMin    !if(iter<NewtIterMin) dT0=dT0*Time_Increase
  
  real(dp):: &
  & NewtTOLF,       & ! convergence criterion on function values
  & NewtTolF_Equil, & ! idem for aqu'species equilibrium
  & NewtTOLMIN,     & ! criterion for spurious convergence
  & NewtTOLX          ! convergence criterion on dx
  !---------------------------------------------------------------------
  !
  logical:: &
  & DebNewt, &
  & DebJacob, &
  & TestJacob, &
  & TestMax, &
  & bFinDif
  !
  integer,public::   &
  & Dynam_nStep,     &
  & Dynam_nNewtIter, &
  & Dynam_nTotalNewtIter
  !
end module M_Dynam_Vars
