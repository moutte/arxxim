module M_T_MixPhase
  use M_Kinds
  use M_Trace,only: fTrc, iDebug, T_, Stop_
  use M_T_MixModel
  implicit none
  !
  private
  !
  public:: T_MixPhase
  !
  public:: MixPhase_Index
  public:: MixPhase_Zero
  public:: MixPhase_CalcActivs
  public:: MixPhase_NormCompo
  public:: MixPhase_Weit
  public:: MixPhase_GibbsRT
  public:: MixPhase_Constraint
  public:: MixPhase_Volume
  !
  public:: MixPhase_CalcMixing   !used in Dtb_Test_Solution

  type:: T_MixPhase
  ! implements a mixture phase (= symetric model)(= )
  ! includes a pointer to the mixing model,
  ! the composition data, and stores activity data
    !----------------------------------------------------- fixed data --
    character(len=23):: Name
    integer          :: iModel !index of mixing model in vMixModel
    !-----------------------------------------------------/fixed data --
    !
    !-------------------------------------------------- variable data --
    !----------------------------------------------- composition data --
    logical :: vLPole(1:MaxPole)
    ! vLPole(I)= T/F == end member I is present /absent
    real(dp):: vXPole(1:MaxPole)
    ! comp'n in mole fractions of the end members
    !real(dp):: vXAtom(1:MaxAtom) !comp'n in atom fractions, size=Model%NAtom
    !
    !----------------- values at current (T,P,X), for the end-members --
    real(dp),dimension(1:MaxPole):: & !
    & vLnAct,  & !Ln(Activities of end-members 1:Model%nPole)
    & vLIdeal, & !Ln(ideal Activities of end-members 1:Model%nPole)
    & vLGam      !Ln(activ.coeff's)
    !! & vLMarg     !Ln(activ.coeff's related to Margules parameters)
    !
    !------------------------- values at current (T,P), for the phase --
    real(dp):: Grt,H,S,V,Cp
    !--------------------------------------------------/variable data --
    !
  end type T_MixPhase

  type:: T_MixPhaseData
  !-- variable data characteristic of a mixture phase (- symetric model)
  !-- includes composition data, activity data, etc.
    !----------------------------------------------- composition data --
    logical :: vLPole(1:MaxPole) ! end member is present/absent
    real(dp):: vXPole(1:MaxPole) ! comp'n in mole fractions of the end members
    !real(dp):: vXAtom(1:MaxAtom) ! comp'n in atom fractions, size=Model%NAtom
    !
    !----------------- values at current (T,P,X), for the end-members --
    real(dp),dimension(1:MaxPole):: & !
    & vLnAct,  & ! Ln(Activities of end-members 1:Model%nPole)
    & vLIdeal, & ! Ln(ideal Activities of end-members 1:Model%nPole)
    & vLGam,   & ! Ln(activ.coeff's)
    & vLMarg     ! Ln(activ.coeff's related to Margules parameters)
    !
    !------------------------- values at current (T,P), for the phase --
    real(dp):: Grt,H,S,V,Cp
    !
  end type T_MixPhaseData

contains

integer function MixPhase_Index(V,Str)
!--
!-- position of mixture phase named Str in V(1:size(V))
!--
  type(T_MixPhase),intent(in)::V(:)
  character(*),    intent(in)::Str
  !
  integer::I
  !
  MixPhase_Index=0
  if(size(V)==0) return
  !
  I=0
  do
    I=I+1 !; if(iDebug>0) write(fTrc,'(A)') vCpn(I)%SpName
    if(trim(Str)==trim(V(I)%Name)) then
      MixPhase_Index= I
      exit
    end if
    if(I==size(V)) exit
  end do
  !if Str not found, I=0
  !
  return
end function MixPhase_Index

subroutine MixPhase_Zero(S)
  type(T_MixPhase),intent(out):: S
  S%Name=      "Z"
  S%iModel=  0
  !
  S%vLPole=  .false.
  S%vXPole=  Zero
  !S%vXAtom=  Zero
  !~ S%vLnAct=  Zero
  S%vLIdeal= Zero
  S%vLGam=   Zero
  !! S%vLMarg=  Zero
  !
end subroutine MixPhase_Zero

subroutine MixPhase_NormCompo(Mix)
!--
!-- normalize end-member fractions to tot-100%
!-- should not always normalize ???
!--
  type(T_MixPhase),intent(inout):: Mix
  !
  real(dp),dimension(1:MaxPole)::vX
  ! = composition: mole fractions of the end members
  real(dp):: Y,Eps
  !
  Eps= EPSILON(Y)
  !
  where(.not. Mix%vLPole) Mix%vXPole=Zero
  !
  vX(1:MaxPole)=Mix%vXPole(1:MaxPole)
  Y= SUM(vX(1:MaxPole),MASK=Mix%vLPole(1:MaxPole))
  !
  if(Y>Eps) Mix%vXPole(1:MaxPole)=vX(1:MaxPole)/Y
  !
end subroutine MixPhase_NormCompo

real(dp) function MixPhase_Weit( &
& vSpc,      & !IN
& MM,        & !IN, mixing model
& Mix)         !IN, mixture phase, composition normalized
  use M_T_Species, only: T_Species
  !
  type(T_Species), intent(in),dimension(:)::vSpc
  type(T_MixModel),intent(in)   :: MM
  type(T_MixPhase),intent(in)   :: Mix !phase
  !
  real(dp):: Weit
  integer :: I
  !
  Weit=Zero
  do I=1,MM%NPole
    if(Mix%vLPole(I)) &
    & Weit= Weit &
    &     + Mix%vXPole(I) *vSpc(MM%vIPole(I))%WeitKg
  end do
  !
  MixPhase_Weit= Weit
  !
end function MixPhase_Weit

real(dp) function MixPhase_Volume( & !not complete,
!& TdgK,Pbar, & !IN
& vSpc,      & !IN, database
& SM,        & !IN, solution model
& Fas)         !IN, solution phase, composition normalized
  use M_T_Species, only: T_Species
  !
  !real(dp),        intent(in)   :: TdgK,Pbar
  type(T_Species), intent(in):: vSpc(:)
  type(T_MixModel),intent(in):: SM
  type(T_MixPhase),intent(in):: Fas !phase
  !
  real(dp),dimension(1:MaxPole):: vX !composition: mole fractions of the end members
  real(dp):: V, V_XS
  integer :: I
  !
  vX(1:SM%NPole)=Fas%vXPole(1:SM%NPole)
  V= Zero
  !
  do I=1,SM%NPole
    if(Fas%vLPole(I)) &
    & V= V + vX(I) *vSpc(SM%vIPole(I))%V0
  end do
  !
  !V= dot_product(vX(1:SM%NPole),vSpc(SM%vIPole(1:SM%NPole))%V0)
  !------------------------------ XS volume calc. not implemented !!! --
  V_XS=Zero
  !
  MixPhase_Volume= V + V_XS
  !
end function MixPhase_Volume

subroutine MixPhase_CalcMixing( & !
& TdgK,Pbar, & !IN
& MM,        & !IN
& Fas,       & !IN
& GMix,G_IdMix,G_XsMix)  !OUT
!--
!-- Output -
!--   G_IdMix: ideal part of free energy of mixing of mixture Phase at T,P
!--   GMix: free energy of mixing of phase Phase at T,P
!--   Gibbs free energy of the solution at T,P - GMeca + GMix
!--
  use M_Dtb_Const,only: R_jk
  !
  real(dp),        intent(in) :: Pbar, TdgK
  type(T_MixModel),intent(in) :: MM
  type(T_MixPhase),intent(in) :: Fas
  !real(dp),        intent(in) :: vLPole(:),vXPole(:)
  real(dp),        intent(out):: GMix,G_IdMix,G_XsMix
  !
  ! real(dp),allocatable:: vXAtom(:)
  !
  G_IdMix= MixModel_GibbsIdeal( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing model
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  !
  GMix= MixModel_GibbsMixRT( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing model
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  GMix= GMix *R_jk*TdgK
  !
  G_XsMix= GMix - G_IdMix
  !
  return
end subroutine MixPhase_CalcMixing

real(dp) function MixPhase_GibbsRT( & !
& TdgK,Pbar, & !IN
& vSpc,      & !IN, database
& MM,        & !IN, mixing model
& Fas)         !IN, mixture phase, composition normalized
!--
!-- not complete, check the SITE models !!!
!--
  use M_Dtb_Const,only: R_jk
  use M_T_Species,only: T_Species
  !
  real(dp),        intent(in):: TdgK,Pbar
  type(T_Species), intent(in):: vSpc(:)
  type(T_MixModel),intent(in):: MM
  type(T_MixPhase),intent(in):: Fas !phase
  !
  !~ real(dp),allocatable:: vXAtom(:)
  real(dp):: Gmix,Gmeca
  integer :: I
  !
  Gmix= MixModel_GibbsMixRT( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing MM
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  !
  !---------------------------------------------- "mechanical" mixing --
  Gmeca= Zero
  do I=1,MM%NPole
    if(Fas%vLPole(I)) &
    & Gmeca= Gmeca &
    &      + Fas%vXPole(I) *vSpc(MM%vIPole(I))%G0rt
  end do
  !----------------------------------------------/"mechanical" mixing --
  !
  !~ if(allocated(vXAtom)) deallocate(vXAtom)
  !
  MixPhase_GibbsRT= Gmix +Gmeca
  !
  return
end function MixPhase_GibbsRT

subroutine MixPhase_CalcActivs( & !
& TdgK,Pbar, & ! in
& MM,        & ! in:    mixing model
& F)           ! inout: mixture phase
!--
!-- calculate activities of end-members in phase F at given T,P
!--
  use M_Dtb_Const,only: R_jk
  !
  real(dp),        intent(in)   :: Pbar, TdgK
  type(T_MixModel),intent(in)   :: MM
  type(T_MixPhase),intent(inout):: F
  !
  logical :: Ok
  integer :: N
  character(len=80):: Msg
  !
  real(dp),allocatable:: vLGam(:),vLIdeal(:),vLnAct(:) !!,vLMarg(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)' ) "< MixPhase_CalcActivs"
  !
  N= size(F%vLnAct)
  allocate(vLGam(N),vLIdeal(N),vLnAct(N)) !!,vLMarg(N)
  !
  call MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MM,        & ! in: mixing model
  & F%vXPole,  & ! in
  & F%vLPole,  & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & ! out
  & vLIdeal,   & ! out
  & vLnAct)      ! out
  !
  F%vLIdeal(:)= vLIdeal(:)
  !! F%vLMarg(:)=  vLMarg(:)
  F%vLGam(:)=   vLGam(:)
  F%vLnAct(:)=  vLnAct(:)
  !
  deallocate(vLGam,vLIdeal,vLnAct) !!,vLMarg
  !
  if(iDebug>0) write(fTrc,'(A,/)' ) "</ MixPhase_CalcActivs"
  !
end subroutine MixPhase_CalcActivs

subroutine MixPhase_Constraint( &
& S,         & ! in: mixing model
& F,         & ! in: mixture phase
& Ok, Msg,   & ! out
& vC)          !
  type(T_MixModel),intent(in) :: S
  type(T_MixPhase),intent(in) :: F
  logical,         intent(out):: Ok
  character(*),    intent(out):: Msg
  real(dp),        intent(out):: vC(:)
  !
  integer:: I
  character:: c
  !
  select case(S%Model)

  case(Mix_Molecular,Mix_Felspar) ! ("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    vC(1)= One
    do I=1,S%NPole
      if(F%vLPole(I)) &
      & vC(1)= vC(1) - F%vXPole(I)
    end do
    !
  case(Mix_Site) ! ("SITE")
    vC(1)= One
    do I=1,S%NPole
      if(F%vLPole(I)) &
      & vC(1)= vC(1) - F%vXPole(I)
    end do
    !
  case default
    Ok= .false.
    c= char(ichar('0')+S%Model)
    Msg= c//"= invalid S%Model in MixPhase_ConstraintGrad"
    ! Msg= trim(S%Model)//"= invalid S%Model in MixPhase_Constraint"

  end select
  !
end subroutine MixPhase_Constraint

subroutine MixPhase_ConstraintGrad( & !
& S,         & ! in: mixing model
& F,         & ! in: mixture phase
& Ok, Msg,   & ! out
& vGC)         !
  type(T_MixModel),intent(in) :: S
  type(T_MixPhase),intent(in) :: F
  logical,         intent(out):: Ok
  character(*),    intent(out):: Msg
  real(dp),        intent(out):: vGC(:,:)
  !
  integer:: I
  character:: c
  !
  vGC= Zero
  !
  select case(S%Model)

  case(Mix_Molecular,Mix_Felspar) ! ("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    do I=1,S%NPole
      if(F%vLPole(I)) vGC(1,I)= -One
    end do
    !
  case(Mix_Site) ! ("SITE")
    do I=1,S%NPole
      if(F%vLPole(I)) vGC(1,I)= -One
    end do
    !
  case default
    Ok= .false.
    c= char(ichar('0')+S%Model)
    Msg= c//"= invalid S%Model in MixPhase_ConstraintGrad"

  end select
  !
end subroutine MixPhase_ConstraintGrad

end module M_T_MixPhase

! subroutine MixPhase_CalcActivities( & !
! & TdgK,Pbar, & ! in
! & MM,         & ! in: mixing model
! & F,         & ! in: mixture phase
! & Ok, Msg,   & ! out
! & vLGam,     & !
! !! & vLMarg,    & !
! & vLIdeal,   & !
! & vLnAct)      !
! !--
! !-------- calculate activities of end-members in phase F at given T,P --
! !--
!   use M_Dtb_Const,only: R_jk
!   use M_MixModel_Special
!   !
!   real(dp),        intent(in) :: Pbar, TdgK
!   type(T_MixModel),intent(in) :: MM
!   type(T_MixPhase),intent(in) :: F
!   logical,         intent(out):: Ok
!   character(*),    intent(out):: Msg
!   real(dp),        intent(out):: vLGam(:)
!   !! real(dp),        intent(out):: vLMarg(:)
!   real(dp),        intent(out):: vLIdeal(:)
!   real(dp),        intent(out):: vLnAct(:)
!   !
!   integer :: iP,iM
!   real(dp):: P
!   real(dp),allocatable:: vMonome(:)
!   !
!   P=Pbar !for future use ??
!   !
!   Ok= .true.
!   Msg= "Ok"
!   !
!   !F%vLPole(1:MM%NPole)= vX(1:MM%NPole)>Zero
!   !
!   vLGam(:)=   Zero !default
!   vLIdeal(:)= Zero !default
!   !
!   if(trim(MM%Model)=="SPECIAL") then
!     !
!     call MixModel_Special_Activities( &
!     & MM%Name,      &
!     & TdgK, Pbar,   &
!     & F%vXpole,     &
!     & F%vLPole,     &
!     & vLIdeal,      &
!     & vLGam         )
!     !
!   else
! 
!   select case(trim(MM%Model))
!     !
!     case("IDEAL","POLE","MOLECULAR","FELSPAR")
!       call MixModel_Pole_LnActivsIdeal(MM,F%vXpole,F%vLpole,vLIdeal)
! 
!     case("SITE")
!       call MixModel_Site_LnActivsIdeal(MM,F%vXatom,F%vLpole,Ok,vLIdeal)
! 
!     case default
!       Ok= .false.
!       Msg= trim(MM%Model)//"= invalid MM%Model in MixPhase_CalcActivities"
! 
!     end select
!     !
!     vLGam(1:MM%NPole)=Zero
!     !
!     !-------------------------- activ coeff related to Margules Terms --
!     if(MM%NMarg>0) then
!       allocate(vMonome(MM%NMarg))
!       !
!       select case(trim(MM%Model))
!       !
!       case("IDEAL","POLE","MOLECULAR","FELSPAR")
!         do iM=1,MM%NMarg
!           vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),F%vXPole)
!         end do
!         do iP=1,MM%NPole
!           if(F%vLPole(iP)) then
!             vLGam(iP)= MixModel_Pole_LnGammaMargules(MM,iP,vMonome,F%vXPole) /R_jk/TdgK
!           end if
!         end do
!       !
!       case("SITE")
!         do iM=1,MM%NMarg
!           vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),F%vXatom)
!         end do
!         do iP=1,MM%NPole
!           if(F%vLPole(iP)) then
!             vLGam(iP)= MixModel_Site_LnGammaMargules(MM,iP,vMonome,F%vXAtom) /R_jk/TdgK
!           end if
!         end do
!       !
!       case default
!         Ok= .false.
!         Msg= trim(MM%Model)//"= invalid MM%Model in MixPhase_CalcActivities"
!       !
!       end select
!       !
!       deallocate(vMonome)
!       !
!     end if
!     !-------------------------/ activ coeff related to Margules Terms --
!     !
!   end if
!   !
!   vLnAct(1:MM%NPole)= vLIdeal(1:MM%NPole) + vLGam(1:MM%NPole)
!   !
!   if(iDebug>0) then !------------------------------------------ trace --
!     write(fTrc,'(A)') "MixPhase_CalcActivs -> X,ActIdeal,Gamma,Activ"
!     write(fTrc,'(4A)') "Phase=",F%Name, "MixModel=",MM%Name
!     do iP=1,MM%NPole
!       if(F%vLPole(iP)) then
!         write(fTrc,'()')
!         write(fTrc,'(A,I2,A1,A15,A1,4(A4,G11.6,A1))') &
!         & "POLE",iP,           T_,&
!         & trim(MM%vNamPole(iP)),T_,&
!         & "Frc=",F%vXPole(iP), T_,&
!         & "XId=",exp(vLIdeal(iP)),T_,&
!         & "Gam=",exp(vLGam(iP)),  T_,&
!         & "Act=",exp(vLnAct(iP)), T_
!       end if
!     end do
!   end if !-----------------------------------------------------/ trace --
!   !
! end subroutine MixPhase_CalcActivities

! subroutine MixPhase_CalcMixing_( & !
! & TdgK,Pbar, & !IN
! & MM,        & !IN
! & Fas,       & !IN
! & GMix,G_IdMix,G_XsMix)  !OUT
! !--
! !-- Output -
! !--   G_IdMix: ideal part of free energy of mixing of mixture Phase for a given P,T
! !--   GMix: free energy of mixing of phase Phase for a given P,T
! !--   Gibbs free energy of the solution at P,T - GMeca + GMix
! !--
!   use M_Dtb_Const,only: R_jk
!   !
!   real(dp),        intent(in) :: Pbar, TdgK
!   type(T_MixModel),intent(in) :: MM
!   type(T_MixPhase),intent(in) :: Fas
!   !real(dp),        intent(in) :: vLPole(:),vXPole(:)
!   real(dp),        intent(out):: GMix,G_IdMix,G_XsMix
!   !
!   real(dp),allocatable:: vXAtom(:)
!   !
!   G_XsMix= Zero
!   G_XsMix= Zero
!   !
!   if(trim(MM%Model)=="SPECIAL") then
! 
!   else
!     !
!     select case(trim(MM%Model))
!     !
!     case("IDEAL","POLE","MOLECULAR","FELSPAR")
!       G_IdMix= -TdgK *MixModel_Pole_SConf(MM,Fas%vLPole,Fas%vXPole)
!       !
!       if(MM%NMarg>0) G_XsMix= &
!       & MixModel_Pole_XsMargules(MM, Fas%vLPole, Fas%vXpole)
!     !
!     case("SITE")
!       allocate(vXAtom(MM%NAtom))
!       call MixModel_XPoleToXSite(MM,Fas%vXPole,vXAtom)
!       !
!       G_IdMix= -TdgK *MixModel_Site_SConf(MM,vXAtom)
!       !
!       if(MM%NMarg>0) G_XsMix= &
!       & MixModel_Site_XsMargules(MM,vXAtom)
!     !
!     end select
!     !
!   end if
!   !
!   if(allocated(vXAtom)) deallocate(vXAtom)
!   !
!   GMix= G_IdMix +G_XsMix
!   !
!   return
! end subroutine MixPhase_CalcMixing_

