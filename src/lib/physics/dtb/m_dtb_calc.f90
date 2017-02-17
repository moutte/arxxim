module M_Dtb_Calc
!--
!-- routines to compute species properties from databases
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,Stop_,T_,Pause_
  implicit none
  !
  private
  !
  public:: Species_TP_Update_fromDtb
  public:: SpeciesMin_TP_Update_fromDtb
  public:: Dtb_TP_Check
  public:: DtbSpc_Table_Build
  public:: SpeciesDtb_ToSpecies
  !
contains

subroutine SpeciesDtb_ToSpecies( &
& vSpcDtb, & ! IN
& vSpc)      ! OUT

  use M_T_Species,only: T_Species,T_SpeciesDtb,Species_Index
  !
  use M_Dtb_Vars, only: &
  & vDtbMinHkf,vDtbMinThr,vDtbAquHkf,vDtbLogKtbl,vDtblogKanl
  !
  type(T_SpeciesDtb),intent(in)   :: vSpcDtb(:)
  type(T_Species),   intent(inout):: vSpc(:)
  !
  type(T_Species):: S
  integer:: I,J,N
  !
  N= size(vSpcDtb)
  !
  do J= 1,N

    I= vSpcDtb(J)%Indx

    select case(vSpcDtb(J)%DtbModel)

    case("H2O_HGK")  ! HGK= Haar-Gallagher-Kell
      S%NamSp=    "H2O"
      S%WeitKg=   0.0180152D0
      S%Formula=  "H(2)O(1)"
      S%Typ=      "AQU"
      S%AquSize=  0.D0

    case("AQU_HKF")
      S%NamSp=    trim(vDtbAquHkf(I)%Name)
      S%Formula=  trim(vDtbAquHkf(I)%Formula)
      S%WeitKg=   vDtbAquHkf(I)%WeitKg
      !S%AquSize=  vDtbAquHkf(I)%AquSize
      S%Typ=      "AQU"

    !! case("AQU_THR")

    case("MIN_HKF","GAS_HKF")
      S%NamSp=    trim(vDtbMinHkf(I)%Name)
      S%Formula=  trim(vDtbMinHkf(I)%Formula)
      S%WeitKg=   vDtbMinHkf(I)%WeitKg
      S%Typ=      vDtbMinHkf(I)%Typ

    case("MIN_THR","GAS_THR")
      S%NamSp=    trim(vDtbMinThr(I)%Name)
      S%Formula=  trim(vDtbMinThr(I)%Formula)
      S%WeitKg=   vDtbMinThr(I)%WeitKg
      S%Typ=      vDtbMinThr(I)%Typ

    case("LOGKTBL")
      S%NamSp=    trim(vDtbLogKtbl(I)%Name)
      S%Formula=  trim(vDtbLogKtbl(I)%Formula)
      S%WeitKg=   vDtbLogKtbl(I)%WeitKg
      S%AquSize=  vDtbLogKtbl(I)%AquSize
      S%Typ=      vDtbLogKtbl(I)%Typ

    case("LOGKANL")
      S%NamSp=    trim(vDtbLogKanl(I)%Name)
      S%Formula=  trim(vDtblogKanl(I)%Formula)
      S%WeitKg=   vDtbLogKanl(I)%WeitKg
      S%AquSize=  vDtbLogKanl(I)%AquSize
      S%Typ=      vDtbLogKanl(I)%Typ

    end select

    S%iDiscret= 0
    S%iDtb=     J
    vSpc(J)=    S

    ! print *,S%NamSp
    
  end do

  ! pause_
  
  !to prevent some bugs,
  !put solvent species H2O as species 1 in vSpc
  N= Species_Index("H2O",vSpc)
  if(N/=0 .and. N/=1) then
    S=       vSpc(1)
    vSpc(1)= vSpc(N)
    vSpc(N)= S
  end if

  if(iDebug<4) return

  print *,"SpeciesDtb_ToSpecies"
  do i=1,N
    print *,vSpc(i)%NamSp,vSpc(i)%Typ,"=",vSpcDtb(vSpc(i)%iDtb)%DtbModel
  end do

  return
end subroutine SpeciesDtb_ToSpecies

subroutine Species_TP_Update_fromDtb( &
& TdgK,Pbar,PropsH2O,vSpcDtb, &
& S)
!--
!-- compute thermodyn.prop's of species S at TdgK,Pbar,PropsH2O
!-- according to S%Model and S%iDtb
!-- should be called only when S%iDtb/=0
!--
  use M_T_Species,  only: T_Species,T_SpeciesDtb
  !
  use M_T_DtbH2OHkf, only: T_DtbH2OHkf
  use M_T_DtbAquHkf, only: DtbAquHkf_Calc,DtbAquHkf_CalcThr
  use M_T_DtbMinHkf, only: DtbMinHkf_Calc
  use M_T_DtbMinThr, only: DtbMinThr_Calc
  use M_T_DtbLogKtbl,only: DtbLogKtbl_Calc
  use M_T_DtbLogKanl,only: DtbLogKanl_Calc
  !
  use M_Dtb_Vars,   only: vDtbAquHkf,vDtbMinHkf,vDtbMinThr  ! data-
  use M_Dtb_Vars,   only: vDtbLogKtbl,vDtbLogKanl           ! bases
  !
  use M_Fluid_Calc, only: EosFluid_Calc
  !!use M_Fluid_Calc, only: CalcGH2O_Supcrt !,Eos_H2O_Haar
  
  !---------------------------------------------------------------------
  real(dp),          intent(in) :: TdgK,Pbar
  type(T_DtbH2OHkf), intent(in) :: PropsH2O
  type(T_SpeciesDtb),intent(in) :: vSpcDtb(:)
  !
  type(T_Species),   intent(inout):: S
  !---------------------------------------------------------------------
  integer :: i,j
  real(dp):: LnFug
  logical :: Ok
  !---------------------------------------------------------------------
  
  i= S%iDtb
  if(i==0) return
  !
  j= vSpcDtb(i)%Indx
  !
  select case(trim(vSpcDtb(i)%DtbModel))

  case("H2O_HGK")
    call EosFluid_Calc("H2O","HAAR.GHIORSO",TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0,LnFug,Ok)
  ! case("H2O_HGK") ; call CalcGH2O_Supcrt(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
  !                                    !IN                         OUT
  case("AQU_HKF") ; call DtbAquHkf_Calc(vDtbAquHkf(j),PropsH2O,    S)
  case("AQU_THR") ; call DtbAquHkf_CalcThr(vDtbAquHkf(j),PropsH2O, S)

  case("MIN_HKF") ; call DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,   S)
  case("GAS_HKF") ; call DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,   S)

  case("MIN_THR") ; call DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,   S)
  case("GAS_THR") ; call DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,   S)

  case("LOGKTBL") ; call DtbLogKtbl_Calc(vDtbLogKtbl(j),TdgK,Pbar, S)
  case("LOGKANL") ; call DtbLogKanl_Calc(vDtbLogKanl(j),TdgK,Pbar, S)

  end select
  
  ! if(trim(S%NamSp)=="WATER") &
  ! & write(94,'(A,3(G15.6,1X))') "T,P,V_H2O_cm3=",TdgK,Pbar,S%V0*1.0D6

  !if(trim(S%Typ)=="MIN") then
  !  print *,"S%Name,S%V0= ",S%NamSp,S%V0
  !  call pause_
  !end if
 
  return
end subroutine Species_TP_Update_fromDtb

subroutine SpeciesMin_TP_Update_fromDtb(TdgK,Pbar,vSpcDtb, S)
!--
!-- calculate thermodyn.properties of species S at TdgK,Pbar
!-- same as DtbSpc_TP_Update, for non'aqu'species only -> no need for PropsH2O
!-- should be called only when S%iDtb/-0
!--
  use M_T_Species,  only: T_Species,T_SpeciesDtb
  !
  use M_T_DtbH2OHkf, only: T_DtbH2OHkf
  use M_T_DtbMinHkf, only: DtbMinHkf_Calc
  use M_T_DtbMinThr, only: DtbMinThr_Calc
  use M_T_DtbLogKtbl,only: DtbLogKtbl_Calc
  use M_T_DtbLogKanl,only: DtbLogKanl_Calc
  !
  use M_Dtb_Vars,   only: &
  & vDtbAquHkf,vDtbMinHkf,vDtbMinThr,vDtbLogKtbl, vDtbLogKanl !-> the databases
  !
  use M_Fluid_Calc, only: EosFluid_Calc
  !!use M_Fluid_Calc, only: CalcGH2O_Supcrt !,Eos_H2O_Haar
  !---------------------------------------------------------------------
  real(dp),          intent(in) :: TdgK,Pbar
  type(T_SpeciesDtb),intent(in) :: vSpcDtb(:)
  !
  type(T_Species),   intent(inout):: S
  !---------------------------------------------------------------------
  real(dp):: LnFug
  integer :: i,j
  logical :: Ok
  !---------------------------------------------------------------------
  i= S%iDtb
  !
  if(i==0) return
  !
  j= vSpcDtb(i)%Indx
  !
  select case(vSpcDtb(i)%DtbModel)
    case("H2O_HGK")
      call EosFluid_Calc( &
      & "H2O","HAAR.GHIORSO",TdgK,Pbar, &
      & S%G0rt,S%H0,S%S0,S%V0,LnFug,Ok)    
      ! case("H2O_HGK"); call CalcGH2O_Supcrt(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
    case("MIN_HKF")  ;  call DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)
    case("MIN_THR")  ;  call DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)
    case("GAS_HKF")  ;  call DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)
    case("GAS_THR")  ;  call DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)
    case("LOGKTBL")  ;  call DtbLogKtbl_Calc(vDtbLogKtbl(j),TdgK,Pbar,S)
    case("LOGKANL")  ;  call DtbLogKanl_Calc(vDtbLogKanl(j),TdgK,Pbar,S)
  end select
  !
  return
end subroutine SpeciesMin_TP_Update_fromDtb

subroutine Dtb_TP_Check( &
& DtbFormat,DtbLogK_vTPCond,Psat_Auto, &
& TdgK,Pbar,Ok,Msg)
!--
!-- check that TdgK is in the validity range of the database,
!-- and (in case of logK database) compute corresponding pressure
!--
  use M_T_Tpcond,  only: T_TPCond
  use M_Dtb_Const, only: T_CK, PminHSV, PmaxHSV, TCminHSV, TCmaxHSV, Pref
  use M_Fluid_Calc,only: Eos_H2O_Psat
  !
  character(len=*),intent(in)   :: DtbFormat
  type(T_TPCond),  intent(in)   :: DtbLogK_vTPCond(:)
  logical,         intent(in)   :: Psat_Auto
  real(dp),        intent(in)   :: TdgK
  real(dp),        intent(inout):: Pbar
  logical,         intent(out)  :: Ok
  character(len=*),intent(out)  :: Msg
  !
  real(dp):: Tmin,Tmax,Psat !,Pmin,Pmax
  !
  Ok=  .true.
  Msg= "OK"
  !
  !----------------------------------------------------check Temperature
  select case(DtbFormat)

  case("LOGKTBL")
    Tmin= DtbLogK_vTPCond(1)%TdgC +T_CK
    Tmax= DtbLogK_vTPCond(size(DtbLogK_vTPCond))%TdgC +T_CK

  case default
    Tmin= TCminHSV + T_CK
    Tmax= TCmaxHSV + T_CK

  end select
  !
  if(TdgK <Tmin .or. TdgK>Tmax) then
    Ok= .false.
    if(iDebug>2) print *,"TdgK,Tmin,Tmax",TdgK,Tmin,Tmax
    Msg= "Temperature outside validity range of the database"
    return !-------------------------------------------------------error
  end if
  !--------------------------------------------------/ check Temperature 
  !
  !------------------------------------------------------ check Pressure 
  select case(DtbFormat)
  !-------------------------------------------------------logK databases
  case("LOGKTBL","LOGKANL")
    if(TdgK <= 100.0D0+T_CK) then
      if(Pbar<Pref) then
        if(Psat_Auto) then
          Pbar= Pref
        else
          Ok= .false.
          Msg= "Pressure < Psat(T) = not valid for HKF model"
          return !-------------------------------------------------error
        end if
      end if
    else
      call Eos_H2O_Psat(TdgK,Psat)
      if(Pbar<Psat) then
        if(Psat_Auto) then
          Pbar= Psat
        else
          Ok= .false.
          Msg= "Pressure < Psat(T) = not valid for HKF model"
          return
        end if
      end if
    end if
  !------------------------------------------------------/logK databases
  !------------------------------------------------------other databases
  case default
    if(Pbar <PminHSV .or. Pbar>PmaxHSV) then
      Ok= .false.
      Msg= "Pressure outside validity range of the database"
      return !-----------------------------------------------------error
    end if
    !
    !---------------------------------------- if P-PsatH2O(T) then error
    call Eos_H2O_Psat(TdgK,Psat)
    !! Pbar= max(Pbar,Psat)
    if(Pbar < Psat) then
      if(Psat_Auto) then
        Pbar= Psat
      else
        Ok= .false.
        Msg= "Pressure < Psat(T) = not valid for HKF model"
        return !---------------------------------------------------error
      end if
    end if
  !-----------------------------------------------------/other databases
  end select
  !------------------------------------------------------/check Pressure
  !
end subroutine Dtb_TP_Check

subroutine DtbSpc_Table_Build( &
& vTPCond,vSpcDtb,vSpc, &
& tGrt,tVol)
  use M_T_Species,only: T_Species,T_SpeciesDtb
  use M_T_Tpcond, only: T_TPCond
  use M_Dtb_Const,only: T_CK
  use M_T_DtbH2OHkf,only: DtbH2OHkf_Calc,T_DtbH2OHkf
  
  type(T_TPCond),    intent(in) :: vTPCond(:)
  type(T_SpeciesDtb),intent(in) :: vSpcDtb(:)
  type(T_Species),   intent(in) :: vSpc(:)
  real(dp),          intent(out):: tGrt(:,:)
  real(dp),          intent(out):: tVol(:,:)
  
  integer :: iTP,jSp
  real(dp):: TdgK,Pbar
  type(T_DtbH2OHkf)::PropsH2O
  type(T_Species):: S
  
  do iTP=1,size(vTPCond)
    !
    TdgK= vTPcond(iTP)%TdgC + T_CK
    Pbar= vTPcond(iTP)%Pbar
    !--- solvent properties, for aqu'species
    call DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
    !---
    !
    do jSp=1,size(vSpc)
      !
      S= vSpc(jSp)
      !
      if(S%iDtb>0) then
        call Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,S)
        !-> update vSpc(jSp)%G0rt !-> =G/RT
        tGrt(jSp,iTP)= S%G0rt
        tVol(jSp,iTP)= S%V0
      !elseif(vSpc(jSp)%iDtb>0) then
      end if
      !
    end do
    !
  end do
  
end subroutine DtbSpc_Table_Build

end module M_Dtb_Calc

