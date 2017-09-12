module M_Dtb_Read_DtbLogKTbl
  use M_Kinds
  use M_Trace
  use M_Dtb_Read_Vars,only: T_LisLogKTbl,LisLogKTbl,nLogKTbl
  !
  implicit none
  !
  private
  !
  public:: DtbLogKTbl_Read
  public:: DtbLogKTbl_Build
  public:: DtbLogK_TPCond_Read
  !
contains

subroutine DtbLogKTbl_Build
!--
!-- transfer linked list to array of T_DtbLogKTbl
!--
  use M_Dtb_Vars,only: vDtbLogKTbl
  !
  type(T_LisLogKTbl),pointer:: LisCur, LisPrev
  integer:: I
  !
  if(iDebug==5) print '(A)',"< DtbLogKTbl_Build"
  !
  if(allocated(vDtbLogKTbl)) deallocate(vDtbLogKTbl)
  allocate(vDtbLogKTbl(nLogKTbl))
  !
  I=0
  LisCur=> LisLogKTbl
  do while (associateD(LisCur))
    I=I+1
    vDtbLogKTbl(I)=LisCur%Value
    !
    if(iDebug==5) print '(2A)', vDtbLogKTbl(I)%Name,trim(vDtbLogKTbl(I)%Num)
    !
    LisPrev=>LisCur
    LisCur=> LisCur%next
    deallocate(LisPrev)
  end do
  !
  if(iDebug==5) print '(A)',"< DtbLogKTbl_Build"
  if(iDebug==5) call Pause_
  return
end subroutine DtbLogKTbl_Build

logical function LnkSpc_Found(L,W,M)
!--
!-- find index of string W in linked list
!-- return corresponding species Spc
!--
  use M_T_DtbLogKTbl, only: T_DtbLogKTbl
  !
  type(T_LisLogKTbl),pointer    :: L
  character(len=*),  intent(in) :: W
  type(T_DtbLogKTbl),intent(out):: M
  !
  type(T_LisLogKTbl),pointer:: P,pPrev
  integer::I
  !
  P=> L
  I=  0
  LnkSpc_Found=.false.
  do while (associateD(P))
    I=I+1
    if(trim(w)==trim(P%Value%Name)) then
      LnkSpc_Found=.true.
      M= P%Value
      exit;
    end if
    pPrev=> P
    P=>     P%next
  end do
end function LnkSpc_Found

subroutine DtbLogKTbl_Read(F,vEle,N)
!--
!-- build linked list from database in logK format
!--
  use M_Dtb_Const,only: T_CK
  use M_T_Element,only: T_Element,Formula_Read,Formula_Build,Element_Index
  use M_Files,    only: DirDtbLog,Files_Index_Write
  use M_IoTools,  only: DimV,getunit,lintowrd,appendtoend,str_upper
  use M_IoTools,  only: wrdtoreal,readrvalsv
  use M_Dtb_Read_Tools
  !
  use M_T_DtbLogKTbl, only: T_DtbLogKTbl,DimLogK_Max,DtbLogKTbl_Calc_Init
  use M_Dtb_Vars,     only: DtbLogK_Dim,DtbLogK_vTPCond
  !
  integer,        intent(in)   :: F !input file
  type(T_Element),intent(in)   :: vEle(:)
  integer,        intent(inout):: N
  !
  type(T_LisLogKTbl),pointer,save:: LisCur
  !
  character(len=512):: L,W
  type(T_DtbLogKTbl) :: M
  real(dp),dimension(dimV):: vX
  real(dp):: X,Rho
  logical :: EoL, fOk, EcformIsOk
  integer :: I,ff
  integer :: ios,ZSp,Div,mDum
  !
  integer,dimension(:),allocatable::vStoik !for Formula_Read
  !
  !--for formula translation
  character(len=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),dimension(:),allocatable:: vElement  !
  integer:: fFormula
  !--/
  !
  !--------------------------------------------var for header processing
  logical :: IsHeader
  integer,parameter:: nField= 13
  character(len=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  integer :: vifield(1:nField)
  !
  ! character(len=12):: vStrUnit(1:nField)
  !-------------------------------------------/var for header processing
  !
  if(idebug>1) write(fTrc,'(/,A)') "< DtbLogKTbl_Read"
  !
  !----------------------------------------------------------------trace
  if(iDebug>2) then
    call GetUnit(ff)
    open(ff,file="dtb_logktbl_check.tab")
  end if
  !---------------------------------------------------------------/trace
  !
  fFormula= 0
  CodFormula= "ECFORM"
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  !
  allocate(vStoik(1:size(vEle)))
  !
  !--initialize vStrField(:), vifield(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "ABBREV      ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
  &   "SIZE        ","VOLUME      ","DENSITY     " /)
  !
  L= "TYPE INDEX NAME ECFORM SIZE PARAMETERS" != default field list
  if(iDebug==4) print *,"< Default values"
  call FieldList_Read(L,vStrField,vifield)
  if(iDebug==4) print *,"</ Default values"
  !
  if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
  ! -> files contains compact formulas
    call GetUnit(fFormula)
    open(fFormula,file="debug_formula.log")
    write(fFormula,'(A,/)') "resuts of formula conversion"
  end if
  !---------------------------------/initialize vStrField(:), vifield(:)
  !
  !-----------------------------------build a linked list of all species
  !---------------------------------consistent with current element list
  DoFile: do
    !
    read(f,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile
    call AppendToEnd(L,W,EoL) ! if W=="END" append next word
    !
    !--------------------------------------------read header, if present
    !----------if the line begins with any member of array vStrField(:),
    !----------------------------------then it contains the line headers
    IsHeader= .false.
    do I=1,nField
      if(trim(W)==trim(vStrField(I))) then
        IsHeader= .true.
        exit
      end if
    end do
    !
    if(IsHeader) then
      L= trim(W)//" "//trim(L)
      if(iDebug>2) print *,"< values from file"
      call FieldList_Read(L,vStrField,vifield)
      if(iDebug>2) print *,"</ values from file"

      if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
        call GetUnit(fFormula)
        open(fFormula,file="debug_formula.log")
        write(fFormula,'(A,/)') "results of formula conversion"
      end if

      cycle DoFile
    end if
    !-------------------------------------------/read header, if present
    !
    !-----------------------------------------process first word of line
    select case(W)
    !
    case("ENDINPUT")
      exit DoFile
    !
    case("END","ENDSPECIES")
      exit DoFile
    !
    !----------------------------------------------------for old formats
    case("FORMULA")
      call LinToWrd(L,W,EoL)
      CodFormula= trim(W)

      if(CodFormula=="SCFORM") then
        if(vifield(4)==0) then
          vifield(4)= vifield(5)
          vifield(5)= 0
        end if
        if(iDebug>2) then
          call GetUnit(fFormula)
          open(fFormula,file="debug_formula.log")
          write(fFormula,'(A,/)') "resuts of formula conversion"
        end if
      end if

      cycle DoFile !-----------------------------------------------cycle
    !
    !---------------------------------------------/ for old formats --
    !old! case default
    !old!   !--- ignore line beginning with unknown keyword !!! ??? --
    !old!   cycle DoFile !-----------------------------------------cycle
    !old! !
    case default
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= trim(W)//" "//trim(L)
    !
    end select
    !----------------------------------------/process first word of line
    !
    M%Num= "0"
    !---------------------------------------------------scan a data line
    do I= 1,vifield(10)-1
    ! read words up to position before the parameters
      !
      call LinToWrd(L,W,EoL,"NO")
      !
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SIZE        ","VOLUME      ","DENSITY     ","FITTING     ","PARAMETERS  ", &
      !&   "ABBREV      ","SOURCE      "/)
      !
      if(I==vifield(1)) then  ! type
        !
        ! call Str_Upper(W)  ;  M%Typ= trim(W)
        call Str_Upper(W)
        select case(trim(W))
        case("AQU")  ;  M%Typ="AQU"
        case("MIN")  ;  M%Typ="MIN"
        case("GAS")  ;  M%Typ="GAS"
        case("LIQ")  ;  M%Typ="LIQ"
        case("XCH")  ;  M%Typ="XCH"
        case default ;  call Stop_(trim(W)//"<< unknown type in database !!...")          
        end select
        !
      end if

      if(I==vifield(2)) then  ! INDEX
        call Str_Upper(W)  ;  M%Num= trim(W)
      end if

      if(I==vifield(3)) then  ! NAME
        call Str_Upper(W)  ;  M%Name= trim(W)
      end if

      if(I==vifield(4)) then  ! SCFORM
        !
        call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        if(.not.EcformIsOk) cycle DoFile !-------------------------cycle
        !
        call Str_Upper(W)  ;  M%Formula=trim(W)
      end if

      if(I==vifield(5)) then  ! ECFORM
        call Str_Upper(W)  ;  M%Formula=  trim(W)
      end if

      !if(I==vifield(9)) then
      !  call Str_Upper(W)  ;  CodFitting= trim(W)
      !end if

      if(I==vifield(11)) call WrdToReal(W,X) ! size
      if(I==vifield(12)) call WrdToReal(W,X) ! VOLUME
      if(I==vifield(13)) call WrdToReal(W,X) ! DENSITY

    end do
    !
    call ReadRValsV(L,mDum,vX)
    !--------------------------------------------------/scan a data line
    !
    if(mDum<DtbLogK_Dim) &
    & call Stop_("Dimension of LogK array < Dimension of TP table !!")
    !
    !-------------------------------------------------------read formula
    call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    if(.not. fOk) cycle DoFile !-----------------------------------cycle
    !------------------------------------------------------/read formula
    !
    !-- compute Molecular Weit
    M%WeitKg= dot_product(vStoik(:),vEle(:)%WeitKg) /real(Div)
    M%Div= Div
    M%Chg= Zsp
    !
    !---------size parameter for aqu'species, or density for min'species
    if(    vifield(11)/=0 & ! size
    & .or. vifield(12)/=0 & ! VOLUME, m^3 /mole ?
    & .or. vifield(13)/=0 & ! DENSITY, Kg /m^3
    & ) then
      !! print *,"M%Name=",M%Name
      select case(M%Typ)
      case("AQU")
        M%AquSize=X
      case("MIN","GAS")
        if(vifield(12)/=0) then
          M%V0R= X *1.0D-6 ! for THERMODDEM database, CM^3 to M^3
          Rho= M%WeitKg /M%V0R
        else
          Rho=X
          M%V0R= M%WeitKg /Rho
        end if
      end select
    end if
    !--default value of size parameter for charged aqu'species
    if(M%Typ=="AQU" .and. M%Chg/=0 .and. M%AquSize<=Zero) M%AquSize= 3.72D0
    !----/read size parameter for aqu'species, or volume for min'species
    !
    !write(fff,'(A,A1)') trim(L) !-> write the logK's
    !
    !--- read list of logK's --
    M%DimLogK=  DtbLogK_Dim
    M%vLogK(1:DimLogK_Max)= vX(1:DimLogK_Max)
    !---/ read list of logK's --
    !
    M%vTdgK(1:DtbLogK_Dim)= DtbLogK_vTPCond(1:DtbLogK_Dim)%TdgC + T_CK
    !
    call DtbLogKTbl_Calc_Init(M)
    !
    !! if(iDebug==5) write(ff,'(A3,A1,A15,A1,A39,A1,I3)') &
    !! & M%Typ,T_,M%Name,T_,M%Formula,T_,N
    !
    if(iDebug>2) call Test(ff,M) !--trace
    !
    !---------------------------------------------"fill" the linked list
    N=N+1
    if(N==1) then
      allocate(LisLogKTbl)
      nullify(LisLogKTbl%next)
      LisLogKTbl%Value= M
      LisCur=> LisLogKTbl
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    end if
    !--------------------------------------------/"fill" the linked list
    !
  end do DoFile
  !
  deallocate(vStoik)
  !
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(iDebug>2) close(ff)!--trace
  !
  if(idebug>1) write(fTrc,'(A,/)') "< DtbLogKTbl_Read"
  if(N==0) call Stop_("NO SPECIES FOUND ....")
  !
end subroutine DtbLogKTbl_Read

subroutine Test(ff,M)
  use M_T_Species
  use M_Numeric_Const,only: Ln10
  use M_T_DtbLogKTbl
  
  integer,intent(in):: ff
  type(T_DtbLogKtbl),intent(in) :: M
  
  type(T_Species):: S
  real(dp):: vT(1:5),vLogK(1:5),Pbar
  integer :: I
  
  vT(1:5)= (/0.01,25.,50.,75.,100./)
  Pbar= 1.0D0
  do I=1,5
    call DtbLogKtbl_Calc(M,vT(I)+273.15D0,Pbar,S)
    vLogK(I)= - S%G0RT /Ln10
  end do
  
  write(ff,'(2(A,1X))',advance="NO") M%Typ,M%Name
  write(ff,'(5(G15.6,1X))',advance="NO") (vLogK(I),I=1,5)
  if(M%V0R>Zero) then
    write(ff,'(G15.6,1X)',   advance="NO") M%WeitKg /M%V0R
  else
    write(ff,'(A,1X)',   advance="NO") "_"
  end if
  write(ff,'(A,1X)') M%Formula
  
end subroutine Test

subroutine DtbLogK_TPCond_Read(N)
!--
!-- reads block TP.TABLE, build DtbLogK_vTPCond
!--
  use M_IoTools !,only:GetUnit,dimV
  use M_Files,     only: NamFInn
  use M_Dtb_Const, only: T_CK
  use M_T_Tpcond,  only: T_TPCond
  use M_Fluid_Calc,only: Eos_H2O_psat
  !
  use M_Dtb_Vars,    only: DtbLogK_vTPCond,DtbLogK_Dim
  use M_T_DtbLogKTbl,only: DimLogK_Max
  !
  integer,intent(out):: N
  !
  logical :: EoL,Ok
  integer :: i,mDum,ios,f
  character(len=512):: L,W,W1
  !
  real(dp):: TdgK,Pbar,PSatBar
  real(dp):: vX(dimV)
  type(T_TPCond):: vCond(dimV)
  !
  if(idebug>1) write(fTrc,'(/,A)') "< DtbLogK_TPCond_Read"
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  N= 0
  Ok=.false.
  !
  vCond(:)%TdgC= Zero
  vCond(:)%Pbar= Zero
  vCond(:)%Name= "NONE"
  !
  DoFile: do
    !
    read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!')   cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit  DoFile
    !
    case("TP.TABLE","TPTABLE")
      Ok=  .true.
      N= dimV
      !
      !---------------------------------------------- read table loop --
      DoTPtable: do
        !
        read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,Eol)
        if(W(1:1)=='!') cycle DoTPtable
        call AppendToEnd(L,W,EoL)
        !
        select case(trim(W))
        !
        case("ENDINPUT"); exit DoFile
        case("END","ENDTPTABLE","ENDTP.TABLE"); exit DoTPtable
        !
        case("TDGC","TDGK", "PBAR","PMPA","TDEGC","TDEGK")
          !
          !!call ReadRValsV(L,mDum,vX)
          !
          i=0
          do
            call LinToWrd(L,W1,EoL)
            !print *,W1
            i=i+1
            if(i>DimV) exit
            if(trim(W1)=="PSAT") then
              if(vCond(i)%TdgC>Zero) then
                call Eos_H2O_psat(vCond(i)%TdgC+T_CK,vX(i))
                !print *,"TdgC= ", vCond(i)%TdgC
                !print *,"psat= ", vX(i)  ;  call Pause_
              else
                call Stop_("error in TP.TABLE: invalid TdgC for Psat computation")
              end if
            else
              call WrdToReal(W1,vX(i))
            end if
            if(EoL) exit
          end do
          !
          mDum= i
          N=min(N,mDum)
          !
          if(idebug>1) write(fTrc,'(A,A1,A,2I3)') trim(W),T_," DIM= ", mDum, N
          !
          select case(trim(W))
          case("TDGC");  vCond(:)%TdgC=vX(:)
          case("TDGK");  vCond(:)%TdgC=vX(:)-T_CK
          !
          case("PBAR");  vCond(:)%Pbar=vX(:)       !pressure in Bar
          case("PMPA");  vCond(:)%Pbar=vX(:)*1.D-1 !MegaPascal to Bar
          !
          case("TDEGC"); vCond(:)%TdgC=vX(:)       ; call Warning_("Depreciated Keyword TDEGC")
          case("TDEGK"); vCond(:)%TdgC=vX(:)-T_CK  ; call Warning_("Depreciated Keyword TDEGK")
          !
          end select
        !endcase("TDGC",..
        !
        case("NAME")
          i=0
          do
            call LinToWrd(L,W,EoL)
            i=i+1
            if(I>DimV) exit
            vCond(i)%Name=trim(W)
            if(idebug>1) write(fTrc,'(I3,1X,A15)') i,vCond(i)%Name
            if(EoL) exit
          end do
        !endcase("NAME")
        !
        case default; call Stop_(trim(W)//"<<unknown Keyword...")
        !
        end select !end if
        !end if
      end do DoTPtable
      !---------------------------------------------/ read table loop --
    !end case("ENDTP.TABLE")
    end select
    !
  end do DoFile
  !
  close(f)
  !
  if(Ok) then
    !! if(any(vCond(1:N)%TdgC<Zero)) call Stop_("all temperatures should be >0 !!!")
    if(any(vCond(1:N)%Pbar<1.D-9)) call Stop_("all pressures should be >0 !!!")
    !
    do I=1,N
      TdgK= vCond(I)%TdgC +T_CK
      Pbar= vCond(I)%Pbar
      if(TdgK <Zero) TdgK= 200.0D0 ! Zero
      !
      !-- check that pressure is not below
      !-- the vapor saturation curve for pure water
      !-- critical point H2O- 647.25D0 TdgK / 220.55D0 Pbar
      if (TdgK<=647.25D0) then; call Eos_H2O_psat(TdgK,PSatBar)
      else                    ; PSatBar= 220.55D0
      end if
      !-- if pressure is too low, adjust to saturation pressure at TdgK
      if(Pbar<PSatBar) Pbar=PSatBar
      !
      vCond(I)%TdgC= TdgK -T_CK
      vCond(I)%Pbar= Pbar
    end do
    !
    if(N>DimLogK_Max) N= DimLogK_Max
    !
    DtbLogK_Dim= N
    !
    if(allocated(DtbLogK_vTPCond)) deallocate(DtbLogK_vTPCond)
    allocate(DtbLogK_vTPCond(N))
    DtbLogK_vTPCond(1:N)= vCond(1:N)
    !
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ DtbLogK_TPCond_Read"
  !
end subroutine DtbLogK_TPCond_Read

end module M_Dtb_Read_DtbLogKTbl
