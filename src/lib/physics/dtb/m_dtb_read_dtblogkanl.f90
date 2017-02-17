module M_Dtb_Read_DtbLogKAnl
  use M_Kinds
  use M_Trace
  use M_Dtb_Read_Vars,only: T_LisLogKAnl,LisLogKAnl,nLogKAnl
  !
  implicit none
  !
  private
  !
  public:: DtbLogKAnl_Read
  public:: DtbLogKAnl_Build
  !
contains

subroutine DtbLogKAnl_Build
!--
!-- transfer linked list to array of T_DtbLogKAnl
!--
  use M_Dtb_Vars,only: vDtbLogKAnl
  !
  type(T_LisLogKAnl),pointer:: LisCur, LisPrev
  integer:: I
  !
  if(iDebug==5) print '(A)',"< DtbLogKAnl_Build"
  !
  if(allocated(vDtbLogKAnl)) deallocate(vDtbLogKAnl)
  allocate(vDtbLogKAnl(nLogKAnl))
  !
  I=0
  LisCur=> LisLogKAnl
  do while (associateD(LisCur))
    I=I+1
    vDtbLogKAnl(I)= LisCur%Value
    if(iDebug==5) &
    & print '(2A)', vDtbLogKAnl(I)%Name,trim(vDtbLogKAnl(I)%Num)
    LisPrev=>LisCur
    LisCur=> LisCur%next
    deallocate(LisPrev)
  end do
  !
  if(iDebug==5) print '(A)',"</ DtbLogKAnl_Build"
  if(iDebug==5) call Pause_
  return
end subroutine DtbLogKAnl_Build

logical function LnkSpc_Found(L,W) !,M)
!--
!-- find index of string W in linked list
!-- return corresponding species Spc
!--
  use M_T_DtbLogKAnl, only: T_DtbLogKAnl
  !
  type(T_LisLogKAnl),pointer    :: L
  character(len=*),  intent(in) :: W
  ! type(T_DtbLogKAnl),intent(out):: M
  !
  type(T_LisLogKAnl),pointer:: P,pPrev
  integer::I
  !
  P=> L
  I=  0
  LnkSpc_Found=.false.
  do while (associateD(P))
    I=I+1
    if(trim(W)==trim(P%Value%Name)) then
      LnkSpc_Found=.true.
      !M= P%Value
      exit
    end if
    pPrev=> P
    P=>     P%next
  end do
  !
end function LnkSpc_Found

subroutine DtbLogKAnl_Read(F,vEle,N)
!--
!-- build linked list from database in logK/analytic format
!--
  use M_Dtb_Const,only: T_CK
  use M_T_Element,only: T_Element,Formula_Read,Formula_Build,Element_Index
  use M_Files,    only: DirDtbLog,Files_Index_Write
  use M_IoTools
  use M_Dtb_Read_Tools
  !
  use M_T_DtbLogKAnl, only: T_DtbLogKAnl,DtbLogKAnl_New
  use M_Dtb_Vars,     only: DtbLogK_Dim,DtbLogK_vTPCond
  !
  integer,        intent(in)   :: F !input file
  type(T_Element),intent(in)   :: vEle(:)
  integer,        intent(inout):: N
  !
  type(T_LisLogKAnl),pointer,save:: LisCur
  !
  character(len=512):: L
  character(len=80):: W
  type(T_DtbLogKAnl) :: M
  logical :: EoL, fOk
  logical :: EcformIsOk
  integer :: ios,ZSp,Div,mDum,ff,I
  real(dp):: X,Rho
  real(dp):: vX(dimV)
  !
  !! character(len=80):: vWord(1:12)
  !
  !--- for header processing --
  integer,parameter:: nField= 13
  character(len=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  integer:: vifield(1:nField)
  logical:: IsHeader
  !
  ! character(len=12):: vStrUnit(1:nField)
  !---/ for header processing --
  !
  !~ character(len=512):: FieldList_Default= &
  !~ & "TYPE SOURCE NAME SCFORM size parameterS"
  !
  integer,dimension(:),allocatable::vStoik !for Formula_Read
  !
  !--- for formula translation --
  character(len=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),dimension(:),allocatable:: vElement  !
  integer:: fFormula
  !---/
  !
  character(len=15):: CodFitting ! type of the fitting function, default: PHREEQ
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbLogKAnl_Read"
  !
  !--------------------------------------------------------------- trace
  if(iDebug>2) then
    call GetUnit(ff)
    open(ff,file="dtb_logKAnl_check.tab")
  end if
  !------------------------------------------------------------- / trace
  !
  fFormula= 0
  CodFormula= "ECFORM"
  !
  !------------------------------ initialize vStrField(:), vifield(:) --
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
  &   "SIZE        ","VOLUME      ","DENSITY     " /)
  !
  L= "TYPE INDEX NAME SCFORM size parameterS" != default field list
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
  CodFitting= "PHREEQC"
  !
  !--for formula translation
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  !--/for formula translation
  !
  allocate(vStoik(1:size(vEle)))
  !
  !------------------------------- build a linked list of all species --
  !----------------------------- consistent with current element list --
  DoFile: do

    read(f,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile
    call AppendToEnd(L,W,EoL) ! if W=="END" append next word
    !
    call DtbLogKAnl_New(M)
    !
    IsHeader= .false.
    !--------------------------------------------read header, if present 
    !--------- if the line begins with any member of array vStrField(:), 
    !--------------------------------- then it contains the line headers 
    do I=1,nField
      if(trim(W)==trim(vStrField(I))) then
        IsHeader= .true.
        exit
      end if
    end do
    if(IsHeader) then
      L= trim(W)//" "//trim(L)
      if(iDebug>2) print *,"< values from file"
      call FieldList_Read(L,vStrField,vifield)
      if(iDebug>2) print *,"</ values from file"

      if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
        call GetUnit(fFormula)
        open(fFormula,file="debug_formula.log")
        write(fFormula,'(A,/)') "resuts of formula conversion"
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
    !---------------------------------------------- for old formats --
    case("FITTING")
      call LinToWrd(L,W,EoL)
      CodFitting= trim(W)
      cycle DoFile !---------------------------------------------cycle
    !
    case("FORMULA")
      call LinToWrd(L,W,EoL)
      CodFormula= trim(W)
      !
      if(CodFormula=="SCFORM") then
        if(vifield(4)==0) then
          vifield(4)= vifield(5)
          vifield(5)= 0
        end if
        if(iDebug>2 .and. fFormula==0) then
          call GetUnit(fFormula)
          open(fFormula,file="debug_formula.log")
          write(fFormula,'(A,/)') "resuts of formula conversion"
        end if
      end if
      !
      cycle DoFile !---------------------------------------------cycle
    !
    !old! !---------------------------------------------/ for old formats
    !old! case default
    !old!   !--- ignore line beginning with unknown keyword !!! ??? --
    !old!   cycle DoFile !----------------------------------------- cycle
    !old! !
    !
    case default
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= trim(W)//" "//trim(L)
    !
    end select
    !----------------------------------------/process first word of line
    !
    M%Num=     "0"
    M%Name=    "0"
    M%Formula= "0"
    !---------------------------------------------------scan a data line
    do I= 1,vifield(10)-1
    ! read words up to position before the parameters
      !
      call LinToWrd(L,W,EoL,"NO")
      !
      !vStrField(1:12)= &
      !!!!!"____________","____________","____________","____________","____________",
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
      !&   "SIZE        ","VOLUME      ","DENSITY     "/)
      !
      if(I==vifield(1)) then  !type
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
      end if

      if(I==vifield(2)) then  !INDEX
        call Str_Upper(W)  ;  M%Num= trim(W)
      end if

      if(I==vifield(3)) then  !NAME
        call Str_Upper(W)  ;  M%Name= trim(W)
      end if

      if(I==vifield(4)) then  !SCFORM
        !
        call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        if(.not.EcformIsOk) cycle DoFile !-------------------------cycle
        !
        call Str_Upper(W)  ;  M%Formula=trim(W)
      end if

      if(I==vifield(5)) then  !ECFORM
        call Str_Upper(W)  ;  M%Formula=  trim(W)
      end if

      if(I==vifield(9)) then  !FITTING
        call Str_Upper(W)  ;  CodFitting= trim(W)
      end if

      if(I==vifield(11)) call WrdToReal(W,X) ! size
      if(I==vifield(12)) call WrdToReal(W,X) ! VOLUME
      if(I==vifield(13)) call WrdToReal(W,X) ! DENSITY

    end do
    !
    call ReadRValsV(L,mDum,vX)
    !-------------------------------------------------/ scan a data line
    !
    ! if(vifield(2)/=0)  print *, "M%Num=      ",trim(M%Num)
    ! if(vifield(3)/=0)  print *, "M%Name=     ",trim(M%Name)
    ! if(vifield(4)/=0)  print *, "M%Formula=  ",trim(M%Formula)
    ! if(vifield(5)/=0)  print *, "M%Formula=  ",trim(M%Formula)
    ! if(vifield(9)/=0)  print *, "CodFitting= ",trim(CodFitting)
    ! call pause_
    !
    !---------------------------- check whether Name not already read --
    if(LnkSpc_Found(LisLogKAnl,M%Name)) then
      if(iDebug>0) write(fTrc,'(A)') &
      & trim(M%Name)//"= name already used -> SKIPPED"
      !print *, trim(M%Name)//" <-this name already used by another species ???"
      !call pause_
      cycle DoFile
    end if
    !---/
    !
    !------------------------------------------------------ read formula
    call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    if(.not. fOk) cycle DoFile !--============================< cycle ==
    !-----------------------------------------------------/ read formula
    !
    !-- compute Molecular Weit
    M%WeitKg= dot_product(vStoik(:),vEle(:)%WeitKg) /real(Div)
    M%Div= Div
    M%Chg= Zsp
    !
    !------------size parameter for aqu'species, density for min'species
    if(    vifield(11)/=0 & ! size
    & .or. vifield(12)/=0 & ! VOLUME, m^3 /mole ?
    & .or. vifield(13)/=0 & ! DENSITY, Kg /m^3
    & ) then
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
    !-------- default value of size parameter for charged aqu'species --
    if(M%Typ=="AQU" .and. M%Chg/=0 .and. M%AquSize<=Zero) M%AquSize= 3.72D0
    !-----------/size parameter for aqu'species, density for min'species
    !
    !----------------------------------------------read function coeff's
    select case(trim(CodFitting))
      case("FIXED","none")      ;  M%iFitting= 0
      case("PHREEQC")           ;  M%iFitting= 1
      case("CHRISTOV")          ;  M%iFitting= 2
      case default
        call Stop_(trim(W)//" invalid Fitting mode")
    end select

    select case(M%iFitting)
    case(0)
      M%vX(1)= vX(1)
    case(1)
      if(M%Typ=="AQU") then  ;  M%vX(1:6)=  vX(1:6)
      else                   ;  M%vX(1:6)= -vX(1:6)
      ! for PHREEQC or THERMODDEM analytical,
      ! logK of MINERALS is for DISSOCIATION
      end if
    case(2)
      M%vX(1:6)= vX(1:6)
    end select
    !---------------------------------------------/read function coeff's
    !
    !write(fff,'(A,A1)') trim(L) !-> write the logK's
    !
    !M%Fitting= trim(CodFitting)
    !
    if(iDebug>2) call Test(ff,M) !---------------------------------trace
    !
    !---------------------------------------------"fill" the linked list
    N=N+1
    if(N==1) then
      allocate(LisLogKAnl)
      nullify(LisLogKAnl%next)
      LisLogKAnl%Value= M
      LisCur=> LisLogKAnl
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
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(iDebug>2) close(ff)!------------------------------------------trace
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbLogKAnl_Read"
  !
  if(N==0) call Stop_("NO SPECIES FOUND ....")
  !
end subroutine DtbLogKAnl_Read

subroutine FieldList_Read( &
& L,          &
& vStrField, &
& vifield)
  
  use M_IoTools
  
  !---------------------------------------------------------------------
  character(len=*),intent(inout):: L
  character(len=12),intent(in):: vStrField(:)
  integer,intent(out):: vifield(:)
  !---------------------------------------------------------------------
  character(len=80):: W
  logical:: EoL
  integer:: I,J
  !---------------------------------------------------------------------
  
  vifield(:)= 0
  I=0
  do
    call LinToWrd(L,W,EoL)
    I=I+1
    do J=1,size(vStrField)
      if( trim(W)==trim(vStrField(J)) ) then
        vifield(J)= I
        exit
      end if
    end do
    if(EoL) exit
  end do

  if(iDebug==4) then
    do I=1,size(vStrField)
      print *,vifield(I),trim(vStrField(I))
    end do
  end if
  !call pause_
  !
  if(vifield(10)==0) & ! for "PARAMETERS"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(10)))

  if(vifield(3)==0) & ! for "NAME"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(3)))

  if(vifield(4)==0 .and. vifield(5)==0) & ! for ECFORM/SCFORM
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(4))//"_"//trim(vStrField(5)))
  !
  !print *,vifield(1),vifield(2),vifield(3),vifield(4),vifield(5),vifield(9),vifield(10)
  !call pause_
end subroutine FieldList_Read

subroutine Test(ff,M)
  use M_T_Species
  use M_Numeric_Const,only: Ln10
  use M_T_DtbLogKAnl
  integer,intent(in):: ff
  type(T_DtbLogKAnl),intent(in) :: M
  !
  type(T_Species):: S
  real(dp):: vT(1:5),vLogK(1:5),Pbar
  integer :: I
  vT(1:5)= (/0.01,25.,50.,75.,100./)
  Pbar= 1.0D0
  do I=1,5
    call DtbLogKAnl_Calc(M,vT(I)+273.15D0,Pbar,S)
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

end module M_Dtb_Read_DtbLogKAnl







