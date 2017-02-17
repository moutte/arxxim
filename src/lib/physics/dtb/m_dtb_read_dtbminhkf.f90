module M_Dtb_Read_DtbMinHkf
  use M_Kinds
  use M_IoTools
  use M_Trace
  use M_Dtb_Read_Vars
  !
  implicit none
  !
  private
  !
  public:: DtbMinHKF_Build
  public:: DtbMinHKF_Read_Old
  public:: DtbMinHKF_Read
  !
  type(T_LisMinHkf),pointer,save:: LisCur, LisPrev

contains

subroutine DtbMinHKF_Build
!--
!-- transfer linked list to array of T_DtbMinHkf
!--
  use M_Dtb_Vars,only: vDtbMinHkf
  !
  integer:: I
  !
  if(iDebug>5) print '(A)',"< DtbMinHKF_Build"
  !
  if(allocated(vDtbMinHkf)) deallocate(vDtbMinHkf)
  allocate(vDtbMinHkf(nMinHkf))
  !
  LisCur=>LisMinHkf
  !
  I=0
  do while (associateD(LisCur))
    I=I+1
    vDtbMinHkf(I)=LisCur%Value
    if(iDebug>5) print '(2A)', vDtbMinHkf(I)%Name,trim(vDtbMinHkf(I)%Num) !; call pause_
    LisPrev=>LisCur
    LisCur=> LisCur%next
    deallocate(LisPrev)
  end do
  !
  if(iDebug>5) print '(A)',"</ DtbMinHKF_Build"
  if(iDebug>5) call Pause_
  return
end subroutine DtbMinHKF_Build

subroutine DtbMinHKF_Read(F,vEle,N)
!--
!-- reads data for minerals & gases directly from (MODIFIED) slop*.dat, or obigt.dat
!-- "new" formats are without subblocks for minerals / gas:
!-- SPECIES HSV.HKF
!--   MIN M001 FORSTERITE ../..
!--   MIN ../..
!--   GAS G001 CO2,G ../..
!--   GAS ../..
!-- end SPECIES
!--
!-- "old" formats have separate blocks for MINERAL and GAS:
!-- SPECIES MIN.HKF
!-- MINERAL
!--   M001 FORSTERITE ../..
!--   ../..
!-- end MINERAL
!-- GAS
!--   G001 CO2,G ../..
!--   ../..
!-- end GAS
!-- end SPECIES
!--
!-- NB !!! parameters for phase transition are not read !!!
!--
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_Dtb_Const,  only: CalToJoule,Tref,DtbConv_Benson
  use M_T_Element,  only: T_Element,Formula_Read,Formula_Build,Element_Index
  use M_Dtb_Read_Vars,only: LisMinHkf
  use M_T_DtbMinHKF
  use M_Dtb_Read_Tools
  !
  use M_Dtb_Vars,   only: DtbFormat
  !
  integer,intent(in):: F !input file
  type(T_Element),intent(in):: vEle(:)
  integer,intent(inout):: N
  !
  character(len=255):: L,W,sFormul
  character(len=4)  :: ICode
  character(len=10) :: sFormat
  type(T_DtbMinHkf) :: M
  real(dp),dimension(dimV)::vX
  logical :: EoL, fOk, SubBlock
  ! logical :: bMin, bGas
  integer :: I,ios,K,ZSp,FF,Div,iEl,ieO_
  real(dp):: H_Calc
  !
  integer,allocatable::vStoik(:) ! for Formula_Read
  !
  !--- for header processing --
  logical :: IsHeader
  integer,parameter:: nField= 10
  character(len=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  integer :: vifield(1:nField)
  ! character(len=12):: vStrUnit(1:nField)
  !---/ for header processing --
  !
  !--- for formula translation --
  character(len=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),allocatable:: vElement(:)
  logical:: EcformIsOk
  integer:: fFormula
  !---/
  !
  ! if (N>1) LisCur => LisMinHkf
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbMinHKF_Read"
  !
  !----------------------------------------------------------------trace
  FF= 0
  if(iDebug>2) then
    write(fTrc,'(A)') &
    & "stoikio in "//trim(DirDtbLog)//"minhkf_stoik.log"
    !
    call GetUnit(FF)
    open(FF,file=trim(DirDtbLog)//"minhkf_stoik.log")
    !
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"minhkf_stoik.log",&
    & "check species stoikio of min.species from supcrt-like database")
    !
    write(FF,"(A15,A1)",advance="NO") "NAME",T_
    do iEl=1,size(vEle)
      write(FF,"(A3,A1)",advance="NO") vEle(iEl)%NamEl,T_
    end do
    write(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  end if
  !---------------------------------------------------------------/trace
  !
  SubBlock= .false.
  sFormat= "HKF_CAL"
  !
  !--- for formula translation --
  CodFormula= "ECFORM"
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  !~ if(iDebug>3) write(fTrc,'(6(A,A1))') &
  !~ & "M%Name",T_, "M%G0R",T_, "M%S0Ele",T_, "M%H0R",T_, "H_Calc",T_, "M%H0R-H_Calc",T_
  !
  allocate(vStoik(1:size(vEle)))
  ieO_= Element_Index("O__",vEle)
  !
  FilCode=trim(DtbFormat) !default, also read from file
  !
  !----------------------------------initialize vStrField(:), vifield(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  !~ if (FilCode(1:5)=="OBIGT") then
    !~ L= "SOURCE NAME SKIP SCFORM SKIP SKIP SKIP SKIP parameterS"
  !~ else ! case of SLOP98.DAT
    !~ L= "SOURCE NAME SKIP SKIP ECFORM SKIP SKIP parameterS"
  !~ end if
  !~ !
  !~ !if(iDebug==4) print *,"< Default values"
  !~ call FieldList_Read(L,vStrField,vifield)
  !~ !if(iDebug==4) print *,"</ Default values"
  !---/ scan default field list
  !
  if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
  ! -> files contains compact formulas
    call GetUnit(fFormula)
    open(fFormula,file="debug_formula.log")
    write(fFormula,'(A,/)') "resuts of formula conversion"
  end if
  !---------------------------------/initialize vStrField(:), vifield(:)
  !
  DoFile: do
    !
    call DtbMinHkf_Zero(M)
    !
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!") cycle DoFile
    call AppendToEnd(L,W,EoL)
    !
    !------------------------------------------- read header, if present
    !--------- if the line begins with any member of array vStrField(:),
    !--------------------------------- then it contains the line headers
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
      !if(iDebug==4) print *,"< values from file"
      call FieldList_Read(L,vStrField,vifield)
      !if(iDebug==4) print *,"</ values from file"

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
    case("ENDSPECIES")
      exit DoFile
    !
    case("END")
      exit DoFile
    !
    !----------------------------------------------------for old formats
    !old! case("END")
    !old!   if(SubBlock) then
    !old!     SubBlock= .false.
    !old!     cycle DoFile
    !old!   else
    !old!     exit DoFile
    !old!   end if
    !old! case("ENDMINERAL","ENDGAS")
    !old!   SubBlock= .false.
    !old!   cycle DoFile
    !old! !
    !old! case("MINERAL")
    !old!   SubBlock= .true.
    !old!   M%Typ= "MIN"
    !old!   cycle DoFile
    !old! case("GAS")
    !old!   SubBlock= .true.
    !old!   M%Typ= "GAS"
    !old!   cycle DoFile
    !old! case("CODE")
    !old!   !-> can change the CODE, i.e. the data source inside the block
    !old!   call LinToWrd(L,W,EoL)
    !old!   FilCode= trim(W)
    !old!   !
    !old!   ! if FilCode changes, must re-initialize vifield ...!
    !old!   if (FilCode(1:5)=="OBIGT") then
    !old!     L= "INDEX NAME SKIP SCFORM type SKIP SKIP SKIP parameterS"
    !old!   else ! case of SLOP98.DAT
    !old!     L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP parameterS"
    !old!   end if
    !old!   call FieldList_Read(L,vStrField,vifield)
    !old!   !
    !old!   cycle DoFile
    !----------------------------------------------/for old formats --
    !
    case default
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= trim(W)//" "//trim(L)
      !
    end select
    !------------------------------------/ process first word of line --
    !
    M%Num=trim(FilCode)
    !
    !--------------------------------------------- scan the data line --
    !--------------------------------- up to column before parameterS --
    do I= 1,vifield(10)-1
      !
      !vStrField(1:nField)= &
      !!!!!"____________","____________","____________","____________","____________",
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
      !
      call LinToWrd(L,W,EoL,"NO")
      !
      if(EoL) cycle DoFile ! data line, contains no sufficient data -> skip !!
      !
      if(I==vifield(7)) then  !CODE / SOURCE
        call Str_Upper(W)  ;  FilCode= trim(W)
      end if

      if(I==vifield(8)) then  !format
        call Str_Upper(W)  ;  sFormat= trim(W)
      end if

      if(I==vifield(1)) then  !type
        !old! call Str_Upper(W)  ;  M%Typ= trim(W)
        call Str_Upper(W)
        select case(trim(W))
        ! case("AQU")  ;  M%Typ="AQU"
        case("MIN")  ;  M%Typ="MIN"
        case("GAS")  ;  M%Typ="GAS"
        case("LIQ")  ;  M%Typ="LIQ"
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

    end do
    !------------------------ scan the data line (left to parameterS) --
    !
    !--------------------------------------------------- read formula --
    call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    if(.not. fOk) cycle DoFile !--============================< cycle ==
    !--------------------------------------------------/ read formula --
    !
    M%S0Ele=  dot_product(vStoik(:), vEle(:)%S0) /real(Div) !in Joule
    M%WeitKg= dot_product(vStoik(:), vEle(:)%WeitKg) /real(Div)
    M%Div= Div
    !
    call ReadRValsV(L,K,vX)
    !
    if (FilCode(1:5)=="OBIGT") then
      ! in OBIGT database, Cp is tabulated fot min'species
      ! but not used as input parameters -> vX(4) ignored
      M%G0R= vX(1)
      M%H0R= vX(2)
      M%S0_= vX(3)
      M%V0R= vX(5)
      M%MK1(1)= vX(6)
      M%MK1(2)= vX(7) *1.0D-3
      M%MK1(3)= vX(8) *1.0D5
      M%NTran= 0
    else
      ! default FilCode is SLOP
      M%G0R=    vX(1)
      M%H0R=    vX(2)
      M%S0_=    vX(3)
      M%V0R=    vX(4)
      M%MK1(1)= vX(5)
      M%MK1(2)= vX(6) *1.0D-3
      M%MK1(3)= vX(7) *1.0D5
      M%NTran= 0
    end if
    !
    if(sFormat(1:7)=="HKF_CAL") then
      M%G0R= M%G0R*CalToJoule
      M%H0R= M%H0R*CalToJoule
      M%S0_= M%S0_*CalToJoule
      M%MK1(1)= M%MK1(1)
      M%MK1(2)= M%MK1(2)*CalToJoule
      M%MK1(3)= M%MK1(3)*CalToJoule
    end if
    !
    N=N+1
    !
    call IntToStr4(N,ICode)
    M%Num= trim(FilCode)//"_"//trim(ICode)
    !
    if(.not. DtbConv_Benson) M%G0R= M%G0R - Tref *M%S0Ele !!!Berman Convention!!!
    !
    !----------------------------------------- "fill" the linked list --
    !call Save_Record(N,M)
    if(N==1) then
      allocate(LisMinHkf)
      nullify(LisMinHkf%next)
      LisMinHkf%Value= M
      LisCur=> LisMinHkf
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    end if
    !-----------------------------------------/"fill" the linked list --
    !
    H_Calc= M%G0R + Tref *M%S0_ - Tref*M%S0Ele
    !~ if(iDebug>3) write(fTrc, '(A,A1,5(G15.6,A1))') &
    !~ & M%Name,T_, M%G0R,T_, M%S0Ele,T_, M%H0R,T_, H_Calc,T_, M%H0R-H_Calc,T_ ! M%H0R-H_Calc
    !
    if(FF>0) then
      write(FF,"(A15,A1)",advance="NO") M%Name,T_
      do iEl=1,size(vEle)
        write(FF,"(I7,A1)",advance="NO") vStoik(iEl),T_
      end do
      write(FF,'(4(G15.6,A1))',advance="NO") &
      & M%G0R,T_,H_Calc,T_,M%H0R,T_,M%S0_,T_
      !
      ! build new formula, with fixed element order
      call Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!do iEl=1,size(vEle)
      !!  if(vStoik(iEl)>0) sFormul=trim(sFormul)//trim(vEle(iEl)%NamEl)//""
      !!end do
      write(FF,"(I3,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,trim(sFormul)
    end if
    !
  end do DoFile
  !
  deallocate(vStoik)
  !
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(FF>0) close(FF)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbMinHKF_Read"
  !
  return
end subroutine DtbMinHKF_Read

subroutine Save_Record(N,M)
  integer,          intent(in):: N
  type(T_DtbMinHkf),intent(in):: M
  !
    if(N==1) then
      allocate(LisMinHkf)
      nullify(LisMinHkf%next)
      LisMinHkf%Value= M
      LisCur=> LisMinHkf
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    end if
  !
end subroutine Save_Record

subroutine DtbMinHKF_Read_Old(F,vEle,N)
!--
!-- reads data for minerals directly from (MODIFIED) slop*.dat, or obigt.dat
!-- this is for "old" formats, with separate blocks for MINERAL and GAS:
!-- SPECIES MIN.HKF
!-- MINERAL
!--   ..
!-- end MINERAL
!-- GAS
!--   ..
!-- end GAS
!-- end SPECIES
!-- NB !!! parameters for phase transition are not read !!!
!--
  use M_Dtb_Const,  only: CalToJoule,Tref,DtbConv_Benson
  use M_T_Element,  only: T_Element,Formula_Read,Formula_Build,Element_Index
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_Dtb_Vars,   only: DtbFormat
  use M_Dtb_Read_Tools
  !
  integer,intent(in):: F !input file
  type(T_Element),intent(in):: vEle(:)
  integer,intent(inout):: N
  !
  !type(T_LisMinHkf),pointer,save:: LisCur
  !
  character(len=255):: L,W,sFormul
  character(len=4)  :: ICode
  type(T_DtbMinHkf) :: M
  real(dp),dimension(dimV)::vX
  logical :: EoL, fOk
  logical :: bMin, bGas
  integer :: ios,K,ZSp,FF,Div,iEl,ieO_
  real(dp):: H_Calc
  !
  integer,allocatable::vStoik(:) !for Formula_Read
  !
  !--- for formula translation --
  character(len=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),allocatable:: vElement(:)
  logical:: EcformIsOk
  integer:: fFormula
  !---/
  !
  ! if (N>1) LisCur => LisMinHkf
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbMinHKF_Read"
  !
  !----------------------------------------------------------------trace
  FF= 0
  if(iDebug>2) then
    write(fTrc,'(A)') &
    & "stoikio in "//trim(DirDtbLog)//"minhkf_stoik.log"
    !
    call GetUnit(FF)
    open(FF,file=trim(DirDtbLog)//"minhkf_stoik.log")
    !
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"minhkf_stoik.log",&
    & "check species stoikio of min.species from supcrt-like database")
    !
    write(FF,"(A15,A1)",advance="NO") "NAME",T_
    do iEl=1,size(vEle)
      write(FF,"(A3,A1)",advance="NO") vEle(iEl)%NamEl,T_
    end do
    write(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  end if
  !---------------------------------------------------------------/trace
  !
  !--- for formula translation --
  CodFormula= "ECFORM"
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  !~ if(iDebug>3) write(fTrc,'(6(A,A1))') &
  !~ & "M%Name",T_, "M%G0R",T_, "M%S0Ele",T_, "M%H0R",T_, "H_Calc",T_, "M%H0R-H_Calc",T_
  !
  allocate(vStoik(1:size(vEle)))
  ieO_= Element_Index("O__",vEle)
  !
  FilCode=trim(DtbFormat) !default, also read from file
  !
  DoFile: do

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!") cycle DoFile
    call AppEndToEnd(L,W,EoL)
    !
    select case(W)
      !
      case("END","ENDSPECIES","ENDINPUT")
        exit DoFile
        !
      case("CODE")
        call LinToWrd(L,W,EoL)
        FilCode=trim(W)
        cycle DoFile
        !
      case("FORMULA")
        call LinToWrd(L,W,EoL)
        CodFormula= trim(W)
        if(iDebug>2 .and. CodFormula=="SCFORM") then
          call GetUnit(fFormula)
          open(fFormula,file="debug_formula.log")
          write(fFormula,'(A,/)') "resuts of formula conversion"
        end if
        cycle DoFile
        !
      case("MINERAL","GAS")
        !
        select case(W)
          case("MINERAL"); bGas=.false.; bMin=.true.
          case("GAS");     bGas=.true.;  bMin=.false.
        end select
        !
        DoReadMin: do
          !
          read(F,'(A)',iostat=ios) L
          if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=="!") cycle DoReadMin
          call AppendToEnd(L,W,EoL)
          !
          select case(W)
            case("ENDINPUT","ENDSPECIES")
              exit DoFile
            case("END","ENDMINERAL","ENDGAS")
              cycle DoFile !exit DoReadMin
          end select
          !
          M%Num=trim(W)
          if(bGas) M%Typ="GAS"
          if(bMin) M%Typ="MIN"
          !
          if (FilCode(1:5)=="OBIGT") then
            !
            !L= "SOURCE NAME SKIP SCFORM SKIP SKIP SKIP SKIP parameterS"
            CodFormula="SCFORM"
            !
            call LinToWrd(L,W,EoL); M%Name=trim(W) !if(iDebug>0) write(fTrc,"(A)") M%Name
            call LinToWrd(L,W,EoL) !skip ABBRV
            call LinToWrd(L,W,EoL,"NO") !-> compact formula, character case is conserved :!!
            !
            if(CodFormula=="SCFORM") then
              call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
              if(.not.EcformIsOk) then
                if(iDebug>3) call Warning_("!!! Cannot translate "//trim(W))
                cycle DoReadMin
              end if
            end if
            call Str_Upper(W)
            !
            M%Formula=trim(W)
            !
            call LinToWrd(L,W,EoL) !skip STATE
            call LinToWrd(L,W,EoL) !skip SOURCE1
            call LinToWrd(L,W,EoL) !skip SOURCE2
            call LinToWrd(L,W,EoL) !DATE
            !
            call ReadRValsV(L,K,vX)
            !
          else ! default FilCode is SLOP98
            !
            !L= "SOURCE SKIP NAME SKIP ECFORM SKIP SKIP parameterS"
            !
            call LinToWrd(L,W,EoL) !skip abredged name ; M%Abbr=trim(W)
            call LinToWrd(L,W,EoL) ; M%Name=trim(W)
            call LinToWrd(L,W,EoL) !skip scform
            call LinToWrd(L,W,EoL) ; M%Formula=trim(W)
            call LinToWrd(L,W,EoL) !skip REF
            call LinToWrd(L,W,EoL) !skip DATE
            !
            call ReadRValsV(L,K,vX)
            !
          end if
          !
          !-------------------------------------------------read formula
          call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
          if(.not. fOk) cycle DoReadMin !--------------------------cycle
          !------------------------------------------------/read formula
          call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
          !
          if (FilCode(1:5)=="OBIGT") then
            ! in OBIGT database, Cp is tabulated fot min'species
            ! but not used as input parameters -> vX(4) ignored
            M%G0R= vX(1) *CalToJoule
            M%H0R= vX(2) *CalToJoule
            M%S0_= vX(3) *CalToJoule
            M%V0R= vX(5)
            M%MK1(1)= vX(6) *CalToJoule
            M%MK1(2)= vX(7) *CalToJoule *1.0D-3
            M%MK1(3)= vX(8) *CalToJoule *1.0D5
            M%NTran= 0
          else ! default FilCode is SLOP98
            M%G0R=    vX(1) *CalToJoule
            M%H0R=    vX(2) *CalToJoule
            M%S0_=    vX(3) *CalToJoule
            M%V0R=    vX(4)
            M%MK1(1)= vX(5) *CalToJoule
            M%MK1(2)= vX(6) *CalToJoule *1.0D-3
            M%MK1(3)= vX(7) *CalToJoule *1.0D5
            M%NTran= 0
          end if

          !if(fOk) then
          N=N+1
          !
          call IntToStr4(N,ICode)
          M%Num= trim(FilCode)//"_"//trim(ICode)
          !
          M%S0Ele=  dot_product(vStoik(:),vEle(:)%S0) /real(Div) ! is in Joule !!
          M%WeitKg= dot_product(vStoik(:),vEle(:)%WeitKg) /real(Div)
          M%Div=    Div
          !
          if(.not. DtbConv_Benson) M%G0R= M%G0R - Tref *M%S0Ele !!!Berman Convention!!!
          !
          !--- "fill" the linked list
          if(N==1) then
            allocate(LisMinHkf)
            nullify(LisMinHkf%next)
            LisMinHkf%Value= M
            LisCur=> LisMinHkf
          else
            allocate(LisCur%next)
            nullify(LisCur%next%next)
            LisCur%next%Value=M
            LisCur=>LisCur%next
          end if
          !---/
          !
          H_Calc= M%G0R + Tref *M%S0_ - Tref*M%S0Ele
          !~ if(iDebug>3) write(fTrc, '(A,A1,5(G15.6,A1))') &
          !~ & M%Name,T_, M%G0R,T_, M%S0Ele,T_, M%H0R,T_, H_Calc,T_, M%H0R-H_Calc,T_ ! M%H0R-H_Calc
          !
          if(FF>0) then
            write(FF,"(A15,A1)",advance="NO") M%Name,T_
            do iEl=1,size(vEle)
              write(FF,"(I3,A1)",advance="NO") vStoik(iEl),T_
            end do
            write(FF,'(4(G15.6,A1))',advance="NO") &
            & M%G0R,T_,H_Calc,T_,M%H0R,T_,M%S0_,T_
            !
            ! build new formula, with fixed element order
            call Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
            !!sFormul=""
            !!do iEl=1,size(vEle)
            !!  if(vStoik(iEl)>0) sFormul=trim(sFormul)//trim(vEle(iEl)%NamEl)//""
            !!end do
            write(FF,"(I3,A1,A39,A1,A15,A1,A)") &
            & M%Div,T_,M%Formula,T_,M%Name,T_,trim(sFormul)
          end if
          !  !
          !end if
        end do DoReadMin
      !endcase("MINERAL")
    end select !case(W0)
  end do DoFile
  !
  deallocate(vStoik)
  !
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(FF>0) close(FF)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbMinHKF_Read"
  !
  return
end subroutine DtbMinHKF_Read_Old

subroutine FieldList_Read( &
& L,          &
& vStrField, &
& vifield)
  use M_IoTools
  !
  character(len=*),intent(inout):: L
  character(len=12),intent(in):: vStrField(:)
  integer,intent(out):: vifield(:)
  !
  character(len=80):: W
  logical:: EoL
  integer:: I,J
  !
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
  !
  !! if(iDebug==4) then
  !!   do I=1,size(vStrField)
  !!     print *,vifield(I),trim(vStrField(I))
  !!   end do
  !! end if
  !call pause_
  !
  ! (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  !   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  if(vifield(10)==0) & ! for "PARAMETERS"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(10)))

  if(vifield(1)==0) & ! for "TYPE"  !!MIN/GAS/AQU
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(1)))

  if(vifield(3)==0) & ! for "NAME"
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(3)))

  if(vifield(4)==0 .and. vifield(5)==0) & ! for ECFORM/SCFORM
  call Stop_( &
  & "in FieldList_Read: keyword not found for "//trim(vStrField(4))//"_"//trim(vStrField(5)))
  !
end subroutine FieldList_Read

end module M_Dtb_Read_DtbMinHkf
