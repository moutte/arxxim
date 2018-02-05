module M_Dtb_Read_DtbMinThr
  use M_Kinds
  use M_IoTools
  use M_Trace
  use M_Dtb_Read_Vars
  !
  implicit none
  !
  private
  !
  public:: DtbMinThr_Build
  public:: DtbMinThr_Read
  public:: DtbMinThr_Read_Line
  public:: DtbMinThr_ToLine
  !
  type(T_LisMinThr),pointer,save:: LisCur, LisPrev

contains

subroutine DtbMinThr_Build
!--
!-- -> build vDtbMinThr
!-- transform LinkedList LisMinThr to mineral table vDtbMinThr
!--
  use M_T_Element,  only: T_Element,Formula_Read
  use M_Dtb_Vars,   only: vDtbMinThr
  !
  integer:: I
  !
  if(iDebug==5) print '(A)',"< DtbMinThr_Build"
  !
  if(allocated(vDtbMinThr)) deallocate(vDtbMinThr)
  allocate(vDtbMinThr(nMinThr))
  LisCur=>LisMinThr
  I=0
  do while (associateD(LisCur))
    I=I+1
    vDtbMinThr(I)=LisCur%Value
    if(iDebug==5) &
    & print '(2A)', vDtbMinThr(I)%Name,trim(vDtbMinThr(I)%Num) !; pause
    LisPrev=>LisCur
    LisCur=> LisCur%next
    deallocate(LisPrev)
  end do
  !
  if(iDebug==5) print '(A)',"</ DtbMinThr_Build"
  if(iDebug==5) call Pause_
  !
  return
end subroutine DtbMinTHR_Build

subroutine DtbMinThr_Read(F,vEle,N)
!--
!-- -> build or Append LisMinThr
!-- reads thermodyn. data from (MODIFIED) theriak-compatible datafile
!--
  use M_T_Element,  only: T_Element,Element_Index
  use M_T_DtbMinThr,only: DtbMinThr_Zero
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_Dtb_Read_Tools
  use M_Dtb_Vars,   only: DtbFormat
  !
  integer,intent(in):: F !input file
  type(T_Element),intent(in):: vEle(:)
  integer,intent(inout):: N
  !
  character(len=512):: L,W,Wb
  logical :: EoL
  logical :: bMin, bGas
  real(dp):: X
  integer :: ios,I,iEl,FF,ieO_
  type(T_DtbMinThr) :: M
  !
  character(len=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),dimension(:),allocatable:: vElement
  integer:: fFormula
  logical:: EcformIsOk
  !
  if(idebug>1) write(fTrc,"(/,A)") "< DtbMinThr_Read"
  !
  CodFormula= "ECFORM" ! is default value
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  !
  fFormula= 0
  !
  !--------------------------------------------------------------- trace
  if(iDebug>2) then
    write(fTrc,'(A)') &
    & "stoikio in "//trim(DirDtbLog)//"min_thr_stoik.log"
    call GetUnit(FF)
    open(FF,file=trim(DirDtbLog)//"min_thr_stoik.log")
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"min_thr_stoik.log",&
    & "check species stoikio of min.species from theriak-type database")
    write(FF,"(A15,A1)",advance="NO") "NAME",T_
    do iEl=1,size(vEle)
      write(FF,"(A3,A1)",advance="NO") vEle(iEl)%NamEl,T_
    end do
    write(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
  end if
  !--------------------------------------------------------------/ trace
  !
  FilCode=trim(DtbFormat) !-> default value
  ieO_= Element_Index("O__",vEle)
  !
  M%Name= "ZZZ"
  M%Num=  "0"
  !
  DoFile: do
    
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!") cycle DoFile
    call AppendToEnd(L,W,EoL)
    
    select case(W)
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
    case("CODE")
      call LinToWrd(L,W,EoL)
      FilCode=trim(W)
      cycle DoFile
    !_
    case("ENDINPUT"); exit DoFile
    !
    case("END","ENDSPECIES"); exit DoFile
    !
    case("MINERAL","GAS")
    
      select case(W)
      case("MINERAL") ; bGas=.false. ; bMin=.true.
      case("GAS")     ; bGas=.true.  ; bMin=.false.
      end select

      call DtbMinThr_Zero(M)
      
      DoReadMin: do
        
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        
        if(W(1:1)=="!") cycle DoReadMin !skip comment lines
        
        call AppendToEnd(L,W,EoL)
        
        select case(W)
        case("ENDINPUT"); exit DoFile
        case("ENDSPECIES"); exit DoFile
        case("END","ENDMINERAL","ENDGAS"); cycle DoFile
        end select
        
        if(W(1:1)/="&") then
        !if first char is not "&", the line may contain name, formula, etc
        
          if(M%Name/="ZZZ") then !before proceeding,  save current M in list
            call Save_Record
            call DtbMinThr_Zero(M)
          end if
          if(bGas) M%Typ="GAS"
          if(bMin) M%Typ="MIN"
          M%Name=trim(W)
          
          !----------------------------------processing compact formulas
          if(CodFormula=="SCFORM") then
            call LinToWrd(L,W,EoL,"NO")
            !-> compact formula, character case is conserved :!!
            call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
            if(.not.EcformIsOk) then !cycle DoReadMin
              call Stop_("!!! Cannot translate "//trim(W))
            end if
            call Str_Upper(W)
            M%Formula=trim(W)
          else
            call LinToWrd(L,W,EoL)
            M%Formula= trim(W)
          end if
          !----------------------------------------------------------end
          
          call LinToWrd(L,W,EoL) !; M%Abbr=trim(W)
          !
          cycle DoReadMin
          
        end if
        
        if(W(1:1)=="&") then
        ! if first char is "&",
        ! then the line is the continuation of preceding, i.e. contains data
          
          call LinToWrd(L,W,EoL)
          
          select case(trim(W))
          
          case("ST")
            call ReadRVals5(L,M%G0R,M%H0R,M%S0_,M%V0R,X)
          
          case("C1","CP1")
            call ReadRVals5(L,M%K1,M%K4,M%K3,M%K8,X)
            !print *,"M%K8,M%K9",M%K8,M%K9
            !K1,K4,K3,K8
          case("C2","CP2")
            call ReadRVals5(L,M%K6,M%K2,M%K5,M%K7,M%K9)
            !Cp= K6/T + K2*T + K5*TT + K7*SQT + K9*TT*T
          case("C3","CP3")
              call ReadRVals5(L,M%K1,M%K2,M%K3,M%K4,X)
          
          case("V1")
            M%codVol=1
            call ReadRVals5(L,M%VTA,M%VTB,M%VPA,M%VPB,X)
            M%VTA=M%VTA*M%V0R/1.0D5 ; M%VTB=M%VTB*M%V0R/1.0D5
            M%VPA=M%VPA*M%V0R/1.0D5 ; M%VPB=M%VPB*M%V0R/1.0D8
          case("VHP")
            M%codVol=4
            call ReadRVals5(L,M%VTA,M%VTB,M%TKRI,M%SMA,M%VPA)
          case("VH2")
            M%codVol=4
            call ReadRVals5(L,M%D1, M%D2, M%D3,X,X)
          
          case("D1")
            M%DIS=.true.
            call ReadRVals5(L,M%D1, M%D4, M%D3, M%D8, M%D6)
          case("D2")
            call ReadRVals5(L,M%D2, M%D5, M%TD0,M%TDMAX,M%VADJ)
          
          case("TR1","T1","TR3","T3") !for Berman dtbase"s
            if(M%nLanda<2) then
              M%nLanda=M%nLanda+1
            else
              call Stop_("DtbMinThr_Read: nLanda>MaxValue")
            end if
            I=M%nLanda
            call ReadRVals5(L,M%TQ1B(I),M%TRE(I),M%ASPK(I),M%BSPK(I),M%DHTR(I))
          case("TR2","T2")
            if (M%nLanda>0) then
              I=M%nLanda
              call ReadRVals5(L,M%TEQ(I),M%DVTR(I),M%DVDT(I),M%DVDP(I),X)
            end if
          
          case("REDKWON","VDWAALS")
          ! Redlich-Kwong, VanDerWaals
            M%CodGas=trim(W)
            call ReadRVals5(L,M%AA0,M%AAT,M%BB0,M%BBT,X)
          case("SOAVE","PENGROB")
          ! Soave Cubic, Peng Robinson cubic
            M%CodGas=trim(W)
            call ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)
            !!write(fTrc,'(A,3G15.6)') M%CodGas,M%TCrit,M%PCrit,M%ACentric
          case("PRSV")
          ! 
            M%CodGas="PRSV"
            call ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)
          
          case("SPECIAL")
            call LinToWrd(L,Wb,EoL); M%Special= trim(Wb)
          
          !case("V2"); M%codVol=2; call ReadRVals6(L,M%VAA,M%VAB,M%VB,            X,X,X)  !VO2=.true.
          !case("V3"); M%codVol=3; call ReadRVals6(L,M%VL0,M%VLA,M%VLN, M%VL2,     X,X)   !VO3=.true.
          !case("AQ1") K1,K2,K8,K9,K7  AQU=.true. PROVID="AQU"
          !case("AQ2") D1,D2,D3,D4
          !case("AQP") D1,D2           AQU=.true. PROVID="AQP"
          !case("SPC") case            SPC=.true.
          !case("TL1") then; M%TL1=.true.; call ReadRVals6(L,M%TKRI,M%SMA,X,X,X,X); end if !!!obsolete ???
          
          end select
          
        end if !(W(1:1)=="&")
      
      end do DoReadMin
      
    !endcase("MINERAL","GAS")
      
    end select
    
  end do DoFile
  !
  if(M%Name/="ZZZ") call Save_Record !Save last record
  !
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(iDebug>2) close(FF)
  if(idebug>1) write(fTrc,"(A,/)") "</ DtbMinThr_Read"
  
  return
  !
contains

subroutine Save_Record
  use M_T_Element,only: Formula_Read,Formula_Build
  !
  character(len=4)  :: ICode
  character(len=512):: sFormul
  integer,dimension(1:size(vEle))::vStoik !for Formula_Read
  logical :: fOk
  integer :: ZSp,Div
  !
  call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
  if(fOk) then
    !
    N=N+1
    call IntToStr4(N,ICode)
    M%Num=trim(FilCode)//"_"//trim(ICode)
    !
    M%S0Ele= dot_product(vStoik(:), vEle(:)%S0) /real(Div) !in Joule
    M%WeitKg=dot_product(vStoik(:), vEle(:)%WeitKg) /real(Div)
    M%Div= Div
    !
    if(N==1) then
      allocate(LisMinThr)
      nullify(LisMinThr%next)
      LisMinThr%Value=M
      LisCur=> LisMinThr
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=> LisCur%next
    end if
    !
    if(iDebug>2) then
      write(FF,"(A15,A1)",advance="NO") M%Name,T_
      do iEl=1,size(vEle)
        write(FF,"(I5,A1)",advance="NO") vStoik(iEl),T_
      end do
      !build new formula, with fixed element order
      call Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!do iEl=1,size(vEle)
      !!  if(vStoik(iEl)>0) sFormul=trim(sFormul)//trim(vEle(iEl)%NamEl)//""
      !!end do
      write(FF,"(I5,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,trim(sFormul)
    end if
    !
  end if
end subroutine Save_Record

end subroutine DtbMinThr_Read

subroutine DtbMinThr_Read_Line(F,vEle,N)
!--
!---> build or Append LisMinThr
!--
!-- reads thermodyn. data from theriak-compatible single line datafile
!-- this single line reading is called when DtbFormat is "HSV.THR"
!-- (multiline is called for "MIN.THR" or "GAS.THR")
!-- differences from the multiline format:
!--   1/ whether the species is MIN or GAS is indicated for each species,
!--   and not for a whole block
!--   2/ columns on left are free format
!-----------------------------------------------------------------------
!-- first column- AQU/MIN/GAS
!-- then columns for species description (Name, Formula, Data source, ...),
!-- then one or more blocks of 6 columns,
!-- each with a code (ST,C1,...) on column 1 and data on columns 2..6
!-----------------------------------------------------------------------

  use M_T_Element,  only: T_Element,Element_Index
  use M_T_Element,  only: Formula_Read
  use M_Dtb_Const,  only: Tref
  use M_T_DtbMinThr,only: DtbMinThr_Zero
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_Dtb_Read_Tools
  use M_Dtb_Vars,   only: DtbFormat
  !---------------------------------------------------------------------
  integer,        intent(in):: F !input file
  type(T_Element),intent(in):: vEle(:)
  !
  integer,        intent(inout):: N
  !---------------------------------------------------------------------
  character(len=512):: L,W,Wb
  logical :: EoL
  real(dp):: X
  integer :: ios,I,iEl,FF,ieO_,iLine
  type(T_DtbMinThr) :: M
  !type(T_LisMinThr),pointer,save:: LisCur
  !
  !----------------------------------------------- for header processing
  logical :: IsHeader
  integer,parameter:: nField= 10
  character(len=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  integer :: vifield(1:nField)
  ! character(len=12):: vStrUnit(1:nField)
  !---/ for header processing
  !
  !---------------------------------------------------- for Formula_Read
  integer :: vStoik(1:size(vEle))
  logical :: fOk
  integer :: ZSp,Div
  !---/
  !
  !--------------------------------------------- for formula translation
  character(len=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),allocatable:: vElement(:)  !
  integer:: fFormula
  logical:: EcformIsOk
  !---/
  !---------------------------------------------------------------------
  if(idebug>1) write(fTrc,"(/,A)") "< DtbMinThr_Read_line"
  !
  !--- for formula translation --
  CodFormula= "ECFORM" ! is default value
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  ! if(N>0) LisCur=> LisMinThr
  !
  !--------------------------------------------------------------- trace
  if(iDebug<3) then
    FF= 0
  else
    !
    write(fTrc,'(A)') &
    & "stoikio in "//trim(DirDtbLog)//"min_thr_stoik.log"
    !
    call GetUnit(FF)
    open(FF,file=trim(DirDtbLog)//"min_thr_stoik.log")
    !
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"min_thr_stoik.log",&
    & "check species stoikio of min.species from theriak-type database")
    !
    write(FF,"(A15,A1)",advance="NO") "NAME",T_
    do iEl=1,size(vEle)
      write(FF,"(A3,A1)",advance="NO") vEle(iEl)%NamEl,T_
    end do
    write(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  end if
  !--------------------------------------------------------------/ trace
  !
  FilCode= trim(DtbFormat) !-> default value
  iLine= 0
  ieO_= Element_Index("O__",vEle)
  !
  !------------------------------ initialize vStrField(:), vifield(:) --
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "ABBREV      ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  L= "TYPE NAME ECFORM ABBREV SOURCE FORMAT PARAMETERS"
  L= "TYPE NAME ECFORM SKIP   SOURCE SKIP   PARAMETERS"
  !
  !if(iDebug==4) print *,"< Default values"
  call FieldList_Read(L,vStrField,vifield)
  call FieldList_Check(1,vStrField,vifield) !"TYPE"
  !if(iDebug==4) print *,"</ Default values"
  !---/ scan default field list
  !
  if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
  ! -> files contains compact formulas
    call GetUnit(fFormula)
    open(fFormula,file="debug_formula.log")
    write(fFormula,'(A,/)') "resuts of formula conversion"
  end if
  !--------------------------------/ initialize vStrField(:), vifield(:)
  !
  !---------------------------------- build a linked list of all species
  !-------------------------------- consistent with current element list
  DoFile: do

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
      call FieldList_Check(1,vStrField,vifield) !"TYPE"
      !if(iDebug==4) print *,"</ values from file"

      if(vifield(4)/=0 .and. fFormula==0 .and. iDebug>2) then
        call GetUnit(fFormula)
        open(fFormula,file="debug_formula.log")
        write(fFormula,'(A,/)') "results of formula conversion"
      end if

      cycle DoFile
    end if
    !------------------------------------------/ read header, if present
    !
    !---------------------------------------- process first word of line
    select case(W)
    !
    case("ENDINPUT")
      exit DoFile
    !
    case("END","ENDSPECIES")
      exit DoFile
    !
    case default
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= trim(W)//" "//trim(L)
    !
    end select
    !---------------------------------------/ process first word of line --
    !
    call DtbMinThr_Zero(M)
    !
    M%Num=trim(FilCode)
    !
    !------------------------------------------------ scan the data line
    !------------------------------------ up to column before parameters
    do I= 1,vifield(10)-1
      !
      !vStrField(1:nField)= &
      !!!!!"____________","____________","____________","____________","____________",
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "ABBREV      ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
      !
      call LinToWrd(L,W,EoL,"NO")
      !
      if(EoL) cycle DoFile ! data line, contains no sufficient data -> skip !!
      !
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
      !
      if(I==vifield(2)) then  !INDEX
        call Str_Upper(W)  ;  M%Num= trim(W)
      end if
      !
      if(I==vifield(3)) then  !NAME
        call Str_Upper(W)  ;  M%Name= trim(W)
      end if
      !
      if(I==vifield(4)) then  !SCFORM
        call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        if(.not.EcformIsOk) cycle DoFile !-------------------------cycle
        call Str_Upper(W)  ;  M%Formula= trim(W)
      end if
      !
      if(I==vifield(5)) then  !ECFORM
        call Str_Upper(W)  ;  M%Formula= trim(W)
      end if
      !
      if(I==vifield(7)) then ! SOURCE
        call Str_Upper(W)  ;  M%Source= trim(W)
      end if
      !
      if(I==vifield(8)) then ! FORMAT !---------------------NEW--2018-01
        call Str_Upper(W)  ;  M%HSV_Format= trim(W)
      end if
      !
    end do
    !---------------------------- scan the data line (left to PARAMETERS)
    !
    !-------------------------------------------------------read formula
    call Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    if(.not. fOk) cycle DoFile !-----------------------------------cycle
    !-----------------------------------------------------/ read formula
    !
    M%S0Ele=  dot_product(vStoik(:), vEle(:)%S0) /real(Div) !in Joule
    M%WeitKg= dot_product(vStoik(:), vEle(:)%WeitKg) /real(Div)
    M%Div= Div
    !
    !-------------------------------------------- loop on 6-cells groups
    do
    ! read species data
    ! (= organized as groups of 6 cells:
    !  1 group= 1 code (ST,C1,V1,...) followed by 5 numeric data cells)
    !
      call LinToWrd(L,W,EoL)
      if(EoL) exit

      select case(trim(W))

      case("SPECIAL")
        call LinToWrd(L,Wb,EoL)
        M%Special= trim(Wb)
        
      case("ST")
        call ReadRVals5(L,M%G0R,M%H0R,M%S0_,M%V0R,X)
        if(M%HSV_Format=="GSV") then
          M%H0R= M%G0R + Tref*M%S0_ - Tref*M%S0Ele
        end if

      case("C1","CP1")
        call ReadRVals5(L,M%K1,M%K4,M%K3,M%K8,X)
        ! K1,K4,K3,K8

      case("C2","CP2")
        call ReadRVals5(L,M%K6,M%K2,M%K5,M%K7,M%K9)
        ! Cp= K6/T + K2*T + K5*TT + K7*SQT + K9*TT*T

      case("C3","CP3")
        call ReadRVals5(L,M%K1,M%K2,M%K3,M%K4,X)
        ! Cp= K1 + K2*T + K3/TT + K4/SQT -> Maier-Kelley

      case("V1")
        M%codVol=1
        call ReadRVals5(L,M%VTA,M%VTB,M%VPA,M%VPB,X)
        !
        M%VTA=M%VTA*M%V0R/1.0D5  ;  M%VTB=M%VTB*M%V0R/1.0D5
        M%VPA=M%VPA*M%V0R/1.0D5  ;  M%VPB=M%VPB*M%V0R/1.0D8

      case("VHP")
        M%codVol=4
        call ReadRVals5(L,M%VTA,M%VTB,M%TKRI,M%SMA,M%VPA)

      case("VH2")
        M%codVol=4
        call ReadRVals5(L,M%D1, M%D2, M%D3,X,X)

      case("D1")
        M%DIS=.true.
        call ReadRVals5(L,M%D1, M%D4, M%D3, M%D8,   M%D6)

      case("D2")
        call ReadRVals5(L,M%D2, M%D5, M%TD0,M%TDMAX,M%VADJ)

      case("TR1","T1","TR3","T3") !for Berman database"s
        M%nLanda= M%nLanda +1
        I= M%nLanda
        call ReadRVals5(L,M%TQ1B(I),M%TRE(I),M%ASPK(I),M%BSPK(I),M%DHTR(I))

      case("TR2","T2")
        if (M%nLanda>0) then
          I=M%nLanda
          call ReadRVals5(L,M%TEQ(I),M%DVTR(I),M%DVDT(I),M%DVDP(I),X)
        end if

      case("REDKWON","VDWAALS")
        M%CodGas=trim(W)
        call ReadRVals5(L,M%AA0,M%AAT,M%BB0,M%BBT,X)

      case("SOAVE","PENGROB")
        M%CodGas= trim(W)
        call ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)

      case("PRSV")
        M%CodGas= "PRSV"
        call ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)

      case default
        call Stop_(trim(W)//"= unknown code in HSV/Theriak format")

      end select

    end do
    !-------------------------------------------/ loop on 6-cells groups
    !
    call Save_Record
    !
  end do DoFile
  !
  !! print *, "nMinThr=", N  ;  pause
  !
  deallocate(vElement)
  if(fFormula>0) close(fFormula)
  !
  if(FF>0) close(FF)
  if(idebug>1) write(fTrc,"(A,/)") "</ DtbMinThr_Read"
  !
  return
  !
contains

  subroutine Save_Record
  !--
  !-- save record in the linked list LisMinThr
  !--
    use M_T_Element,only: Formula_Build
    !
    character(len=4)  :: ICode
    character(len=512):: sFormul
    !
    N=N+1
    !
    iLine= iLine +1
    call IntToStr4(iLine,ICode)
    M%Num= trim(M%Num)//"_"//trim(ICode)
    !
    if(N==1) then
      allocate(LisMinThr)
      nullify(LisMinThr%next)
      LisMinThr%Value=M
      LisCur=> LisMinThr
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=> LisCur%next
    end if
    !
    if(FF>0) then !--=========================================< trace ==
      write(FF,"(A15,A1)",advance="NO") M%Name,T_
      do iEl=1,size(vEle)
        write(FF,"(I5,A1)",advance="NO") vStoik(iEl),T_
      end do
      ! build new formula, with fixed element order
      call Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!do iEl=1,size(vEle)
      !!  if(vStoik(iEl)>0) sFormul=trim(sFormul)//trim(vEle(iEl)%NamEl)//""
      !!end do
      write(FF,"(I5,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,trim(sFormul)
    end if !==================================================</ trace ==
    !
  end subroutine Save_Record

end subroutine DtbMinThr_Read_Line

subroutine ReadRVals5(Line,x1,x2,x3,x4,x5)
  use M_IOTools,only: LinToWrd,WrdToReal
  !
  character(len=*),intent(inout):: Line
  real(dp),        intent(out)  ::x1,x2,x3,x4,x5
  !
  character(255)   Word
  integer  ::i
  logical  ::EoL
  real(dp),dimension(1:5)::vX
  !
  vX(:)= Zero
  EoL=  .true.
  do i=1,5
    call LinToWrd(Line,Word,EoL)
    call WrdToReal(Word,vX(i))
    if(EoL) exit
  end do
  !
  x1=vX(1)
  x2=vX(2)
  x3=vX(3)
  x4=vX(4)
  x5=vX(5)
  !
end subroutine ReadRVals5

subroutine DtbMinThr_ToLine(f)
!--
!-- write a multiline theriak-compatible datafile to a monoline equivalent
!--
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_Dtb_Vars,   only: DtbFormat
  !
  integer,intent(in):: f
  !
  character(len=512):: L,W
  logical :: EoL,bGas,bMin
  integer :: ios,I
  type(T_DtbMinThr) :: M
  !
  if(idebug>1) write(fTrc,"(/,A)") "< DtbMinThr_ToLine"
  !
  write(fTrc,'(A)') &
  & "stoikio in "//trim(DirDtbLog)//"min_thr_line.log"
  if(fLin==0) then
    call GetUnit(fLin)
    open(fLin,file=trim(DirDtbLog)//"min_thr_line.log")
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"min_thr_line.log",&
    & "multiline datafile from therak converted to monoline")
  end if
  !
  FilCode=trim(DtbFormat) !-> default value
  M%Name="Z"
  !
  DoFile: do

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!")   cycle DoFile
    call AppendToEnd(L,W,EoL)

    select case(W)
    !
    case("CODE")
      call LinToWrd(L,W,EoL)
      FilCode=trim(W)
      cycle DoFile

    case("ENDINPUT")
      exit DoFile

    case("ENDSPECIES")
      exit DoFile

    case("MINERAL","GAS")

      select case(W)
        case("MINERAL"); bGas=.false.; bMin=.true.
        case("GAS");     bGas=.true.;  bMin=.false.
      end select

      DoReadMin: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=="!") cycle DoReadMin !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        select case(W)
          case("ENDINPUT"); exit DoFile
          case("ENDSPECIES"); exit DoFile
          case("END","ENDMINERAL","ENDGAS"); cycle DoFile !exit DoReadMin
        end select

        if(W(1:1)/="&") then
        !if first char is not "&", the line may contain name, formula, etc
          write(fLin,*)
          M%Name=trim(W)
          if(bGas) write(fLin,'(A,A1)',advance="NO") "GAS",T_
          if(bMin) write(fLin,'(A,A1)',advance="NO") "MIN",T_
          call LinToWrd(L,W,EoL); M%Formula=trim(W)
          call LinToWrd(L,W,EoL) !; M%Abbr=trim(W)
          write(fLin,'(3(A,A1))',advance="NO") &
          & trim(FilCode),T_,M%Name,T_,M%Formula,T_
          !
          cycle DoReadMin
        end if

        if(W(1:1)=="&") then
          ! if first char is "&", 
          ! then the line is the continuation of preceding, i.e. contains data
          do I=1,7
            call LinToWrd(L,W,EoL)
            if(.not. EoL) then
              write(fLin,'(A,A1)',advance="NO") trim(W),T_
            else
              write(fLin,'(A,A1)',advance="NO") "_",T_
            end if
          end do
        end if !(W(1:1)=="&")

      end do DoReadMin
    !endcase("MINERAL")
    end select
  end do DoFile
  !
  write(fLin,*)
  !
  if(idebug>1) write(fTrc,"(A,/)") "</ DtbMinThr_ToLine"
  !
  return
end subroutine DtbMinThr_ToLine

end module M_Dtb_Read_DtbMinThr
