module M_Dtb_Read_DtbAquHkf
  use M_Kinds
  use M_IoTools
  use M_Trace
  use M_Dtb_Read_Vars
  !
  implicit none
  !
  private
  !
  public:: DtbAquHKF_Read
  public:: DtbAquHKF_Build
  !
contains

subroutine DtbAquHKF_Build
!--
!-- transfer linked list to array of T_DtbAquHkf
!--
  use M_Dtb_Vars,only: vDtbAquHkf
  !
  type(T_LisAquHkf),pointer:: LisAquCur, LisAquPrev
  integer:: I
  !
  !!if(iDebug>2) print '(A)',"< DtbAquHKF_Build"
  !
  if(allocated(vDtbAquHkf)) deallocate(vDtbAquHkf)
  allocate(vDtbAquHkf(nAquHkf))
  !
  I=0
  LisAquCur=>LisAquHkf
  do while (associateD(LisAquCur))
    I=I+1
    vDtbAquHkf(I)=LisAquCur%Value
    if(iDebug>5) print '(2A)', vDtbAquHkf(I)%Name,trim(vDtbAquHkf(I)%Num)
    LisAquPrev=>LisAquCur
    LisAquCur=> LisAquCur%next
    deallocate(LisAquPrev)
  end do
  !
  !!if(iDebug>2) print '(A)',"</ DtbAquHKF_Build"
  !
  if(iDebug>5) call Pause_
  
  return
end subroutine DtbAquHKF_Build

subroutine DtbAquHKF_Read(F,vEle,N)
!--
!-- reads data directly from MODIFIED SLOP98.DAT
!-- (the newest version of SPRONS.DAT,the text-file for SupCrt)
!-- ("normally", SUPCRT92.FOR reads from DPRONS,
!-- a "direct-access file" that is built from SLOP98.DAT using CPRONS92.FOR)
!--
  use M_Dtb_Const,  only: T_CK,Tref,Pref,S0_Hydrogen
  use M_Files,      only: DirDtbLog,Files_Index_Write
  use M_T_Element,  only: T_Element,Formula_Build,Formula_Read
  use M_Dtb_Read_Tools
  !
  use M_Dtb_Vars,   only: DtbFormat
  
  !---------------------------------------------------------------------
  integer,        intent(in)   :: F !input file
  type(T_Element),intent(in)   :: vEle(:)
  integer,        intent(inout):: N
  !---------------------------------------------------------------------
  type(T_LisAquHkf),pointer,save:: LisCur
  type(T_DtbAquHkf):: M
  !
  character(len=255):: L,W,sFormul
  !
  character(len=4)  :: ICode
  logical :: EoL,fOk,EcformIsOk
  integer :: ios,K,ZSp,Div_,iEl
  integer :: I,FF
  real(dp):: X
  real(dp),dimension(dimV)::vX
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
  integer,allocatable::vStoik(:) !for Formula_Read
  !
  !--- for formula translation --
  character(len=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  character(len=2),allocatable:: vElement(:)  !
  integer:: fFormula
  !---/
  !---------------------------------------------------------------------
  ! if (N>1) LisCur => LisAquHkf
  !
  if(idebug>1) write(fTrc,"(/,A)") "< DtbReadAquHKF"
  !
  !--------------------------------------------------------------- trace
  FF= 0
  if(iDebug>2) then
    write(fTrc,'(A)') &
    & "stoikio in "//trim(DirDtbLog)//"aqu_hkf_stoik.log"
    !
    call Files_Index_Write(fHtm,&
    & trim(DirDtbLog)//"aquhkf_stoik.log",&
    & "check species stoikio of aqu.species from hkf-type database")
    !
    call GetUnit(FF)
    open(FF,file=trim(DirDtbLog)//"aquhkf_stoik.log")
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
  !--- for formula translation --
  CodFormula= "ECFORM"
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  allocate(vStoik(1:size(vEle)))
  !
  FilCode= trim(DtbFormat)
  !
  !--------------------------------- initialize vStrField(:), vifield(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "ABBREV      ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  if (FilCode(1:5)=="OBIGT") then
    L= "INDEX NAME SKIP SCFORM TYPE SKIP SKIP SKIP PARAMETERS"
  else ! case of SLOP98.DAT
    L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
  end if
  !
  if(iDebug==4) print *,"< Default values"
  call FieldList_Read(L,vStrField,vifield)
  if(iDebug==4) print *,"</ Default values"
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
    !
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=="!") cycle DoFile
    call AppendToEnd(L,W,EoL)
    !
    !---------------------------------------- read header, if present --
    !------ if the line begins with any member of array vStrField(:), --
    !------------------------------ then it contains the line headers --
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
      if(iDebug==4) print *,"< values from file"
      call FieldList_Read(L,vStrField,vifield)
      if(iDebug==4) print *,"</ values from file"

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
    !------------------------------------------------ for old formats --
    case("CODE")
      !-> can change the CODE, i.e. the data source inside the block
      call LinToWrd(L,W,EoL)
      FilCode=trim(W)
      !
      ! if FilCode changes, must re-initialize vifield ...!
      if (FilCode(1:5)=="OBIGT") then
        L= "INDEX NAME SKIP SCFORM type SKIP SKIP SKIP PARAMETERS"
      else ! case of SLOP98.DAT
        L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
      end if
      call FieldList_Read(L,vStrField,vifield)
      !
      cycle DoFile

    case("FORMULA")
      call LinToWrd(L,W,EoL)
      CodFormula= trim(W)
      !
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
      !
      cycle DoFile !-----------------------------------------------cycle
    !--------------------------------------------------/ for old formats
    case default
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= trim(W)//" "//trim(L)
      !
    end select
    !---------------------------------------/ process first word of line
    !
    M%Num=     "0"
    M%Name=    "0"
    M%Abbr=    "0"
    M%Formula= "0"
    !-------------------------------------------------- scan a data line
    do I= 1,vifield(10)-1
    ! read words up to position before the parameters
      !
      call LinToWrd(L,W,EoL,"NO")
      !
      ! (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !   "ABBREV      ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
      !
      if(I==vifield(2)) then  ! INDEX
        call Str_Upper(W)  ;  M%Num= trim(W)
      end if

      if(I==vifield(3)) then  ! NAME
        call Str_Upper(W)  ;  M%Name= trim(W)
      end if

      if(I==vifield(4)) then ! SCFORM
        !
        call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        if(.not.EcformIsOk) cycle DoFile !-------------------------cycle
        !
        call Str_Upper(W)  ;  M%Formula=trim(W)
      end if

      if(I==vifield(5)) then ! ECFORM
        call Str_Upper(W)  ;  M%Formula=  trim(W)
      end if

      !if(I==vifield(9)) then
      !  call Str_Upper(W)  ;  CodFitting= trim(W)
      !end if

      !! if(I==vifield(6)) call WrdToReal(W,X) ! size
      !! if(I==vifield(7)) call WrdToReal(W,X) ! VOLUME
      !! if(I==vifield(8)) call WrdToReal(W,X) ! DENSITY

    end do
    !
    call ReadRValsV(L,K,vX)
    !-------------------------------------------------/ scan a data line
    !!
    !! select case(W)
    !!   case("ENDINPUT"); exit DoFile !DoFile
    !!   case("END","ENDSPECIES"); exit DoFile !DoFile
    !!   !!case("END"); exit DoFile
    !!   case("FORMULA")
    !!     call LinToWrd(L,W,EoL)
    !!     CodFormula= trim(W)
    !!     if(iDebug>2 .and. CodFormula=="SCFORM") then
    !!       call GetUnit(fFormula)
    !!       open(fFormula,file="debug_formula.log")
    !!       write(fFormula,'(A,/)') "resuts of formula conversion"
    !!     end if
    !!     cycle DoFile
    !! end select
    !! !
    !! M%Num=trim(W)
    !! !
    !! if (FilCode(1:5)=="OBIGT") then
    !!   CodFormula="SCFORM"
    !!   !
    !!   call LinToWrd(L,W,EoL); M%Name=trim(W) !if(idebug>1) write(fTrc,"(A)") M%Name
    !!   call LinToWrd(L,W,EoL) !skip ABBRV
    !!   call LinToWrd(L,W,EoL,"NO") !-> compact formula, character case is conserved :!!
    !!   !
    !!   if(CodFormula=="SCFORM") then
    !!     call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
    !!     if(.not.EcformIsOk) cycle DoFile
    !!   end if
    !!   !
    !!   call Str_Upper(W)
    !!   !
    !!   M%Formula=trim(W)
    !!   !
    !!   call LinToWrd(L,W,EoL) !skip STATE
    !!   call LinToWrd(L,W,EoL) !skip SOURCE1
    !!   call LinToWrd(L,W,EoL) !skip SOURCE2
    !!   call LinToWrd(L,W,EoL) !DATE
    !!   !
    !!   call ReadRValsV(L,K,vX)
    !!   !
    !! else ! default FilCode is SLOP98
    !!   !
    !!   call LinToWrd(L,W,EoL) !skip ABREV
    !!   call LinToWrd(L,W,EoL) !skip NAME
    !!   !
    !!   call LinToWrd(L,W,EoL); M%Name=trim(W) !if(idebug>1) write(fTrc,"(A)") M%Name
    !!   call LinToWrd(L,W,EoL); M%Formula=trim(W) !if(idebug>1) write(fTrc,"(A)") M%Formula
    !!   call LinToWrd(L,W,EoL) !DATE
    !!   call LinToWrd(L,W,EoL) !REF
    !!   !
    !!   call ReadRValsV(L,K,vX)
    !!   !
    !! end if
    !
    !------------------------------------------------------ read formula
    call Formula_Read(M%Formula,vEle,ZSp,Div_,fOk,vStoik)
    if(.not. fOk) cycle DoFile !-----------------------------------cycle
    !-----------------------------------------------------/ read formula
    !
    if (FilCode(1:5)=="OBIGT") then
      ! in OBIGT database, Cp and V are tabulated fot aqu'species
      ! but not used as input parameters -> vX(4:5) ignored
      !
      M%G0R=vX(1);   M%H0R=vX(2);  M%S0_= vX(3)
      M%A1= vX(6);   M%A2= vX(7);  M%A3=  vX(8); M%A4=vX(9);
      M%C1= vX(10);  M%C2= vX(11)
      M%wref=vX(12)
      M%Chg= INT(vX(13))
      !
    else
      !
      M%G0R=vX(1);  M%H0R=vX(2); M%S0_= vX(3)
      M%A1= vX(4);  M%A2= vX(5); M%A3=  vX(6); M%A4=vX(7);
      M%C1= vX(8);  M%C2= vX(9)
      M%wref= vX(10)
      M%Chg=  INT(vX(11))
      !
    end if
    !
    M%A1=    M%A1*1.0D-1
    M%A2=    M%A2*1.0D02
    M%A4=    M%A4*1.0D04
    M%C2=    M%C2*1.0D04
    M%Wref=  M%wref*1.0D05
    !
    N=N+1
    call IntToStr4(N,ICode)
    M%Num= trim(FilCode)//"_"//trim(ICode)

    M%S0Ele= dot_product(vStoik(:),vEle(:)%S0) - ZSp *S0_Hydrogen ! is in Joule
    M%S0Ele= M%S0Ele /real(Div_)

    M%WeitKg= dot_product(vStoik(:), vEle(:)%WeitKg) /real(Div_) !-> in Kg
    !
    !-------------------------------------------- "fill" the linked list
    if(N==1) then
      allocate(LisAquHkf)
      nullify(LisAquHkf%next)
      LisAquHkf%Value= M
      LisCur=> LisAquHkf
    else
      allocate(LisCur%next)
      nullify(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    end if
    !-------------------------------------------/ "fill" the linked list
    !
    if(FF>0) then
      write(FF,"(A15,A1)",advance="NO") M%Name,T_
      do iEl=1,size(vStoik)
        write(FF,"(I3,A1)",advance="NO") vStoik(iEl),T_
      end do
      X= M%H0R -Tref*M%S0_
      write(FF,'(3(G15.6,A1))',advance="NO") &
      & M%G0R,T_, X,T_, (M%G0R-X)/Tref,T_
      !
      ! build new formula, with fixed element order
      call Formula_Build(vEle,vStoik,Zsp,Div_,sFormul)
      !!!build new formula (without coeff"s !!), with fixed element order
      !!sFormul=""
      !!do iEl=1,size(vEle)
      !!  if(vStoik(iEl)>0) sFormul=trim(sFormul)//trim(vEle(iEl)%NamEl)//""
      !!end do
      write(FF,"(I3,A1,3(A,A1))") &
      & M%Div,T_,trim(M%Formula),T_,trim(M%Name),T_,trim(sFormul)
      !
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
  if(idebug>1) write(fTrc,"(A,/)") "</ DtbReadAquHKF"
  !
  return
end subroutine DtbAquHKF_Read

end module M_Dtb_Read_DtbAquHkf
