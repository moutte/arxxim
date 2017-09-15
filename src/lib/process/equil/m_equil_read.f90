module M_Equil_Read
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  use M_Kinds
  implicit none
  !
  private
  !
  public:: Equil_Read_YesList
  public:: Equil_Read_Debug
  public:: Equil_Read_Numeric
  public:: Equil_Read_PhaseAdd
  
contains

subroutine Equil_Read_PhaseAdd(vFas)
!--
!-- read block SYSTEM.ROCK,
!-- -> read the non-fluid part of the system,
!-- -> update vFas%Mole
!--
  use M_IOTools
  use M_Files,  only: NamFInn
  use M_T_Phase,only: T_Phase,Phase_Index
  !
  type(T_Phase), intent(inout) :: vFas(:)
  !
  character(len=512):: L,W,W1
  logical :: EoL
  integer :: F,ios,I
  real(dp):: X
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_PhaseAdd"
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  !
  DoFile: do
    ! 
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
      !
      case("ENDINPUT"); exit DoFile
      !
      case("SYSTEM.ROCK") !
        DoBlock: do
          !
          read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=='!') cycle DoBlock !skip comment lines
          call AppendToEnd(L,W,EoL)
          !
          select case(W)
          case("ENDINPUT")              ;  exit DoFile
          case("END","ENDSYSTEM.ROCK")  ;  exit DoBlock
          case("MOLE","GRAM")           ;  call LinToWrd(L,W1,EoL)
          case default
            call Stop_(trim(W)//" <- unknown keyword in SYSTEM.ROCK !!!") 
          end select
          !
          I= Phase_Index(W1,vFas)
          !
          if(I==0) &
          & call Stop_(trim(W1)//" <-phase unknown (SYSTEM.ROCK)")
          !
          if(vFas(I)%iSpc==0) &
          & call Stop_(trim(W1)//" <-not valid, use only pure non'aqu'phases (SYSTEM.ROCK)")
          !
          call LinToWrd(L,W1,EoL)
          call WrdToReal(W1,X)
          !
          if(X<=Zero) call Stop_(trim(W1)//" <-not valid as species amount  (SYSTEM.ROCK)")
          if(W=="GRAM") X= X /vFas(I)%WeitKg /1.0D3
          !
          vFas(I)%MolFs= X
          !
        end do DoBlock
      !endcase("SYSTEM.ROCK")
      !
    end select
    !
  end do DoFile
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_PhaseAdd"
  !
end subroutine Equil_Read_PhaseAdd

subroutine Equil_Read_YesList(nFas,vFas,vYesList)
!--
!-- reads lists of included or excluded phases -> update vYesList
!--
  use M_T_Phase,only: T_Phase
  !
  integer,      intent(in) :: nFas
  type(T_Phase),intent(in) :: vFas(:)
  logical,      intent(out):: vYesList(:)
  !
  logical:: vInclud(nFas),vExclud(nFas)
  !vInclude= included phases; vExclude= excluded phases
  integer:: I
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_YesList"
  !
  !-------------------------------- read list of "equilibrium phases" --
  vInclud=.false.
  vExclud=.false.
  !
  call Equil_Read_EquilPhase(vFas,vInclud,vExclud) !,InOk,ExOk)
  !
  vYesList=.true. !-> default= all phases from database taken into account
  if(count(vInclud)>0) vYesList(:)= vInclud(:) !-> only included phases
  if(count(vExclud)>0) vYesList(:)= vYesList(:) .and. (.not. vExclud(:))
  !
  !-------------------------------/ read list of "equilibrium phases" --
  !
  !! vYesList(1)=.false. !-> ??? to exclude species H2O as phase ???
  !
  !------------------------------------------------------------ trace --
  if(iDebug==4) then
    print '(/,A)',"< -- List of Phases --"
    do I=1,nFas
      if(vYesList(I)) print '(A)',vFas(I)%NamFs
    end do
    print '(A,/)',"</-- List of Phases --"
  end if
  !
  if(idebug>1) then
    do I=1,nFas
      if(vYesList(I)) write(fTrc,'(A)') vFas(I)%NamFs
    end do
  end if
  !-----------------------------------------------------------/ trace --
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_YesList"
end subroutine Equil_Read_YesList
  !  
subroutine Equil_Read_Debug(DebFormula,DebNewt,DebJacob)
!--
!-- read block CONDITIONS from arxim.inn to retrive DEBUG options
!--
  use M_Trace,only: LWarning
  use M_IoTools
  use M_Files,only: NamFInn
  !
  logical,intent(out):: DebFormula,DebNewt,DebJacob
  !
  character(len=255):: L
  character(len=80) :: W
  logical           :: EoL,Ok
  integer           :: f,N,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_Debug"
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  Ok=.false.
  !
  DoFile: do 
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    
    select case(trim(W))
    !
    case("ENDINPUT"); exit  DoFile 
    !
    case("CONDITIONS") !logical global variables
      Ok=.true.
      !
      DoDebug: do
        
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoDebug
        call AppendToEnd(L,W,EoL)
        
        select case(trim(W))
          !
          case("ENDINPUT"); exit DoFile
          
          case("END","ENDCONDITIONS"); exit DoDebug
          
          case("DEBUG")
            call LinToWrd(L,W,EoL); call WrdToInt(W,N)
            if(N>0 .and. iDebug==0) iDebug=N
            !iDebug is 0 if it has not been already assigned in command line
          !
        end select
      
      end do DoDebug
      !
    end select
  end do DoFile
  close(f)
  !
  LWarning= (idebug>1)
  !
  DebFormula= (iDebug==4)
  DebJacob=   (iDebug==4)
  DebNewt=    (iDebug==4 .or. iDebug==9)
  !
  if(.not.Ok) then 
    if(idebug>1) write(fTrc,'(A)') &
    & "Block CONDITIONS not Found -> using default values for DEBUG"
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_Debug"
end subroutine Equil_Read_Debug

subroutine Equil_Read_OutputOptions(OutDistrib)
!--
!-- provisional, may be used for other output options ...
!-- read block CONDITIONS from arxim.inn for OUTPUT options
!--
  use M_IoTools
  use M_Files,only: NamFInn
  !
  logical,intent(out):: OutDistrib
  !
  character(len=255):: L
  character(len=80) :: W
  logical           :: EoL,Ok
  integer           :: f,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_OutputOptions"
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  Ok=.false.
  !
  OutDistrib= .false.
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(trim(W))
      !
      case("ENDINPUT")
        exit  DoFile 
      !
      case("CONDITIONS")
        Ok=.true.
        !
        DoDebug: do
          read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=='!') cycle DoDebug
          call AppendToEnd(L,W,EoL)
          select case(trim(W))
            case("ENDINPUT"); exit DoFile
            case("END","ENDCONDITIONS"); exit DoDebug
            !
            case("OUTPUT.DISTRIB"); OutDistrib=.true. 
            !
          end select
        end do DoDebug
        !
    end select
    !
  end do DoFile
  !
  close(f)
  !
  if(.not.Ok) then 
    if(idebug>1) write(fTrc,'(A)') &
    & "Block CONDITIONS not Found -> using default values for OUTPUT options"
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_OutputOptions"
  !
end subroutine Equil_Read_OutputOptions

subroutine Equil_Read_EquilPhase(vFas,vInclude,vExclude) !,InOk,ExOk)
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_Files,only: NamFInn
  use M_T_Phase,only: T_Phase,Phase_Index
  !
  type(T_Phase),intent(in) :: vFas(:)
  logical,      intent(out):: vInclude(:)
  logical,      intent(out):: vExclude(:)
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,I,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_EquilPhase"
  !
  vInclude=.false. !InOk=.false.; 
  vExclude=.false. !ExOk=.false.; 
  !
  !!vInclude(:)= vFas(:)%iSpc/=0
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  DoFile: do
  
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    select case(trim(W))
      
      case("ENDINPUT"); exit DoFile
      
      case("EQUIL.INCLUDE","EQU.INCLUDE")
        
        DoInclude: do
          
          read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=='!') cycle DoInclude !skip comment lines
          call AppendToEnd(L,W,EoL)
          select case(trim(W))
            case("ENDINPUT"); exit DoFile
            case("END","ENDEQU.INCLUDE","ENDEQUIL.INCLUDE"); exit DoInclude
          end select
          
          I= Phase_Index(W,vFas)
          
          if(I>0) then
            if(vFas(I)%iSpc/=0) then !restriction to pure phases, provisionally
              vInclude(I)=.true.
              if(iDebug>2) print      '(A,A1,G15.6)', vFas(I)%NamFs,T_,vFas(I)%Grt
              if(idebug>1) write(fTrc,'(A,A1,G15.6)') vFas(I)%NamFs,T_,vFas(I)%Grt
            else
              !!vInclude(I)=.true.
              !!if(iDebug>2) print      '(A,A1,G15.6)', vFas(I)%NamFs,T_,vFas(I)%Grt
              !!if(idebug>1) write(fTrc,'(A,A1,G15.6)') vFas(I)%NamFs,T_,vFas(I)%Grt
              if(iDebug>2) print      '(2A)', trim(W)," is not pure -> not included"
              if(idebug>1) write(fTrc,'(2A)') trim(W)," is not pure -> not included"
            end if
          else
            if(iDebug>2) print '(3A)',"WARNING, ",trim(W)," is not in current database"
            if(idebug>1) write(fTrc,'(3A)') "WARNING !!! ",trim(W)," is not in current database"
          end if
        
        end do DoInclude
        
        !!InOk= count(vInclude)>0
      !endcase
      
      case("EQUIL.EXCLUDE","EQU.EXCLUDE")
        
        DoExclude: do
          
          read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=='!') cycle DoExclude !skip comment lines
          call AppendToEnd(L,W,EoL)
          select case(trim(W))
            case("ENDINPUT"); exit DoFile
            case("END","ENDEQUIL.EXCLUDE","ENDEQU.EXCLUDE"); exit DoExclude
          end select
          
          I= Phase_Index(W,vFas)
          if(I>0) vExclude(I)=.true.
          
          if(I>0) then !--====================================< trace ==
            if(idebug>1) &
            & write(fTrc,'(2A,A1,G15.6)') "EXCLUDED, Grt=",vFas(I)%NamFs,T_,vFas(I)%Grt
            if(iDebug>1) &
            & print '(2A,1X,G15.6)',"EXCLUDED, Grt=",vFas(I)%NamFs,vFas(I)%Grt
          end if !--==========================================</ trace ==
          
        end do DoExclude
        
        !!ExOk= count(vExclude)>0
      !endcase
    end select
    
  end do DoFile
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_EquilPhase"
  !
end subroutine Equil_Read_EquilPhase

subroutine Equil_Read_Numeric( &
& initMolNu, &
& LogForAqu, &
& DirectSub, &
& cMethod,   &
& bFinDif,   &
& Error,     &
& ErrorMsg   &
& )
!--
!-- read block EQU.NUMERIC
!--
  use M_IoTools
  use M_Files,only: NamFInn
  !
  real(dp),    intent(out)  :: initMolNu
  logical,     intent(inout):: LogForAqu,DirectSub
  character(*),intent(out)  :: cMethod
  logical,     intent(inout):: bFinDif
  logical,     intent(out)  :: Error
  character(*),intent(out)  :: ErrorMsg
  !
  character(len=255):: L
  character(len=80) :: W
  logical           :: EoL,Ok
  integer           :: f,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Read_Numeric"
  !
  Error= .false.
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  Ok=.false.
  
  DoFile: do 
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit DoFile 
    !
    case("EQU.NUMERIC")
      Ok=.true.
      DoCond: do
      
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoCond
        call AppendToEnd(L,W,EoL)
        
        select case(W)
        
        case("ENDINPUT")                ;  exit DoFile
        case("END","ENDEQU.NUMERIC")    ;  exit DoCond
        
        case("INITIAL")
          call LinToWrd(L,W,EoL)
          call WrdToReal(W,initMolNu)
        
        case("LOGFORAQU")
          call LinToWrd(L,W,EoL)
          select case(trim(W))
          case("TRUE", "T", "YES")  ;  LogForAqu= .true.
          case("FALSE", "F", "NO")  ;  LogForAqu= .false.
          case default
            Error= .true.
            ErrorMsg= trim(W)//" is invalid keyword in EQU.NUMERIC"
          end select
        
        !!! obsolete : replace by LOGFORAQU= YES/NO !!
        !! case("AQUEOUS")
        !!   call LinToWrd(L,W,EoL)
        !!   select case(trim(W))
        !!   case("MOLE")  ;  LogForAqu= .false.
        !!   case("LOG")   ;  LogForAqu= .true.
        !!   case default
        !!   Error= .true.
        !!   ErrorMsg= trim(W)//" is invalid keyword in EQU.NUMERIC"
        !!   end select
          
        case("DIRECT")
          call LinToWrd(L,W,EoL)
          select case(trim(W))
          case("TRUE", "T", "YES")  ;  DirectSub= .true.
          case("FALSE", "F", "NO")  ;  DirectSub= .false.
          case default
            Error= .true.
            ErrorMsg= trim(W)//" is invalid keyword for DIRECT in EQU.NUMERIC"
          end select
        
        case("JACOBIAN")
          call LinToWrd(L,W,EoL)
          select case(trim(W))
          case("NUMERIC")   ;  bFinDif= .true.
          case("ANALYTIC")  ;  bFinDif= .false.
          case default
            Error= .true.
            ErrorMsg= trim(W)//" is invalid keyword in EQU.NUMERIC"
          end select
        
        case("METHOD")
          call LinToWrd(L,W,EoL)
          select case(trim(W))
          case("NEWTON")       ;   cMethod= trim(W)
          case("NEWTLNSRCH")   ;   cMethod= trim(W)
          case("NEWTONPRESS")  ;   cMethod= trim(W)
          case("BROYDEN")      ;   cMethod= trim(W)
          case("NEWTONCHESS")  ;   cMethod= trim(W)
          case("TENSOLVE_1")   ;   cMethod= trim(W)
          case("NEWTONKELLEY") ;   cMethod= trim(W)
          case("NEWTONWALKER") ;   cMethod= trim(W)
          case default
            Error= .true.
            ErrorMsg= trim(W)//" is invalid METHOD in EQU.NUMERIC"
          end select
        !_
        case default
          Error= .true.
          ErrorMsg= trim(W)//" is invalid keyword in EQU.NUMERIC"
          !print *,trim(W)//" = unknown keyword in EQU.NUMERIC !!"
        
        end select
      
      end do DoCond
    !endcase("CONDITIONS") 
    end select
    
  end do DoFile
  close(f)
  !
  if(.not.Ok) then 
    if(idebug>1)  write(fTrc,'(A)') "Block EQU.NUMERIC not Found ...!!!"
    !if(idebug>1)  print '(A)',"WARNING: Block CONDITIONS not Found, using default values"
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Read_Numeric"
  !
end subroutine Equil_Read_Numeric

end module M_Equil_Read

