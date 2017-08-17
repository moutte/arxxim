module M_Test_Read
!==================================================================
! Read Block Test and Options
!==================================================================
  use M_Trace,only: iDebug,fTrc,Stop_,T_
  use M_Kinds
  
  implicit none
  
  private
  
  integer :: OptComputeSize
  character(len=10), allocatable :: OptCompute(:)
  logical :: Ok_Initialized = .false.
  
  public:: Test_Read_Init
  public:: Test_Read_OptCompute
  
contains

  function Test_Read_OptCompute(I, Value) result (Ok)
  !================================================== 
  ! Get the COMPUTE Field number i from block TEST
  !==================================================
    implicit none
    logical :: Ok
    integer, intent(in) :: I
    character(len=*), intent(out) :: Value
    !---
    if (I>OptComputeSize) then
      Ok = .false.
      Value = "NONE"
    else
      Ok = .true.
      Value = OptCompute(I)
    end if
    
  end function Test_Read_OptCompute

  !---
  
  subroutine Test_Read_Init(Ok)
  !================================================== 
  ! Init Test Block => Reinit if Necessary
  !==================================================
    implicit none
    logical,intent(out) :: Ok
    !--
    if (.not. Ok_Initialized) call Test_Read_ReInit(Ok)
  
  end subroutine Test_Read_Init
  
  !---

  subroutine Test_Read_ReInit(Ok)
  !================================================== 
  ! Reinit Test block = Read Test Block
  !==================================================
    use M_Files_Vars,only:NamFInn
    use M_IOTools 
    implicit none
    !-----
    logical,intent(out) :: Ok
    !-----
    logical :: sEOL,EoL
    integer :: nCompute
    character(len=512):: L,W
    integer :: ios,F
    integer :: i
    !-----
    Ok = .false.
    nCompute =0

    if(idebug>1) write(fTrc,'(/,A,/)') "< Test_ReadInput"

    !// Estimate the maximal number of values to store
    call GetUnit(f)
    open(f,file=trim(NamFInn),STATUS='OLD')

    nCompute = 0
    do 
      read(f,'(A)',iostat=ios) L
      if(ios/=0) exit
      call LinToWrd(L,W,EoL)
      if(W(1:1)=='!') cycle
      if(trim(W)=="ENDINPUT") exit
      if(trim(W)=="COMPUTE") nCompute = nCompute +1
    end do
    close(f)
    
    !! print *,'nCompute=',nCompute
    !! print *,'..DEBBUGG'  ;  pause
    
    !// Allocate the table
    if(allocated(OptCompute)) deallocate(OptCompute)
    allocate(OptCompute(10))

    !// Read the File and store the values
    call GetUnit(f)
    open(f,file=trim(NamFInn),STATUS='OLD')
    !
    nCompute = 0
    !
    DoFile: do
      !
      read(f,'(A)',iostat=ios) L
      if(ios/=0) exit DoFile
      call LinToWrd(L,W,EoL)
      if(W(1:1)=='!') cycle DoFile
      call AppendToEnd(L,W,EoL)
      
      select case(trim(W))
          
      case("ENDINPUT") ; exit DoFile
      
      case("TEST")
        Ok= .true.
        !
        DoTest: do
          !
          read(F,'(A)',iostat=ios) L
          if(ios/=0) exit DoFile
          call LinToWrd(L,W,EoL)
          if(W(1:1)=='!') cycle DoTest
          call AppendToEnd(L,W,EoL)
          
          select case(trim(W))
          
          case("ENDINPUT")
            exit DoFile
          case("END","ENDTEST")
            exit DoTest
          
          case("COMPUTE")
            call LinToWrd(L,W,sEol) 
            nCompute = nCompute + 1
            OptCompute(nCompute)= trim(W)
          
          case default
            Ok= .false.
            print *,trim(W)//" =INVALID keyword in TEST block"
          
          end select
          !
        end do Dotest
        !
      end select
      !
    end do DoFile
    close(f)
    !
    OptComputeSize = nCompute
    Ok_Initialized = .true.
    
    do nCompute=1,OptComputeSize
      write(fTrc,'(2X,A)') OptCompute(nCompute)
    end do
    !
    if(idebug>1) write(fTrc,'(/,A,/)') "</ Test_ReadInput"

  end subroutine Test_Read_ReInit

end module M_Test_Read
