module M_T_EosIfp_Param
  
  use M_Kinds,ONLy: dp
  implicit none

  !// Substances specific parameters
  type:: T_EosIfp_Param
    integer :: NCONS
    character(len=30), pointer :: NAME(:)
    real(dp), pointer :: TC(:), PC(:) ! (NCONS) !critical temperature and pressure
    real(dp), pointer :: OMEGA(:)     ! (NCONS) ! acentric factor
    real(dp), pointer :: COINT(:,:)   ! (NCONS, NCONS) ! interaction coefficients
    character(len=4) :: MODEL             ! type OF MODEL 
    logical :: LOPT_OMGA, LOPT_OMGB       ! ACTIVATION OF optional parameterS    
    real(dp), pointer :: OMGA(:), OMGB(:)      ! (NCONS) ! optional parameterS
    real(dp) :: DLTA1, DLTA2
  end type T_EosIfp_Param

  !// Public functions
  ! public :: EosIfp_Param_delete
  public :: EosIfp_Param_CreateFromFile

contains

subroutine EosIfp_Param_new( T, NCONS)
  type(T_EosIfp_Param), intent(inout) :: T
  integer, intent(in) :: NCONS
  integer :: I
  !-------------

  ! sizes
  T% NCONS = NCONS

  ! arrays
  allocate( T% NAME(NCONS) )
  allocate( T% TC(NCONS) )
  allocate( T% PC(NCONS) )
  allocate( T% OMEGA(NCONS) )
  allocate( T% COINT(NCONS,NCONS) )
  allocate( T% OMGA(NCONS) )
  allocate( T% OMGB(NCONS) )
  !---- init
  T%LOPT_OMGB = .false.
  T%LOPT_OMGA = .false.
  T%OMGA= 0.D0
  T%OMGb= 0.D0

  !----------------
  do I=1, NCONS
    T% NAME(I) = ""
  end do

  T% TC   = 0.D0
  T% PC   = 0.D0
  T% OMEGA   = 0.D0
  T% COINT   = 0.D0

end subroutine EosIfp_Param_new

subroutine EosIfp_Param_Print(T, funit)

  type(T_EosIfp_Param), intent(inout) :: T
  integer, intent(in) :: funit
  integer :: i,j, NCONS

  !-----
  NCONS = T%NCONS
  write(funit,*) ">> Print Param "

  do I=1, NCONS
    write(funit,*) "------------------------------"
    write(funit,*) "NAME  = ", T% Name(i)
    write(funit,*) "TC    = ", T%TC(i)
    write(funit,*) "PC    = ", T%PC(i)
    write(funit,*) "OMEGA = ", T%OMEGA(i)
    write(funit,*) "------------------------------"
    write(funit,*) "OMGA  = ", T%OMGA(i)
    write(funit,*) "OMGB  = ", T%OMGB(i)
  end do

  write(funit,*) "-----------------------------------------"
  write(funit,*) "DLTA1 = ", T%DLTA1
  write(funit,*) "DLTA2 = ", T%DLTA2

  write(funit,*) "-----------------------------------------"
  write(funit,*) "BINARY INTERACTION COEFFICIENTS MATRIX : "

  do J=1, NCONS
    write(funit,*) (T%COINT(i,j), i=1, NCONS)
  end do
  write(funit,*) "-----------------------------------------"

end subroutine EosIfp_Param_Print


!// reading parameters 

subroutine EosIfp_Param_CreateFromFile( T, T_Const, FileInput, NCONS )

  use M_T_EosIfp_Const
  
  character(len=*), intent(in) :: FileInput
  type(T_EosIfp_Param), intent(inout) :: T
  type(T_EosIfp_Const), intent(in) :: T_Const
  !-----
  integer :: NCONS_FILE
  integer :: NCONS 
  integer :: funit = 103

  real(dp) :: OMGA0, OMGB0
  real(dp) :: OMGA(NCONS)
  real(dp) :: OMGB(NCONS)
  integer :: i, j
  real(dp) :: temp

  !-------- Reading Compositional system and Parameters

  open ( UNIT= FUNIT, FILE = FILEINPUT)

  !// get size
  read( funit,'(I8)') NCONS_FILE

  call EosIfp_Param_new(T, NCONS)

  !// get critical parameters 
  do I=1, NCONS
    read(funit,*) T%NAME(i)
    read(funit,*) T%TC(i)
    read(funit,*) T%PC(i)
    read(funit,*) T%OMEGA(i)
  end do

  !// get interaction coefficient matrix
  read(funit,*)
  do i=1, NCONS
    do j=1, NCONS
      read(funit,*) T%COINT(i, j)
    end do
  end do

  !// get model and standard parameters for this model
  read(funit,*) T%MODEL
  select case( T%MODEL)
  case ("PR1", "PR2")
    OMGA0   = T_Const %OMGAPR
    OMGB0   = T_Const %OMGBPR
    T%DLTA1 = T_Const %DLPR1
    T%DLTA2 = T_Const %DLPR2
  case ("SRK1", "SRK2") 
    OMGA0   = T_Const%OMGARK
    OMGB0   = T_Const%OMGBRK
    T%DLTA1 = T_Const %DLRK1
    T%DLTA2 = T_Const %DLRK2
  case default
    write(*,*) "Bad Model : <", T%MODEL, ">"
    stop "Fatal Error"
  end select

  !// set optional parameters OMGA, OMGB
  if ( T% LOPT_OMGA ) then
    do I=1, NCONS
      T%OMGA(I)=OMGA(I)
    end do
  else
    do I=1, NCONS
      T%OMGA(i)=OMGA0
    end do
  end if

  if (T% LOPT_OMGB) then
    do I=1, NCONS
      T%OMGB(I)=OMGB(I)
    end do
  else
    do I=1, NCONS
      T%OMGB(i)=OMGB0
    end do
  end if

end subroutine EosIfp_Param_CreateFromFile

end module M_T_EosIfp_Param
