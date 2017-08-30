module M_T_TimeStepMng

  !======================================================
  ! T_TimeStepMng , Strategy for Time Management
  !======================================================

  use M_Kinds 
  implicit none
  private

  real(DP), parameter :: Default_DTRatio_Increase = 2.*One
  real(DP), parameter :: Default_DTRatio_Reduce   = 0.5*One
  
  public :: T_TimeStepMng
  type T_TimeStepMng 
     real(dp) :: DTInit
     real(dp) :: DTMin
     real(dp) :: DTMax  
     real(dp) :: DTRatio_Increase
     real(dp) :: DTRatio_Reduce     
  end type T_TimeStepMng
  
  !// Public Functions
  public :: TimeStepMng_Zero
  public :: TimeStepMng_Init
  public :: TimeStepMng_Check
  public :: TimeStepMng_Set_DTRatio
  public :: TimeStepMng_Info

contains

!---
  
  subroutine TimeStepMng_Zero(T)
    !===================================================
    ! Purpose : Init TimeStepMng with Zero Values
    !===================================================
    implicit none
    type(T_TimeStepMng) :: T
    !--
    T% DTInit= Zero
    T% DTMax=  Zero
    T% DTMax=  Zero
    
    T% DTRatio_Increase= Zero
    T% DTRatio_Reduce=   Zero
        
  end subroutine TimeStepMng_Zero

  !---

  subroutine TimeStepMng_Init( T, &
       DTInit, DTMin, DTMax, &
       DTRatio_Increase, DTRatio_Reduce )
    !===================================================
    ! Purpose : Init TimeStepMng with User Data
    !===================================================
    implicit none
    type(T_TimeStepMng) :: T
    real(dp), intent(in) :: DTinit
    real(dp), optional, intent(in) :: DTMin, DTMax
    real(dp), optional, intent(in) :: DTRatio_Increase, DTRatio_Reduce    
    !--    
    
    T% DTInit= DTInit  

    T% DTMin= DTInit    
    if (present(DTMin)) T% dTMin = DTMin

    T% DTMax= DTInit
    if (present(DTMax)) T% dTMax = DTMax

    T% DTRATIO_Increase = Default_DTRATIO_Increase
    if (present(DTRATIO_Increase)) T% DTRATIO_Increase = DTRATIO_Increase

    T% DTRATIO_Reduce = Default_DTRATIO_Reduce
    if (present(DTRATIO_Reduce)) T% DTRATIO_Reduce = DTRATIO_Reduce
    
    !// Check TimeMng Validity
    call TimeStepMng_Check(T)

  end subroutine TimeStepMng_Init

  !--

  subroutine TimeStepMng_Set_DTRatio( T, DTRatio_Increase, DTRatio_Reduce ) 
    !===================================================
    ! Purpose : Set TimeStepMng Ratios Only
    !===================================================
    implicit none
    type(T_TimeStepMng) :: T
    real(dp), intent(in) :: DTRatio_Increase, DTRatio_Reduce 
    !-----

    T% DTRatio_Increase= DTRatio_Increase
    T% DTRatio_Reduce=   DTRatio_Reduce
    
    !// Check TimeMng Validity
    call TimeStepMng_Check(T)

  end subroutine TimeStepMng_Set_DTRatio

  !---
  
  subroutine TimeStepMng_Check(T)
    !===============================================
    ! TimeStepMng Minimal Validity Checking
    !----
    ! DTmin <= DTInit <= DTMax
    ! DTRatio_Reduce < 1
    ! DTRatio_Increase > 1
    !===============================================
    use M_Trace
    implicit none
    type(T_TimeStepMng) :: T
    !

    if ( T% DTInit <  T% DTMin ) then
       T% DTMin= T% DTMin
       call Warning_("DTInit < DTmin, TimeStepMng DTMin Corrected")
    end if
    
    if ( T% DTMAx < T% DTInit ) then
       T% DTMax= T% DTInit
       call Warning_("DTMax < DTInit, TimeStepMng DTMax Corrected")
    end if

    if ( T% DTMax < T% DTMin ) then
       call Fatal_("Wrong TimeStepMng, DTMax < DTMin")
    end if

    if ( T% DTRatio_Increase .LE. 1. ) then
       call Fatal_("Wrong TimeStepMng, DTRatio_Increase <= 1")
    end if
    
    if ( T% DTRatio_Reduce .GE. 1. ) then
       call Fatal_("Wrong TimeStepMng, DTRatio_Reduce >= 1")
    end if

  end subroutine TimeStepMng_Check

  !---

  subroutine TimeStepMng_Info(T)
    !===============================================
    ! Information 
    !===============================================
    use M_Trace
    implicit none
    type(T_TimeStepMng) :: T
    !
    call Info_("---- TimeStepMng Object")
    Write(*,*) "DTInit=", T% DTInit
    Write(*,*) "DTMin=",  T% DTMin
    Write(*,*) "DTMax=",  T% DTMAx
    Write(*,*) "DTRatio_Increase=",  T% DTRatio_Increase
    Write(*,*) "DTRatio_Reduce=",  T% DTRatio_Reduce
    
  end subroutine TimeStepMng_Info

end module M_T_TimeStepMng
