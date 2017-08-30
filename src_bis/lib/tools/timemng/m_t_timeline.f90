module M_T_TimeLine
  
  !=========================================================
  ! T_TimeLine, Time Domain of Interest =  [TInit, TFinal]
  !=========================================================

  use M_Kinds
  implicit none
  private

  public ::  T_TimeLine
  type T_TimeLine
     real(dp) :: TInit
     real(dp) :: TFinal
  end type T_TimeLine

  !// public Functions
  public :: TimeLine_Zero
  public :: TimeLine_Init
  public :: TimeLine_Check
  public :: TimeLine_Info
  
contains
 
!---
  
  subroutine TimeLine_Zero( T)
    !===================================================
    ! Purpose : Init TimeLine by using Zero values
    !===================================================
    implicit none
    !
    type(T_TimeLine) :: T
    !-----
    T% TInit  = Zero
    T% TFinal = Zero

  end subroutine TimeLine_Zero

!---

  subroutine TimeLine_Init( T, TInit, TFinal)
    !===================================================
    ! Purpose : Init TimeLine by using given parameters
    !===================================================
    implicit none
    !
    real(dp), intent(in) :: TInit, TFinal
    type(T_TimeLine) :: T
    !-----
    T% TInit  = Tinit
    T% TFinal = TFinal

  end subroutine TimeLine_Init

   !---
  
  subroutine TimeLine_Check( T )
    !===============================================
    ! Check TimeMng Minimal Monotonicity Assumptions
    !----
    ! TInit < TFinal
    !===============================================
    use M_Trace
    implicit none
    type(T_TimeLine) :: T
    !---
    if ( T%TFinal < T% TInit ) then
       call Fatal_("Wrong TimeLine : TFinal < TInit")
    end if

  end subroutine TimeLine_Check

  !---
  
  subroutine TimeLine_Info( T )
    !===============================================
    ! Check TimeMng Minimal Monotonicity Assumptions
    !----
    ! TInit < TFinal
    !===============================================
    use M_Trace
    implicit none
    type(T_TimeLine) :: T
    !---
    call Info_("----- TimeLine Object")
    Write(*,*) "TInit=", T% TInit
    Write(*,*) "TFinal=", T% TFinal

  end subroutine TimeLine_Info

 
end module M_T_TimeLine
