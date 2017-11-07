module M_T_DtbLogKAnl
  !--------------------------------------------------------------
  ! Purpose : data type T_DtbLogKAnl
  ! thermodynamic data given as logK's along a TP.TABLE 
  !--------------------------------------------------------------
  use M_Kinds
  use M_Trace
  
  implicit none
  private
  
  public:: T_DtbLogKAnl
  
  type:: T_DtbLogKAnl
    
    character(len=15):: Num !
    character(len=23):: Name
    character(len=15):: Source
    character(len=5) :: Abbr
    character(len=71):: Formula
    !character(len=6) :: Fitting
    character(len=3) :: Typ       ! AQU:MIN:GAS
    integer  :: iFitting
    integer  :: Div
    real(dp) :: WeitKg
    integer  :: Chg
    real(dp) :: V0R               ! Volume if MIN
    real(dp) :: AquSize           ! Radius if AQU
    
    !----------- Fitting Coeffs 
    real(dp):: vX(8) !!A, B, C, D, E, F
    !! real(dp),allocatable:: vX(:) -> the future ...
    
  end type T_DtbLogKAnl

  !// Public Functions
  public  :: DtbLogKAnl_Calc
  public  :: DtbLogKAnl_New
  
contains

  subroutine DtbLogKAnl_New(M)
    type(T_DtbLogKAnl),intent(out) :: M
    !
    M%Num=     "NONE"
    M%Name=    "NONE"
    M%Source=  "NONE"
    M%Formula= "NONE"
    M%Typ=     "NON"
    !
    M%iFitting= 0
    M%Div= 1
    M%WeitKg= 1.0D0
    M%Chg= 0
    M%V0R= 1.0D0
    M%AquSize= Zero
    M%vX(1:8)= Zero
    !
  end subroutine DtbLogKAnl_New
  
  subroutine DtbLogKAnl_Calc(M,TdgK,Pbar,S)
    !---------------------------------------------------
    ! Purpose : update values of G0RT(S)
    !---------------------------------------------------
    use M_T_Species
    use M_Numeric_Const,only: Ln10
    use M_Dtb_Const,only: R_jK, Pref, T_CK
    !---
    type(T_DtbLogKAnl),intent(in) :: M
    real(dp),          intent(in) :: TdgK,Pbar
    type(T_Species),   intent(out):: S 
    !---
    real(dp):: LogK, T
    !---
    T = TdgK
    !
    !// ref state potential
    select case(M%iFitting)
    
    case(0) ! fixed parameter
      LogK = M%vX(1)
    
    case(1) !"PHREEQC"
      LogK = M%vX(1)               &
      &    + M%vX(2) * T           &
      &    + M%vX(3) / T           &
      &    + M%vX(4) * Log(T)/LN10 &
      &    + M%vX(5) / T**2
    
    case(2) !"CHRISTOV"
      LogK = M%vX(1)               &
      &    + M%vX(2) * T           &
      &    + M%vX(3) / T           &
      &    + M%vX(4) * Log(T)/LN10 &
      &    + M%vX(5) / T**2
    
    end select
    !
    S%G0rt = -LogK *Ln10
    !
    !// molar volume
    select case(trim(M%Typ))
    case("MIN")
      S%V0 = M%V0R                    ! molar vol. in m^3       
    case("GAS")       
      S%V0 =  R_jK *TdgK / Pbar/ 1.D5 ! -> perfect gas molar volume
    !case("AQU")
    !  S%AquSize= M%AquSize
    end select
    
  end subroutine DtbLogKAnl_Calc
  
end module M_T_DtbLogKAnl

