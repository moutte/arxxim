MODULE M_T_DtbLogKAnl
  !--------------------------------------------------------------
  ! Purpose : DATA TYPE T_DtbLogKAnl
  ! thermodynamic data given as logK's along a TP.TABLE 
  !--------------------------------------------------------------
  USE M_Kinds
  USE M_Trace
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC:: T_DtbLogKAnl
  
  TYPE:: T_DtbLogKAnl
    
    CHARACTER(LEN=15):: Num !
    CHARACTER(LEN=23):: Name
    CHARACTER(LEN=71):: Formula
    !CHARACTER(LEN=6) :: Fitting
    CHARACTER(LEN=3) :: Typ       ! AQU:MIN:GAS
    INTEGER  :: iFitting
    INTEGER  :: Div
    REAL(dp) :: WeitKg
    INTEGER  :: Chg
    REAL(dp) :: V0R               ! Volume IF MIN
    REAL(dp) :: AquSize           ! Radius IF AQU
    
    !----------- Fitting Coeffs 
    REAL(dp):: vX(8) !!A, B, C, D, E, F
    !! REAL(dp),ALLOCATABLE:: vX(:) -> the future ...
    
  END TYPE T_DtbLogKAnl

  !// Public Functions
  PUBLIC  :: DtbLogKAnl_Calc
  PUBLIC  :: DtbLogKAnl_New
  
CONTAINS

  SUBROUTINE DtbLogKAnl_New(M)
    TYPE(T_DtbLogKAnl),INTENT(OUT) :: M
    !
    M%Num="NONE"
    M%Name="NONE"
    M%Formula="NONE"
    M%Typ="NON"
    !
    M%iFitting= 0
    M%Div= 1
    M%WeitKg= 1.0D0
    M%Chg= 0
    M%V0R= 1.0D0
    M%AquSize= Zero
    M%vX(1:8)= Zero
    !
  END SUBROUTINE DtbLogKAnl_New
  
  SUBROUTINE DtbLogKAnl_Calc(M,TdgK,Pbar,S)
    !---------------------------------------------------
    ! Purpose : update values of G0RT(S)
    !---------------------------------------------------
    USE M_T_Species
    USE M_Numeric_Const,ONLY: Ln10
    USE M_Dtb_Const,ONLY: R_jK, Pref, T_CK
    !---
    TYPE(T_DtbLogKAnl),INTENT(IN) :: M
    REAL(dp),          INTENT(IN) :: TdgK,Pbar
    TYPE(T_Species),   INTENT(OUT):: S 
    !---
    REAL(dp):: LogK, T
    !---
    T = TdgK
    !
    !// ref state potential
    SELECT CASE(M%iFitting)
    
    CASE(0) ! fixed parameter
      LogK = M%vX(1)
    
    CASE(1) !"PHREEQC"
    !~ IF(TRIM(M%Fitting)=="PHREEQ") THEN
      !~ T = TdgK
      !~ LogK = M%A               &
      !~ &    + M%B * T           &
      !~ &    + M%C * 1.D0/T      &
      !~ &    + M%D * Log(T)/LN10 &
      !~ &    + M%E * 1.D0/T/T
    !~ ENDIF
      LogK = M%vX(1)               &
      &    + M%vX(2) * T           &
      &    + M%vX(3) * 1.D0/T      &
      &    + M%vX(4) * Log(T)/LN10 &
      &    + M%vX(5) * 1.D0/T/T
    
    CASE(2) !"CHRISTOV"
      LogK = M%vX(1)               &
      &    + M%vX(2) * T           &
      &    + M%vX(3) * 1.D0/T      &
      &    + M%vX(4) * Log(T)/LN10 &
      &    + M%vX(5) * 1.D0/T/T
    
    END SELECT
    !
    S%G0rt = -LogK *Ln10
    !
    !// molar volume
    SELECT CASE(TRIM(M%Typ))
    CASE("MIN")
      S%V0 = M%V0R                    ! molar vol. in m^3       
    CASE("GAS")       
      S%V0 =  R_jK *TdgK / Pbar/ 1.D5 ! -> perfect gas molar volume
    !CASE("AQU")
    !  S%AquSize= M%AquSize
    END SELECT
    
  END SUBROUTINE DtbLogKAnl_Calc
  
ENDMODULE M_T_DtbLogKAnl

