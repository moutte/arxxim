MODULE M_T_DtbLogKTbl
  !--------------------------------------------------------------
  ! Purpose : DATA TYPE T_DtbLogKTbl
  ! thermodynamic data given as logK's along a TP.TABLE
  !--------------------------------------------------------------
  USE M_Kinds
  USE M_Trace
  USE M_T_Tpcond,ONLY: T_TPCond
  !
  IMPLICIT NONE
  PRIVATE
  !
  !//Public derived type
  PUBLIC:: T_DtbLogKTbl
  !
  !//Public variables
  PUBLIC:: DimLogK_Max
  INTEGER,PARAMETER:: DimLogK_Max= 16
  !
  !//Public functions
  PUBLIC  :: DtbLogKTbl_Calc_Init
  PUBLIC  :: DtbLogKTbl_Calc
  !
  TYPE:: T_DtbLogKTbl
    !
    CHARACTER(LEN=15):: Num
    CHARACTER(LEN=23):: Name
    CHARACTER(LEN=71):: Formula
    CHARACTER(LEN=3) :: Typ     ! AQU:MIN:GAS
    INTEGER :: Div=1
    INTEGER :: Chg
    REAL(dp):: WeitKg
    REAL(dp):: V0R ! Volume IF MIN
    REAL(dp):: AquSize

    !----------- LogK Table Mod
    INTEGER :: DimLogK
    REAL(dp):: vTdgK(DimLogK_Max)
    REAL(dp):: vLogK(DimLogK_Max)

    !----------- Spline Coeffs
    REAL(dp):: vSplineB(DimLogK_Max)
    REAL(dp):: vSplineC(DimLogK_Max)
    REAL(dp):: vsplineD(DimLogK_Max)

    !---------- Memory of current state
    !! REAL(dp):: G0rt
    !
  END TYPE T_DtbLogKTbl

  !// Private Functions
  PRIVATE :: DtbLogKTbl_LogKSpline_Compute
  PRIVATE :: DtbLogKTbl_LogKSpline_Eval

CONTAINS

  SUBROUTINE DtbLogKTbl_Calc_Init(M)
    !---------------------------------------------------
    ! Purpose : Prepare Spline Interpolator
    !---------------------------------------------------
    TYPE(T_DtbLogKTbl),INTENT(INOUT) :: M

    IF(M%DimLogK>1) CALL DtbLogKTbl_LogKSpline_Compute(M)

  END SUBROUTINE DtbLogKTbl_Calc_Init

  !---

  SUBROUTINE DtbLogKTbl_Calc(M,TdgK,Pbar,S)
    !---------------------------------------------------
    ! Purpose : update values of S%G0RT,S%V0
    !---------------------------------------------------
    USE M_T_Species,    ONLY: T_Species
    USE M_Dtb_Const,    ONLY: R_jK
    USE M_Numeric_Const,ONLY: Ln10
    !! USE M_T_Tpcond,     ONLY: TPcond_IndexTP

    TYPE(T_DtbLogKTbl),INTENT(IN):: M
    REAL(dp),          INTENT(IN):: TdgK,Pbar
    TYPE(T_Species),   INTENT(INOUT):: S

    REAL(dp):: LogK

    S%WeitKg= M%WeitKg

    !// ref state potential
    IF(M%DimLogK>1) THEN
      CALL DtbLogKTbl_LogKSpline_Eval(M,TdgK,LogK)
    ELSE
      LogK= M%vLogK(1)
    ENDIF
    S%G0rt= -LogK *Ln10

    !// molar volume
    SELECT CASE(TRIM(M%Typ))
    CASE("MIN")
      S%V0= M%V0R
    CASE("GAS")
      S%V0= R_jK *TdgK / Pbar/ 1.D5 ! -> perfect gas molar volume
    CASE("AQU")
      S%AquSize= M%AquSize
    END SELECT

  END SUBROUTINE DtbLogKTbl_Calc

  !---

  SUBROUTINE DtbLogKTbl_LogKSpline_Compute(M)
    !---------------------------------------------------
    ! Purpose : compute Spline Coeffs
    !---------------------------------------------------
    USE M_CMM_Spline
    USE M_IoTools,ONLY: GetUnit

    TYPE(T_DtbLogKTbl), INTENT(INOUT) :: M

    INTEGER :: I, n
    REAL(dp),ALLOCATABLE :: B(:), C(:), D(:)
    INTEGER :: ff
    !---
    n= M%DimLogK


    IF(N==1) THEN
      M% vSplineB(1:N)= zero
      M% vSplineC(1:N)= zero
      M% vSplineD(1:N)= zero
      RETURN
    ENDIF

    FF= 0
    IF(iDebug==4) THEN
      CALL GetUnit(FF)
      OPEN(FF,file="debug_spline.log")
    ENDIF

    !// ALLOCATE temporary vectors for output
    ALLOCATE(B(n))
    ALLOCATE(C(n))
    ALLOCATE(D(n))

    !// compute spline coefficients
    CALL CMM_Spline_Compute (N,M%vTdgK(1:N),M%vLogK(1:N),B,C,D)
    !! CALL Spline_Compute_NR(M%vTdgK(1:n),M%vLogK(1:n),B)

    M% vSplineB(1:n)= B(1:n)
    M% vSplineC(1:n)= C(1:n)
    M% vSplineD(1:n)= D(1:n)

    !// Debug
    IF(ff>0) THEN
      DO I=1, n
        WRITE(FF,'(F8.3,4G15.6)') M%vTdgK(I),M%vLogK(I),B(I),C(I),D(I)
      END DO
    ENDIF

    !// clean temporary vectors
    DEALLOCATE(B,C,D)
    !
    IF(FF>0) CLOSE(FF)
  END SUBROUTINE DtbLogKTbl_LogKSpline_Compute

  !---

  SUBROUTINE DtbLogKTbl_LogKSpline_Eval(M,TdgK,LogK)
    !---------------------------------------------------
    ! Purpose : eval LogK with Spline Function
    !---------------------------------------------------
    USE M_CMM_Spline
    !! USE M_Numeric_Interpol
    !
    TYPE(T_DtbLogKTbl), INTENT(IN) :: M
    REAL(dp), INTENT(IN)  :: TdgK
    REAL(dp), INTENT(OUT) :: LogK

    INTEGER :: n
    !---
    n= M%DimLogK

    LogK= Zero

    CALL CMM_Spline_Eval(n, &
    & M%vTdgK(1:n) , M%vLogK(1:n),&
    & M%vSplineB(1:n) , M%vSplineC(1:n), M%vSplineD(1:n),&
    & TdgK, LogK )

    !~ LogK= Spline_Eval_NR( &
    !~ & M%vTdgK(1:n) , M%vLogK(1:n), M%vSplineB(1:n), TdgK)

  END SUBROUTINE DtbLogKTbl_LogKSpline_Eval

ENDMODULE M_T_DtbLogKTbl

