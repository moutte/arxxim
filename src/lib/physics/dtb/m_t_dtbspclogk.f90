MODULE M_T_DtbSpcLogK
!.DATA TYPE T_DtbSpcLogK for thermodynamic DATA as logK's along a TP.TABLE
  USE M_Kinds
  USE M_Trace
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_DtbSpcLogK
  PUBLIC:: DimLogK_Max
  !
  INTEGER,PARAMETER:: DimLogK_Max = 16
  !
  TYPE:: T_DtbSpcLogK 
  ! DATA structure for species following a logK format
    CHARACTER(LEN=15):: Num
    CHARACTER(LEN=23):: Name
    CHARACTER(LEN=71):: Formula
    CHARACTER(LEN=3) :: Typ !MIN:GAS
    !
    INTEGER :: iMode
    INTEGER :: DimLogK
    ! iMode=1 -> table of discrete values = tLogK
    ! iMode=2 -> coeffs of FUNCTION = A,B,C,D
    REAL(dp):: tLogK(1:DimLogK_Max)
    REAL(dp):: A1,A2,A3,A4,A5
  ENDTYPE T_DtbSpcLogK

ENDMODULE M_T_DtbSpcLogK

