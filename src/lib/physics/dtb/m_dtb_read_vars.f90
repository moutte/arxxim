MODULE M_Dtb_Read_Vars
!--
!-- variables for routines reading databases
!--
  USE M_T_DtbAquHkf, ONLY: T_DtbAquHkf
  USE M_T_DtbMinHkf, ONLY: T_DtbMinHkf
  USE M_T_DtbMinThr, ONLY: T_DtbMinThr
  USE M_T_DtbLogKTbl,ONLY: T_DtbLogKTbl
  USE M_T_DtbLogKAnl,ONLY: T_DtbLogKAnl
  !! USE M_T_DtbHSV,    ONLY: T_DtbHSV
  !
  IMPLICIT NONE
  !
  TYPE:: T_LisAquHkf !linked list for "dynamic reading"
    TYPE(T_DtbAquHkf)::Value
    TYPE(T_LisAquHkf),POINTER::Next
  ENDTYPE T_LisAquHkf
  TYPE(T_LisAquHkf),POINTER:: LisAquHkf
  INTEGER:: nAquHkf
  !
  TYPE:: T_LisMinHkf !linked list for "dynamic reading"
    TYPE(T_DtbMinHkf)::Value
    TYPE(T_LisMinHkf),POINTER::Next
  ENDTYPE T_LisMinHkf
  TYPE(T_LisMinHkf),POINTER:: LisMinHkf
  INTEGER:: nMinHkf
  !
  TYPE:: T_LisMinThr
    TYPE(T_DtbMinThr)::Value
    TYPE(T_LisMinThr),POINTER::Next
  ENDTYPE T_LisMinThr
  TYPE(T_LisMinThr),POINTER:: LisMinThr
  INTEGER:: nMinThr
  !
  !~ TYPE:: T_LisHSV
    !~ TYPE(T_DtbHSV)::Value
    !~ TYPE(T_LisHSV),POINTER::Next
  !~ ENDTYPE T_LisHSV
  !~ TYPE(T_LisHSV),POINTER:: LisHSV
  !~ INTEGER:: nHSV
  !
  TYPE:: T_LisLogKTbl
    TYPE(T_DtbLogKTbl)::Value
    TYPE(T_LisLogKTbl),POINTER::Next
  ENDTYPE T_LisLogKTbl
  TYPE(T_LisLogKTbl),POINTER:: LisLogKTbl
  INTEGER:: nLogKTbl
  !
  TYPE:: T_LisLogKAnl
    TYPE(T_DtbLogKAnl)::Value
    TYPE(T_LisLogKAnl),POINTER::Next
  ENDTYPE T_LisLogKAnl
  TYPE(T_LisLogKAnl),POINTER:: LisLogKAnl
  INTEGER:: nLogKAnl
  !
  CHARACTER(LEN=10):: FilCode
  !
  INTEGER:: fLin=0
  !
ENDMODULE M_Dtb_Read_Vars
