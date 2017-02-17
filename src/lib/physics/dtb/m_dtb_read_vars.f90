module M_Dtb_Read_Vars
!--
!-- variables for routines reading databases
!--
  use M_T_DtbAquHkf, only: T_DtbAquHkf
  use M_T_DtbMinHkf, only: T_DtbMinHkf
  use M_T_DtbMinThr, only: T_DtbMinThr
  use M_T_DtbLogKTbl,only: T_DtbLogKTbl
  use M_T_DtbLogKAnl,only: T_DtbLogKAnl
  !! use M_T_DtbHSV,    only: T_DtbHSV
  !
  implicit none
  !
  type:: T_LisAquHkf !linked list for "dynamic reading"
    type(T_DtbAquHkf)::Value
    type(T_LisAquHkf),pointer::Next
  end type T_LisAquHkf
  type(T_LisAquHkf),pointer:: LisAquHkf
  integer:: nAquHkf
  !
  type:: T_LisMinHkf !linked list for "dynamic reading"
    type(T_DtbMinHkf)::Value
    type(T_LisMinHkf),pointer::Next
  end type T_LisMinHkf
  type(T_LisMinHkf),pointer:: LisMinHkf
  integer:: nMinHkf
  !
  type:: T_LisMinThr
    type(T_DtbMinThr)::Value
    type(T_LisMinThr),pointer::Next
  end type T_LisMinThr
  type(T_LisMinThr),pointer:: LisMinThr
  integer:: nMinThr
  !
  !! type:: T_LisHSV
  !!   type(T_DtbHSV)::Value
  !!   type(T_LisHSV),pointer::Next
  !! end type T_LisHSV
  !! type(T_LisHSV),pointer:: LisHSV
  !!   integer:: nHSV
  !
  type:: T_LisLogKTbl
    type(T_DtbLogKTbl)::Value
    type(T_LisLogKTbl),pointer::Next
  end type T_LisLogKTbl
  type(T_LisLogKTbl),pointer:: LisLogKTbl
  integer:: nLogKTbl
  !
  type:: T_LisLogKAnl
    type(T_DtbLogKAnl)::Value
    type(T_LisLogKAnl),pointer::Next
  end type T_LisLogKAnl
  type(T_LisLogKAnl),pointer:: LisLogKAnl
  integer:: nLogKAnl
  !
  character(len=10):: FilCode
  !
  integer:: fLin=0
  !
end module M_Dtb_Read_Vars
