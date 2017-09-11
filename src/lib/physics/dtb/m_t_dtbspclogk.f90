module M_T_DtbSpcLogK
!.data type T_DtbSpcLogK for thermodynamic data as logK's along a TP.TABLE
  use M_Kinds
  use M_Trace
  !
  implicit none
  !
  private
  !
  public:: T_DtbSpcLogK
  public:: DimLogK_Max
  !
  integer,parameter:: DimLogK_Max = 16
  !
  type:: T_DtbSpcLogK 
  ! data structure for species following a logK format
    character(len=15):: Num
    character(len=23):: Name
    character(len=5) :: Abbr
    character(len=71):: Formula
    character(len=3) :: Typ !MIN:GAS
    !
    integer :: iMode
    integer :: DimLogK
    ! iMode=1 -> table of discrete values = tLogK
    ! iMode=2 -> coeffs of function = A,B,C,D
    real(dp):: tLogK(1:DimLogK_Max)
    real(dp):: A1,A2,A3,A4,A5
  end type T_DtbSpcLogK

end module M_T_DtbSpcLogK

