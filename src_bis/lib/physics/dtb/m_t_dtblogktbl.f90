module M_T_DtbLogKTbl
  !--------------------------------------------------------------
  ! Purpose : data type T_DtbLogKTbl
  ! thermodynamic data given as logK's along a TP.TABLE
  !--------------------------------------------------------------
  use M_Kinds
  use M_Trace
  use M_T_Tpcond,only: T_TPCond
  !
  implicit none
  private
  !
  !//Public derived type
  public:: T_DtbLogKTbl
  !
  !//Public variables
  public:: DimLogK_Max
  integer,parameter:: DimLogK_Max= 16
  !
  !//Public functions
  public  :: DtbLogKTbl_Calc_Init
  public  :: DtbLogKTbl_Calc
  !
  type:: T_DtbLogKTbl
    !
    character(len=15):: Num
    character(len=23):: Name
    character(len=71):: Formula
    character(len=3) :: Typ     ! AQU:MIN:GAS:XCH
    integer :: Div=1
    integer :: Chg
    real(dp):: WeitKg
    real(dp):: V0R ! Volume if MIN
    real(dp):: AquSize

    !----------- LogK Table Mod
    integer :: DimLogK
    real(dp):: vTdgK(DimLogK_Max)
    real(dp):: vLogK(DimLogK_Max)

    !----------- Spline Coeffs
    real(dp):: vSplineB(DimLogK_Max)
    real(dp):: vSplineC(DimLogK_Max)
    real(dp):: vsplineD(DimLogK_Max)

    !---------- Memory of current state
    !! real(dp):: G0rt
    !
  end type T_DtbLogKTbl

  !// Private Functions
  private :: DtbLogKTbl_LogKSpline_Compute
  private :: DtbLogKTbl_LogKSpline_Eval

contains

  subroutine DtbLogKTbl_Calc_Init(M)
    !---------------------------------------------------
    ! Purpose : Prepare Spline Interpolator
    !---------------------------------------------------
    type(T_DtbLogKTbl),intent(inout) :: M

    if(M%DimLogK>1) call DtbLogKTbl_LogKSpline_Compute(M)

  end subroutine DtbLogKTbl_Calc_Init

  !---

  subroutine DtbLogKTbl_Calc(M,TdgK,Pbar,S)
    !---------------------------------------------------
    ! Purpose : update values of S%G0RT,S%V0
    !---------------------------------------------------
    use M_T_Species,    only: T_Species
    use M_Dtb_Const,    only: R_jK
    use M_Numeric_Const,only: Ln10
    !! use M_T_Tpcond,     only: TPcond_IndexTP

    type(T_DtbLogKTbl),intent(in):: M
    real(dp),          intent(in):: TdgK,Pbar
    type(T_Species),   intent(inout):: S

    real(dp):: LogK

    S%WeitKg= M%WeitKg

    !// ref state potential
    if(M%DimLogK>1) then
      call DtbLogKTbl_LogKSpline_Eval(M,TdgK,LogK)
    else
      LogK= M%vLogK(1)
    end if
    S%G0rt= -LogK *Ln10

    !// molar volume
    select case(trim(M%Typ))
    case("MIN")
      S%V0= M%V0R
    case("GAS")
      S%V0= R_jK *TdgK / Pbar/ 1.D5 ! -> perfect gas molar volume
    case("AQU")
      S%AquSize= M%AquSize
    end select

  end subroutine DtbLogKTbl_Calc

  !---

  subroutine DtbLogKTbl_LogKSpline_Compute(M)
    !---------------------------------------------------
    ! Purpose : compute Spline Coeffs
    !---------------------------------------------------
    use M_CMM_Spline
    use M_IoTools,only: GetUnit

    type(T_DtbLogKTbl), intent(inout) :: M

    integer :: I, n
    real(dp),allocatable :: B(:), C(:), D(:)
    integer :: ff
    !---
    n= M%DimLogK


    if(N==1) then
      M% vSplineB(1:N)= zero
      M% vSplineC(1:N)= zero
      M% vSplineD(1:N)= zero
      return
    end if

    FF= 0
    if(iDebug==4) then
      call GetUnit(FF)
      open(FF,file="debug_spline.log")
    end if

    !// allocate temporary vectors for output
    allocate(B(n))
    allocate(C(n))
    allocate(D(n))

    !// compute spline coefficients
    call CMM_Spline_Compute (N,M%vTdgK(1:N),M%vLogK(1:N),B,C,D)
    !! call Spline_Compute_NR(M%vTdgK(1:n),M%vLogK(1:n),B)

    M% vSplineB(1:n)= B(1:n)
    M% vSplineC(1:n)= C(1:n)
    M% vSplineD(1:n)= D(1:n)

    !// Debug
    if(ff>0) then
      do I=1, n
        write(FF,'(F8.3,4G15.6)') M%vTdgK(I),M%vLogK(I),B(I),C(I),D(I)
      end do
    end if

    !// clean temporary vectors
    deallocate(B,C,D)
    !
    if(FF>0) close(FF)
  end subroutine DtbLogKTbl_LogKSpline_Compute

  !---

  subroutine DtbLogKTbl_LogKSpline_Eval(M,TdgK,LogK)
    !---------------------------------------------------
    ! Purpose : eval LogK with Spline Function
    !---------------------------------------------------
    use M_CMM_Spline
    !! use M_Numeric_Interpol
    !
    type(T_DtbLogKTbl), intent(in) :: M
    real(dp), intent(in)  :: TdgK
    real(dp), intent(out) :: LogK

    integer :: n
    !---
    n= M%DimLogK

    LogK= Zero

    call CMM_Spline_Eval(n, &
    & M%vTdgK(1:n) , M%vLogK(1:n),&
    & M%vSplineB(1:n) , M%vSplineC(1:n), M%vSplineD(1:n),&
    & TdgK, LogK )

    !! LogK= Spline_Eval_NR( &
    !! & M%vTdgK(1:n) , M%vLogK(1:n), M%vSplineB(1:n), TdgK)

  end subroutine DtbLogKTbl_LogKSpline_Eval

end module M_T_DtbLogKTbl

