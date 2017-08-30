module M_Dtb_Vars
!--
!-- module of thermodynamic data on min./gas species
!--
!-- TOdo:
!-- all databases should include data on validity range,
!-- i.e. Tmin,Tmax + Pmin,Pmax
!-- temperature range should also be included for each T_SpeciesDtb
!--
  use M_Kinds
  use M_T_DtbMinHkf, only: T_DtbMinHkf
  use M_T_DtbMinThr, only: T_DtbMinThr
  use M_T_DtbAquHkf, only: T_DtbAquHkf
  use M_T_DtbH2OHkf, only: T_DtbH2OHkf
  !
  use M_T_DtbLogKTbl,only: T_DtbLogKTbl,DimLogK_Max
  use M_T_DtbLogKAnl,only: T_DtbLogKAnl
  !
  use M_T_Tpcond,   only: T_TPCond
  !
  implicit none
  !
  public
  !
  !! type T_SpeciesDtb
  !!   character(len=7):: EoS !"AQU_HKF","MIN_HKF","LOGKTBL","LOGKANL",etc
  !!   integer:: Indx !index of pure species in the corresponding vDtbEoS
  !! end type
  !! type(T_SpeciesDtb),allocatable:: vSpeciesDtb(:)
  !
  type(T_DtbMinHkf),  allocatable:: vDtbMinHkf(:)
  type(T_DtbMinThr),  allocatable:: vDtbMinThr(:)
  type(T_DtbAquHkf),  allocatable:: vDtbAquHkf(:)
  type(T_DtbLogKTbl), allocatable:: vDtbLogkTbl(:)
  type(T_DtbLogKAnl), allocatable:: vDtbLogkAnl(:)
  !
  type(T_TPCond), allocatable:: DtbLogK_vTPCond(:)
  !
  logical,parameter:: Psat_Auto = .true.
  !
  character(len=7),public:: DtbFormat 
  != "LOGK","LOGKTBL","LOGKANL","HSV.THR",...
  integer:: DtbLogK_Dim
  
  !! logical,public:: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  !! 
  !! type:: T_Spline
  !!   
  !!   !----------- tabulated values (input)
  !!   integer :: Dimm
  !!   real(dp):: vX(DimLogK_Max)
  !!   real(dp):: vY(DimLogK_Max)
  !! 
  !!   !----------- Spline Coeffs 
  !!   real(dp):: vSplineB(DimLogK_Max)
  !!   real(dp):: vSplineC(DimLogK_Max)
  !!   real(dp):: vsplineD(DimLogK_Max)
  !!   
  !! end type T_Spline
  !! 
  !! type(T_Spline):: &
  !! & Rho_Spl, &
  !! & Eps_Spl, &
  !! & dhA_Spl, &
  !! & dhB_Spl, &
  !! & bDot_Spl
  
contains

subroutine Dtb_Vars_Zero
  !
  call Dtb_Vars_Clean
  !
  allocate(vDtbMinHkf(0))
  allocate(vDtbMinThr(0))
  allocate(vDtbAquHkf(0))
  allocate(vDtbLogkTbl(0))
  allocate(vDtbLogkAnl(0))
  !
  allocate(DtbLogK_vTPCond(0))

  ! allocate(vDtb_MeltGhio(0))
  ! allocate(vDtbSpLogK(0))
  !
end subroutine Dtb_Vars_Zero

subroutine Dtb_Vars_Clean
  !
  if(allocated(vDtbMinHkf))  deallocate(vDtbMinHkf)
  if(allocated(vDtbMinThr))  deallocate(vDtbMinThr)
  if(allocated(vDtbAquHkf))  deallocate(vDtbAquHkf)
  if(allocated(vDtbLogkTbl)) deallocate(vDtbLogkTbl)
  if(allocated(vDtbLogkAnl)) deallocate(vDtbLogkAnl)
  !
  if(allocated(DtbLogK_vTPCond)) deallocate(DtbLogK_vTPCond)
  !
  !!if(allocated(vDtb_MeltGhio)) deallocate(vDtb_MeltGhio)
  !!if(allocated(vDtbSpLogK)) deallocate(vDtbSpLogK)
  !
end subroutine Dtb_Vars_Clean

end module M_Dtb_Vars

