MODULE M_Dtb_Vars
!--
!-- module of thermodynamic data on min./gas species
!--
!-- TODO:
!-- all databases should include data on validity range,
!-- i.e. Tmin,Tmax + Pmin,Pmax
!-- temperature range should also be included for each T_SpeciesDtb
!--
  USE M_Kinds
  USE M_T_DtbMinHkf, ONLY: T_DtbMinHkf
  USE M_T_DtbMinThr, ONLY: T_DtbMinThr
  USE M_T_DtbAquHkf, ONLY: T_DtbAquHkf
  USE M_T_DtbH2OHkf, ONLY: T_DtbH2OHkf
  USE M_T_DtbLogKTbl,ONLY: T_DtbLogKTbl,DimLogK_Max
  USE M_T_DtbLogKAnl,ONLY: T_DtbLogKAnl
  !
  USE M_T_Tpcond,   ONLY: T_TPCond
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  !~ TYPE T_SpeciesDtb
    !~ CHARACTER(LEN=7):: EoS !"AQU_HKF","MIN_HKF","LOGKTBL","LOGKANL",etc
    !~ INTEGER:: Indx !index of pure species in the corresponding vDtbEoS
  !~ ENDTYPE
  !~ TYPE(T_SpeciesDtb),ALLOCATABLE:: vSpeciesDtb(:)
  !
  TYPE(T_DtbMinHkf),  ALLOCATABLE:: vDtbMinHkf(:)
  TYPE(T_DtbMinThr),  ALLOCATABLE:: vDtbMinThr(:)
  TYPE(T_DtbAquHkf),  ALLOCATABLE:: vDtbAquHkf(:)
  TYPE(T_DtbLogKTbl), ALLOCATABLE:: vDtbLogkTbl(:)
  TYPE(T_DtbLogKAnl), ALLOCATABLE:: vDtbLogkAnl(:)
  !
  TYPE(T_TPCond), ALLOCATABLE:: DtbLogK_vTPCond(:)
  !
  LOGICAL,PARAMETER:: Psat_Auto = .TRUE.
  !
  CHARACTER(LEN=7),PUBLIC:: DtbFormat 
  != "LOGK","LOGKTBL","LOGKANL","HSV.THR",...
  INTEGER:: DtbLogK_Dim
  
  !! LOGICAL,PUBLIC:: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  !! 
  !! TYPE:: T_Spline
  !!   
  !!   !----------- tabulated values (input)
  !!   INTEGER :: Dimm
  !!   REAL(dp):: vX(DimLogK_Max)
  !!   REAL(dp):: vY(DimLogK_Max)
  !! 
  !!   !----------- Spline Coeffs 
  !!   REAL(dp):: vSplineB(DimLogK_Max)
  !!   REAL(dp):: vSplineC(DimLogK_Max)
  !!   REAL(dp):: vsplineD(DimLogK_Max)
  !!   
  !! END TYPE T_Spline
  !! 
  !! TYPE(T_Spline):: &
  !! & Rho_Spl, &
  !! & Eps_Spl, &
  !! & dhA_Spl, &
  !! & dhB_Spl, &
  !! & bDot_Spl
  
CONTAINS

SUBROUTINE Dtb_Vars_Zero
  !
  CALL Dtb_Vars_Clean
  !
  ALLOCATE(vDtbMinHkf(0))
  ALLOCATE(vDtbMinThr(0))
  ALLOCATE(vDtbAquHkf(0))
  ALLOCATE(vDtbLogkTbl(0))
  ALLOCATE(vDtbLogkAnl(0))
  !
  ALLOCATE(DtbLogK_vTPCond(0))
  !
  ! ALLOCATE(vDtb_MeltGhio(0))
  ! ALLOCATE(vDtbSpLogK(0))
  !
ENDSUBROUTINE Dtb_Vars_Zero

SUBROUTINE Dtb_Vars_Clean
  !
  IF(ALLOCATED(vDtbMinHkf))  DEALLOCATE(vDtbMinHkf)
  IF(ALLOCATED(vDtbMinThr))  DEALLOCATE(vDtbMinThr)
  IF(ALLOCATED(vDtbAquHkf))  DEALLOCATE(vDtbAquHkf)
  IF(ALLOCATED(vDtbLogkTbl)) DEALLOCATE(vDtbLogkTbl)
  IF(ALLOCATED(vDtbLogkAnl)) DEALLOCATE(vDtbLogkAnl)
  !
  IF(ALLOCATED(DtbLogK_vTPCond)) DEALLOCATE(DtbLogK_vTPCond)
  !
  !!IF(ALLOCATED(vDtb_MeltGhio)) DEALLOCATE(vDtb_MeltGhio)
  !!IF(ALLOCATED(vDtbSpLogK)) DEALLOCATE(vDtbSpLogK)
  !
ENDSUBROUTINE Dtb_Vars_Clean

ENDMODULE M_Dtb_Vars

