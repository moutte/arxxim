module M_Solmodel_Vars
  
  use M_Kinds
  use M_T_DtbLogKTbl,only: DimLogK_Max
  
  implicit none
  
  private
  
  !----------------------------------------------------- Solmodel_Vars --
  logical,public:: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot

  type,public:: T_Spline
    
    !----------- tabulated values (input)
    integer :: Dimm
    real(dp):: vX(DimLogK_Max)
    real(dp):: vY(DimLogK_Max)

    !----------- Spline Coeffs 
    real(dp):: vSplineB(DimLogK_Max)
    real(dp):: vSplineC(DimLogK_Max)
    real(dp):: vsplineD(DimLogK_Max)
    
  end type T_Spline

  type(T_Spline),public:: &
  & Rho_Spl, &
  & Eps_Spl, &
  & dhA_Spl, &
  & dhB_Spl, &
  & bDot_Spl
  !----------------------------------------------------/ Solmodel_Vars --
  
end module M_Solmodel_Vars

