MODULE M_Solmodel_Vars
  
  USE M_Kinds
  USE M_T_DtbLogKTbl,ONLY: DimLogK_Max
  
  IMPLICIT NONE
  
  PRIVATE
  
  !----------------------------------------------------- Solmodel_Vars --
  LOGICAL,PUBLIC:: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot

  TYPE,PUBLIC:: T_Spline
    
    !----------- tabulated values (input)
    INTEGER :: Dimm
    REAL(dp):: vX(DimLogK_Max)
    REAL(dp):: vY(DimLogK_Max)

    !----------- Spline Coeffs 
    REAL(dp):: vSplineB(DimLogK_Max)
    REAL(dp):: vSplineC(DimLogK_Max)
    REAL(dp):: vsplineD(DimLogK_Max)
    
  END TYPE T_Spline

  TYPE(T_Spline),PUBLIC:: &
  & Rho_Spl, &
  & Eps_Spl, &
  & dhA_Spl, &
  & dhB_Spl, &
  & bDot_Spl
  !----------------------------------------------------/ Solmodel_Vars --
  
END MODULE M_Solmodel_Vars

