MODULE M_SolGas_Spycher
!---------------------------------------------------------------------
! Spycher Fugacity Model 
! empirical relation from Spycher & Reed,'88, as implemented by TOUGH2
!---------------------------------------------------------------------

USE M_Kinds

PUBLIC::SolGas_GammaSpycher

CONTAINS

REAL(dp) FUNCTION SolGas_GammaSpycher(TdgK,Pbar) !-> ln(fugacity)
  IMPLICIT NONE
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  REAL(dp):: T,P
  REAL(dp),PARAMETER::&
  !coeffs fitted for t=50-350°C, P<=500 bars
  & a=-1430.87_dp, b=3.598_dp,    c=-2.27376D-3, &
  & d=3.47644_dp,  e=-1.04247D-2, f=8.36271D-6
  !
  T=TdgK
  P=Pbar
  SolGas_GammaSpycher= &
  & (a/T/T + b/T + c)*P + (d/T/T + e/T + f)*P*P/Two

ENDFUNCTION SolGas_GammaSpycher

ENDMODULE M_SolGas_Spycher
