module M_Solmodel_Calc_Spycher
!---------------------------------------------------------------------
! Spycher Fugacity Model 
! empirical relation from Spycher & Reed,'88, as implemented by TOUGH2
!---------------------------------------------------------------------

use M_Kinds

public:: SolGas_GammaSpycher

contains

real(dp) function SolGas_GammaSpycher(TdgK,Pbar) !-> ln(fugacity)
  implicit none
  real(dp),intent(in):: TdgK,Pbar
  !
  real(dp):: T,P
  real(dp),parameter::&
  !coeffs fitted for t=50-350°C, P<=500 bars
  & a=-1430.87_dp, b=3.598_dp,    c=-2.27376D-3, &
  & d=3.47644_dp,  e=-1.04247D-2, f=8.36271D-6
  !
  T=TdgK
  P=Pbar
  SolGas_GammaSpycher= &
  & (a/T/T + b/T + c)*P + (d/T/T + e/T + f)*P*P/Two

end function SolGas_GammaSpycher

end module M_Solmodel_Calc_Spycher
