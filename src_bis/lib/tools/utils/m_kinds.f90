!.parameters for size of real types
module M_Kinds !kind defining the size of variables
  integer, parameter :: dp= KIND(1.0D0) 
  integer, parameter :: sp= KIND(1.0)
  real(dp),parameter :: Zero=0.0_dp, One=1.0_dp, Two=2.0_dp
end module M_Kinds
