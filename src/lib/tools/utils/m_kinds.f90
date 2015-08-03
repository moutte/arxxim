!.parameters for size of real types
MODULE M_Kinds !kind defining the size of variables
  INTEGER, PARAMETER :: dp= KIND(1.0D0) 
  INTEGER, PARAMETER :: sp= KIND(1.0)
  REAL(dp),PARAMETER :: Zero=0.0_dp, One=1.0_dp, Two=2.0_dp
ENDMODULE M_Kinds
