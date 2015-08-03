MODULE M_Formula_Parser

  !------------------------------------------------------
  ! Interface to Formula Parsing Utils 
  !------------------------------------------------------

  USE M_Formula_Vector
  USE M_Formula_Arxim
  USE M_Formula_Arxim_Standard
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Formula_Vector_Read
  !---
  PUBLIC :: Formula_Arxim_Build
  PUBLIC :: Formula_Arxim_Read
  !---
  PUBLIC :: Formula_Arxim_Read_Standard
  PUBLIC :: Formula_Arxim_Build_Standard

ENDMODULE M_Formula_Parser

  
