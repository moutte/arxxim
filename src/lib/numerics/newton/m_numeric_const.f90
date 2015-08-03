MODULE M_Numeric_Const !USEful constants
  USE M_Kinds
  IMPLICIT NONE
  PRIVATE
  REAL(dp),PARAMETER,PUBLIC:: &
    Pi=        3.14159265359979_dp, &
    !Pi43=      4._dp*Pi/3._dp, & !4.18879020479972_dp, & != 4.pi/3
    Ln10=      2.302585092994_dp, & != LOG(10.0_dp)
    !MaxExpDP=  200.0_dp, &      !MAXEXPONENT(rDum), bricolage.. to prevent OVERFLOW with EXP in Newton
    !MinExpDP= -200.0_dp, &      !MINEXPONENT(rDum), bricolage.. to prevent OVERFLOW with EXP in Newton
    MaxExpDP=  300.0_dp, &      !MAXEXPONENT(rDum), bricolage.. to prevent OVERFLOW with EXP in Newton
    MinExpDP= -300.0_dp, &      !MINEXPONENT(rDum), bricolage.. to prevent OVERFLOW with EXP in Newton
    TinyDP=    1.0E-300_dp      !TINY(Pi)
  !TINY(x)->smallest number in model of same TYPE and kind PARAMETERs as x
  !  ex:REAL(8) r -> TINY(r)=2.22...E-308
  !HUGE(x)->largest number in model of same TYPE and kind PARAMETERs as x
  !  ex:REAL(4) r -> HUGE(r)=1...E+128
  !EPSILON
  !='a nearly negligible number relative to 1'
  !EPSILON(x)= the smallest REAL>0, of same subTYPE as x, such that 1+ EPSILON(x) > 1
  !for REAL(dp)::x, Epsilon(x)=2.E-16
  !for REAL(sp)::x, Epsilon(x)=1.E-07
ENDMODULE M_Numeric_Const

!!ovlog - the largest INTEGER such that exp(ovlog) DOes not overflow. -> idem MaxExpDP
!!unlog - the largest INTEGER such that exp(-unlog) DOes not underflow. -> idem -MinExpDP
