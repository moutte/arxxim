subroutine test_EosIfp
!=====================================================
! Small Test for EosIfp Module
! Read Data in a File
!=====================================================
use m_Kinds,only: dp
use m_T_EosIfp_Const
use m_T_EosIfp_Param
use m_T_EosIfp_vState
implicit none

integer, parameter :: LVTH = 2
integer, parameter :: LVTH2 = 1

integer, parameter :: NCONS = 2
character(len=10) :: FileInn 
character(len=10) :: FileInn2 

integer :: i

type (T_EosIfp_Const)  :: Tconst
type (T_EosIfp_Param)  :: Tparam
type (T_EosIfp_Param)  :: Tparam2
type (T_EosIfp_vState) :: Tvstate

real(dp) :: vX(LVTH, NCONS), vT(LVTH), vP(LVTH)

write(*,*) "----------------- Prog Test EOS ifP"
FileInn="case1"
FileInn2="case2"

call EosIfp_const_init(Tconst)

!// TEST 1 =================================================================
!// set of parameters
call EosIfp_Param_CreateFromFile( Tparam, Tconst, FileInn, NCONS )
!!call EosIfp_Param_Print( Tparam, 6)
call EosIfp_vState_new( Tvstate, NCONS, LVTH )
vX (1:LVTH,1) = 0.5
vX (1:LVTH,2) = 0.5
do i=1,LVTH
  vT (i) = 300.D0 + 10.d0*i 
  VP (i) = 10.d5  
end do
call EosIfp_vState_set( Tvstate, vT, vP, vX )
call EosIfp_vState_Computefugacity(  Tvstate, Tconst, Tparam, NCONS, LVTH)
call EosIfp_vState_Print( Tvstate , 6)

!// TEST 2 =================================================================

!// new set of parameters + new vState
call EosIfp_Param_CreateFromFile( Tparam2, Tconst, FileInn2, NCONS )
!!call EosIfp_Param_Print( Tparam, 6)
call EosIfp_vState_new( Tvstate, NCONS, LVTH2  )
vX (1:LVTH2,1) = 0.5
vX (1:LVTH2,2) = 0.5
do i=1,LVTH2
  vT (i) = 300.D0 + 10.d0*i 
  VP (i) = 100.d5  
end do
call EosIfp_vState_set( Tvstate, vT, vP, vX )
call EosIfp_vState_Computefugacity(  Tvstate, Tconst, Tparam, NCONS, LVTH2)
call EosIfp_vState_Print( Tvstate , 6)

end subroutine test_EosIfp
