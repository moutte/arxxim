MODULE M_Equil
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Equil_Calc
  PUBLIC:: Equil_Tst
  PUBLIC:: Equil_Get_VMolF

CONTAINS

SUBROUTINE Equil_Calc(Cod)
!--
!-- equilibrium calculation,
!-- switch between the different modes:
!--   SPC>-> speciation: internal equilibrium of the aqueous phase
!--   EQn>-> global equilibrium,
!--          using either "reaction advancement" or direct strategy, or both
!--
  USE M_System_Vars, ONLY: vCpn,TdgK,Pbar
  !
  USE M_Basis,       ONLY: Basis_Change
  USE M_Basis,       ONLY: Basis_Change_Wrk
  USE M_Equil_Vars,  ONLY: Equil_Vars_Clean,LogForAqu
  !
  USE M_Equil_Tools
  USE M_Equil_Write
  USE M_Equil_Specia
  USE M_Equil_1
  USE M_Equil_2
  !
  CHARACTER(LEN=3),INTENT(IN):: Cod
  !
  INTEGER:: iErr
  ! LOGICAL:: Singular
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Calc"
  !
  !Singular= .false.
  !
  CALL Equil_Write_LogK(vCpn,TdgK,Pbar)
  !
  CALL Basis_Change("SPC",vCpn)
  !
  CALL Equil_Zero(Cod)
  !
  CALL Equil_Trace_Init
  !
  CALL Equil_Restart
  !
  iErr= 0
  !
  SELECT CASE(Cod)

  CASE("SPC","BOX","INJ")
  ! calculate "homogeneous" equilibrium speciation of the aqueous phase
  ! without consideration of satur'state versus other phases
    CALL Equil_Specia(iErr)

  CASE("EQM")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using hybrid method 2 and taking mixtures in account

    ! first compute with pure phase only
    CALL Equil_Eq2(.TRUE.,iErr)
    ! then take account of mixtures
    CALL Equil_Eq2(.FALSE.,iErr)
    !~ IF(iErr<0) CALL Equil_Eq2(.FALSE.,iErr,"2")

  CASE("EQ1")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using method 1
    CALL Equil_Eq1(.TRUE.,iErr)

  CASE("EQ2")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using method 2
    CALL Equil_Eq2(.TRUE.,iErr)
    IF(iErr<0) CALL Equil_Eq2(.TRUE.,iErr,"2")

    ! if error do again with method 1
    IF(iErr<0) THEN
      IF(iDebug>2) THEN
        PRINT *,"EQ2 FAILED, TRY WITH EQ1"
        CALL Pause_
      ENDIF
      CALL Equil_Eq1(.TRUE.,iErr)
    ENDIF

  ENDSELECT
  !
  IF(iErr/=0) CALL Equil_Errors(iErr)
  !
  !~ IF(iDebug>2) CALL Basis_Change_Wrk(vCpn)
  !
  CALL Equil_Trace_Close
  !
  CALL Equil_Save
  CALL Equil_Vars_Clean
  !
  !--- screen output
  IF(iDebug>0) CALL Equil_Write_ShoDetail(Cod,TdgK,Pbar,vCpn)
  !
  !--- file output
  CALL Equil_Write_Detail(Cod,TdgK,Pbar,vCpn)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Calc"
ENDSUBROUTINE Equil_Calc

SUBROUTINE Equil_Spl
!calculate stable assemblage excluding solute species
  !~ USE M_Global_Vars
  !~ USE M_System_Vars, ONLY: vCpn
  !~ USE M_Simplex_Vars,ONLY: tStoikSpl,tSimplex,IZROV,iPosV,Simplex_Vars_Alloc,Simplex_Vars_Clean
  !~ USE M_Numeric_Simplex
  !~ !
  !~ INTEGER:: iError,nC,nF
  !~ INTEGER:: i
  !~ !
  !~ WRITE(fTrc,'(/,A,/)') "Equil_Spl, init"
  !~ DO i= 1, SIZE(vCpn)
  !~ !transform the current component set for USE with simplex
    !~ PRINT *,vCpn(i)%Name
  !~ ENDDO
  !~ CALL Pause_
  !~ RETURN
  !~ nC= SIZE(vCpn)
  !~ nF= SIZE(vFas)
  !~ !
  !~ CALL Simplex_Vars_Alloc(nC,nF)
  !~ !
  !~ tSimplex(1:nC,0   )=  vCpn(1:nC)%Mole !first column= bulk compos'n
  !~ tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikSpl(1:nF,1:nC)) !main = stoikiometry matrix
  !~ tSimplex(0,   1:nF)= -vFas(1:nF)%Grt  !first row= Gibbs energy of phases 1:nF at 'point' iTP
  !~ !
  !~ ! CALL Simplex_Calc(iError)
  !~ CALL Simplex_Calc_( &
  !~ & nC,nF, &
  !~ & tSimplex,IZROV,IPOSV, &
  !~ & iError) !,n1,n2)
  !
ENDSUBROUTINE Equil_Spl

SUBROUTINE Equil_Tst
  USE M_Basis,       ONLY: Basis_Change
  USE M_Equil_Vars,  ONLY: Equil_Vars_Clean
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_System_Tools,ONLY: System_TP_Update

  ! update global thermo'parameters (vSpc,vMixModel,vMixFas,vFas)
  ! at (TdgK,Pbar)
  CALL System_TP_Update(TdgK,Pbar)

  CALL Basis_Change("SPC",vCpn) ; IF(iDebug==5) CALL Pause_("Basis_Change Done")

  CALL Equil_Vars_Clean

  RETURN
ENDSUBROUTINE Equil_Tst

SUBROUTINE Equil_Get_vMolF(vX)
  USE M_Global_Vars, ONLY : vSpc
  USE M_Basis_Vars, ONLY : vOrdAq

  REAL(dp), INTENT(OUT)::vX(:)

  INTEGER :: nAq, iAq

  vX = Zero
  nAq = SIZE(vOrdAq)

  DO iAq=1, nAq
    vX(iAq) =  vSpc(vOrdAq(iAq))%Dat%Mole
  END DO

  RETURN
END SUBROUTINE Equil_Get_vMolF

ENDMODULE M_Equil

