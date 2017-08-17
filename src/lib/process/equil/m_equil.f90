module M_Equil
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_

  implicit none

  private

  public:: Equil_Calc
  public:: Equil_Tst
  public:: Equil_Get_VMolF

contains

subroutine Equil_Calc(Cod)
!--
!-- equilibrium calculation,
!-- switch between the different modes:
!--   SPC>-> speciation: internal equilibrium of the aqueous phase
!--   EQn>-> global equilibrium,
!--          using either "reaction advancement" or direct strategy, or both
!--
  use M_Global_Vars, only: vSpcDat
  use M_System_Vars, only: vCpn,TdgK,Pbar
  !
  use M_Basis,       only: Basis_Change
  use M_Basis,       only: Basis_Change_Wrk
  use M_Equil_Vars,  only: Equil_Vars_Clean,LogForAqu
  !
  use M_Equil_Tools
  use M_Equil_Write
  use M_Equil_Specia
  use M_Equil_1
  use M_Equil_2
  !
  character(len=3),intent(in):: Cod
  !
  integer:: iErr
  ! logical:: Singular
  !
  if(idebug>1) flush(fTrc)
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Calc"
  !
  !Singular= .false.
  !
  call Equil_Write_LogK(vCpn,TdgK,Pbar)
  !
  call Basis_Change("SPC",vCpn)
  !
  call Equil_Zero(Cod)
  !
  call Equil_Trace_Init
  !
  call Equil_Restart
  !
  iErr= 0
  !
  select case(Cod)

  case("SPC","BOX","INJ")
  ! calculate "homogeneous" equilibrium speciation of the aqueous phase
  ! without consideration of satur'state versus other phases
    call Equil_Specia(iErr)

  case("EQM")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using hybrid method 2 and taking mixtures in account

    ! first compute with pure phase only
    call Equil_Eq2(.true.,iErr)
    ! then take account of mixtures
    call Equil_Eq2(.false.,iErr)
    !! if(iErr<0) call Equil_Eq2(.false.,iErr,"2")

  case("EQ1")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using method 1
    call Equil_Eq1(.true.,iErr)

  case("EQ2")
    ! calculate "global" equilibrium of (aqueous + rock) system
    ! considering satur'state of aqu'phase versus other phases
    ! using method 2
    call Equil_Eq2(.true.,iErr)
    if(iErr<0) call Equil_Eq2(.true.,iErr,"2")

    ! if error do again with method 1
    if(iErr<0) then
      if(iDebug>2) then
        print *,"EQ2 FAILED, TRY WITH EQ1"
        call Pause_
      end if
      call Equil_Eq1(.true.,iErr)
    end if

  end select
  !
  if(iErr/=0) call Equil_Errors(iErr)
  !
  !! if(iDebug>2) call Basis_Change_Wrk(vCpn)
  !
  call Equil_Trace_Close
  !
  call Equil_Save
  call Equil_Vars_Clean
  !
  !--- screen output
  if(idebug>1) call Equil_Write_ShoDetail(Cod,TdgK,Pbar,vCpn,vSpcDat)
  !
  !--- file output
  call Equil_Write_Detail(Cod,TdgK,Pbar,vCpn)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Calc"
end subroutine Equil_Calc

subroutine Equil_Spl
!calculate stable assemblage excluding solute species
  !! use M_Global_Vars
  !! use M_System_Vars, only: vCpn
  !! use M_Simplex_Vars,only: tStoikSpl,tSimplex,IZROV,iPosV,Simplex_Vars_Alloc,Simplex_Vars_Clean
  !! use M_Numeric_Simplex
  !! !
  !! integer:: iError,nC,nF
  !! integer:: i
  !! !
  !! write(fTrc,'(/,A,/)') "Equil_Spl, init"
  !! do i= 1, size(vCpn)
  !! !transform the current component set for use with simplex
    !! print *,vCpn(i)%Name
  !! end do
  !! call Pause_
  !! return
  !! nC= size(vCpn)
  !! nF= size(vFas)
  !! !
  !! call Simplex_Vars_Alloc(nC,nF)
  !! !
  !! tSimplex(1:nC,0   )=  vCpn(1:nC)%Mole !first column= bulk compos'n
  !! tSimplex(1:nC,1:nF)= -transpose(tStoikSpl(1:nF,1:nC)) !main = stoikiometry matrix
  !! tSimplex(0,   1:nF)= -vFas(1:nF)%Grt  !first row= Gibbs energy of phases 1:nF at 'point' iTP
  !! !
  !! ! call Simplex_Calc(iError)
  !! call Simplex_Calc_( &
  !! & nC,nF, &
  !! & tSimplex,IZROV,IPOSV, &
  !! & iError) !,n1,n2)
  !
end subroutine Equil_Spl

subroutine Equil_Tst
  use M_Basis,       only: Basis_Change
  use M_Equil_Vars,  only: Equil_Vars_Clean
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_System_Tools,only: System_TP_Update

  ! update global thermo'parameters (vSpc,vMixModel,vMixFas,vFas)
  ! at (TdgK,Pbar)
  call System_TP_Update(TdgK,Pbar)

  call Basis_Change("SPC",vCpn) ; if(iDebug==5) call Pause_("Basis_Change Done")

  call Equil_Vars_Clean

  return
end subroutine Equil_Tst

subroutine Equil_Get_vMolF(vSpcDat,vOrdAq,vX)
  use M_T_Species,  only: T_SpcData

  type(T_SpcData), intent(in) :: vSpcDat(:)
  integer,         intent(in) :: vOrdAq(:)
  real(dp),        intent(out):: vX(:)

  integer :: nAq, iAq

  vX=  Zero
  nAq= size(vOrdAq)

  do iAq=1,nAq
    vX(iAq)= vSpcDat(vOrdAq(iAq))%Mole
  end do

  return
end subroutine Equil_Get_vMolF

end module M_Equil

