module M_T_Component
  use M_Kinds
  use M_T_Species,only: nElMax
  ! 
  implicit none 
  private 
  public:: T_Component 
  public:: Component_Zero
  public:: Component_Index
  public:: Component_Print
  public:: Component_Stoikio
  public:: CpnMolMinim
  !
  type:: T_Component !container for storing a thermodynamic component
    character(len=23):: NamCp  ="ZZZ"
    character(len=71):: Formula="ZZZ"     !-> currently used only for simplex !!!
    character(len=7) :: Statut ="INERT" !status of component; INERT,MOBILE,BALANCE,BUFFER,...
    integer :: vStoikCp(0:nElMax)=0
    integer :: iEle= 0
    integer :: iSpc= 0
    real(dp):: Factor= 1.0_dp
    !
    character(len=23):: namMix="ZZZ"
    ! namSol= name of phase controlling this component (used for mobile cpn') 
    !   for an aqueous species or a pure species, namSol="ZZZ" 
    !   for a non aqueous mixture, %namMix is the name of the solid or gas mixture phase 
    !iMix, iPol = used for multi-solution system (e.g. SS-AS)
    !-> iFas is calculated from %namMix, c%iFas=MixPhase_Index(IN=vMixFas,IN=c%namMix) 
    integer :: iMix=0 !index of "controlling" mixture phase in vMixPhase (iFas=0 for aqu. or pure species)
    integer :: iPol=0 !index of the species in the end-member list of that phase 
    ! 
    real(dp):: Mole=  Zero !stores the current mole number (for inert component) 
    real(dp):: LnAct= Zero !stores the current activity (for mobile component) 
    ! 
  end type T_Component
  !
  real(dp):: CpnMolMinim= 1.0D-12 !1.D-16
 
contains 

subroutine Component_Zero(C)
  type(T_Component),intent(out):: C
  !
  C%NamCp=   "ZZZ"
  C%Formula= "ZZZ"
  C%Statut=  "INERT"
  C%vStoikCp(:)= 0
  C%vStoikCp(0)= 1 ! divider
  C%iEle=     0
  C%iSpc=     0
  C%Factor=   1.0_dp
  C%namMix=  "ZZZ"
  C%iMix=     0
  C%iPol=     0
  C%Mole=     Zero
  C%LnAct=    Zero
  !
end subroutine Component_Zero
 
integer function Component_Index(Str,V)
!--
!-- finds component index according to its %NamCp
!--
  character(len=*), intent(in):: Str
  type(T_Component),intent(in):: V(:)
  !
  integer::I 
  !
  Component_Index=0
  if(size(V)==0) return
  !
  I=0
  do 
    I=I+1 
    if(trim(Str)==trim(V(I)%NamCp)) then
      Component_Index=I  ;  exit
    end if 
    if(I==size(V)) exit 
  end do !if Str not found -> I=0
  !
  return
end function Component_Index
 
subroutine Component_Print(f,vEle,vSpc,C) 
  use M_T_Element,only: T_Element 
  use M_T_Species,only: T_Species 
  ! 
  integer,          intent(in):: f 
  type(T_Element),  intent(in):: vEle(:) 
  type(T_Species),  intent(in):: vSpc(:) 
  type(T_Component),intent(in):: C 
  
  write(f,'(3(A,1X),2(A15,1X),2G12.3)') & 
  &  C%NamCp, vEle(C%iEle)%NamEl,C%Statut, &
  &  vSpc(C%iSpc)%NamSp, &
  &  C%namMix, &
  &  C%Mole,C%LnAct
  
  return
end subroutine Component_Print 
 
subroutine Component_Stoikio(vEle,ieOx,Cpn,fOk)
!-- !!! currently used only in simplex routines !!!!!!
!-- build stoichiometry vector of component Cpn 
!-- according to element base vEle
!-- NEW: stoichiometry saved as S%vStoikCp
!--
  use M_T_Element,only: T_Element,Formula_Read
  !---------------------------------------------------------------------
  type(T_Element),  intent(in)::    vEle(:)
  integer,          intent(in)::    ieOx
  type(T_Component),intent(inout):: Cpn
  logical,          intent(out)::   fOk
  !---------------------------------------------------------------------
  integer:: vStoik(1:size(vEle))
  integer:: ZSp,nDiv,nEl,ZExc
  !
  ! print *,"=Component_Stoikio=",trim(Cpn%Formula)  ;  pause
  !
  nEl=  size(vEle)
  !
  call Formula_Read(  &
  & Cpn%Formula,vEle, &
  & ZSp,nDiv,fOk,vStoik)
  !
  Cpn%vStoikCp(:)=     0
  Cpn%vStoikCp(1:nEl)= vStoik(1:nEl)
  Cpn%vStoikCp(0)=     nDiv
  Cpn%vStoikCp(nEl+1)= ZSp  !-> component charge
  !
  return
  !
  ZExc=dot_product(vStoik(1:nEl),vEle(1:nEl)%Z) !-> oxidation state
  !
  if(ieOx/=0) then
    Cpn%vStoikCp(ieOx)= ZSp - ZExc
    fOk=.true.
  end if
  !
  return
end subroutine Component_Stoikio

end module M_T_Component
 
