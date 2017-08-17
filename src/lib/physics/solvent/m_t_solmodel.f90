module M_T_SolModel
!--implementation of asymmetric mixing models,
!--water-like, with a solvent and solutes, 
!--and molality-based concentrations
!--(as opposed to symmetric mixing models, implemented with a T_MixModel)

  use M_Kinds
  implicit none

  private
  
  public:: T_SolModel
  public:: T_SolModelDat
  public:: vSolModelAct
  
  public:: SolModel_Spc_Init
  
  type:: T_SolModelDat
    ! container for current values of some properties of a solution phase
    ! density, dielectric, parameters for the activity model
    ! here, Rho is the density of pure solvent, not the solution !!
    real(dp):: Rho, Eps, dhA, dhB, Bdot
  end type T_SolModelDat
  
  type:: T_SolModel  
    !--- describes an assymetric mixing model (SOLUTION model),
    !--- with a solvent and solutes
    !--- (for symmetric mixing models (MIXTURE model), use a T_MixModel)
    !
    ! CAVEAT PROGRAMMATOR:
    ! whereas T_MixModel is currently a fixed size structure,
    ! with fixed maximum numbers of end members, margules, etc,
    ! T_SolModel is a 'dynamic' structure, with the number of solute 
    ! species determined at run time
    !
    character(len=23):: Name  !
    character(len=3) :: Typ   !"LIQ"|"MIN"|"GAS"|...
    !
    integer:: iActModel
    != index of the activity model in vSolModelAct
    !
    integer :: iSolvent
    !iSolvent= index of the solvent species in current species list
    integer :: isH_,isOH,isO2
    real(dp):: MolWeitSv
    != molar weight solvent
    integer:: nSolute
    != number of solute species, excluding solvent = dim' of vISolute
    integer,pointer:: vISolute(:)
    != indices of solute species involved,
    ! vISolute(i) is index of species i in current species list
    ! (= vSpc in M_Global_Vars)    
    logical,pointer:: vIsNonPolar(:)
    ! for EQ3 option: apply activ'coeff' on neutral species
    ! cf EQ36, v7.0, user's manual (Wolery, 1992)
    ! activ'coeff' of CO2(aq) in a (otherwise pure) NaCl solution
    ! as a function of ionic strength
    ! uses an expression proposed by Drummond, 1981
    ! this should be applied only to non-polar species
    ! (O2(aq), H2(aq), N2(aq), ...),
    ! for which a salting out effect is expected
    !
    type(T_SolModelDat):: Dat
    ! current values of some properties of the solution phase,
    ! density, dielectric, DH parameters ...
    
  end type T_SolModel
  
  ! vSolModelAct= 
  ! list of name codes for the different activity models
  ! implemented for assymetric (molality-based) solutions
  integer,parameter:: nSolModelAct= 12
  character(7),parameter:: vSolModelAct(1:nSolModelAct) = &
  & (/ &
  & "IDEAL  ",   & ! 1
  & "DH1    ",   & ! 2
  & "DH1EQ3 ",   & ! 3
  & "DH2    ",   & ! 4
  & "DH2EQ3 ",   & ! 5
  & "DAV_1  ",   & ! 6
  & "DAV_2  ",   & ! 7
  & "PITZER ",   & ! 8
  & "SAMSON ",   & ! 9
  & "HKF81  ",   & ! 10
  & "SIT    ",   & ! 11
  & "name12 "    & ! 12
  & /)
  
contains

subroutine SolModel_Spc_Init(SolvName,vSpc,    S,Ok,OkMsg)
!-- initialize S%iSolvent, S%nSolute, S%vISolute, S%MolWeitSv, etc
!-- according to a species database vSpc,
!-- and given the name, SolvName, of the solvent species

  use M_Trace
  use M_T_Species,only: T_Species,Species_Index

  character(*),    intent(in)   :: SolvName !name of solvent species
  type(T_Species), intent(in)   :: vSpc(:)
  type(T_SolModel),intent(out)  :: S
  logical,         intent(out)  :: Ok
  character(*),    intent(out)  :: OkMsg

  integer:: N,I

  if(idebug>1) write(fTrc,'(/,A)') "< SolModel_Spc_Init"

  Ok= .true.
  OkMsg= ""
  S%iSolvent= 0

  N= 0
  do I=1,size(vSpc)
    if(trim(SolvName)==trim(vSpc(I)%NamSp)) then
      ! solvent
      S%iSolvent=  I
      S%MolWeitSv= vSpc(S%iSolvent)%WeitKg
    else
      ! solutes
      if(vSpc(I)%Typ=="AQU") N= N+1
    end if
  end do

  if(N==0) then
    Ok= .false.
    OkMsg=                              "!! Solute species not found !!"
    return !------------------------------------------------------------
  end if
  !
  if(S%iSolvent==0) then
    Ok= .false.
    OkMsg=                             "!! Solvent species not found !!"
    return !------------------------------------------------------------
  end if
  !
  S%nSolute= N
  allocate(S%vISolute(N))
  !
  N= 0
  do I=1,size(vSpc)
    if(vSpc(I)%Typ=="AQU" .and. trim(SolvName)/=trim(vSpc(I)%NamSp)) then
      N= N+1
      S%vISolute(N)= I
    end if
  end do
  !
  ! find indexes of solvent, H+, OH-, O2 in vSpc
  ! -> unchanged in basis switch -> called only once for a given vSpc
  S%isOH= Species_Index("OH-",  vSpc)
  S%isH_= Species_Index("H+",   vSpc)
  S%isO2= Species_Index("O2_AQ",vSpc)
  !
  if(S%isH_==0) call Stop_("species H+ not found  !!!") !-----------stop
  if(S%isOH==0) call Stop_("species OH- not found !!!") !-----------stop
  !
  ! print *,"<SolModel_Spc_Init"
  ! print *,S%iSolvent,vSpc(S%iSolvent)%Name
  ! print *,S%vISolute(1),"...",S%vISolute(S%nSolute)
  ! print *,"</SolModel_Spc_Init"
  ! pause
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ SolModel_Spc_Init"
  !
end subroutine SolModel_Spc_Init

end module M_T_SolModel
