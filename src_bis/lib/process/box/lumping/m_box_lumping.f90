module M_Box_Lumping

  use M_Kinds
  use M_Trace

  use M_Box_Lumping_Vars

  implicit none
  private
  public :: Box_Lumping_Read
  public :: Box_Lumping_Output_Simres
  public :: Box_Lumping_Output_Coupler
  public :: Box_Lumping_Output_LumpCompo
  public :: Box_Lumping_Print

contains

subroutine Box_Lumping_Read(Ok)

!=========================================================
! Read Box Lumping Parameters
! BOX.LUMPING
!=========================================================

  use M_Files,only: NamFInn

  use M_IOTools  
  use M_SetOption

  implicit none

  !-- arguments
  logical,       intent(out)  :: Ok

  !-- local variables
  character(len=255):: L
  character(len=80) :: W, W1, W2, W3, W4
  character(len=80) :: cPhase, cCons, cSpecies
  logical           :: EoL
  integer           :: F,ios

  !-------------------------------------------------------
  call Info_("Box_Read_Lumping")

  Ok=.false.

  !// Default Values
  call Info_("Set Default Values")

  call SetOption( "Option Lumping",  .false.", LOption_Lumping )

  !call SetOption( "Lumping Phase",   "WATER ", cPhase )
  !call SetOption( "Lumping Cons",    "H2O   ", cCons )
  !call SetOption( "Lumping Species", "H2O   ", cSpecies )
  !call Box_Lumping_Add( cPhase,cCons, cSpecies) 

  call Info_("Read Bloc")
  !// Open Input File
  call GetUnit(F)
  open(F,file=trim(NamFInn))

  Do01: do 
  read(F,'(A)',iostat=ios) L
  if(ios/=0) exit Do01

  call LinToWrd(L,W,EoL)
  call AppendToEnd(L,W,EoL)

  if(W(1:1)=='!') cycle Do01 !skip comment lines

  select case(trim(W))

  case("ENDINPUT")
  exit  Do01

  case("BOX.LUMPING")
  Ok=.true.
  call SetOption( "Option Lumping",  .true.", LOption_Lumping )          
  Do02: do
  read(F,'(A)',iostat=ios) L
  if(ios/=0) exit Do01

  call LinToWrd(L,W1,EoL)
  call AppendToEnd(L,W1,EoL)

  if(W1(1:1)=='!') cycle Do02

  call LinToWrd(L,W2,EoL)
  call LinToWrd(L,W3,EoL)
  call LinToWrd(L,W4,EoL)

  !-> W2 is second word on line
  !-> W3 is third word on line
  !-> W4 is fourth word on line

  select case(trim(W1))

  case("ENDINPUT")
  exit Do01

  case("END","ENDBOX.LUMPING")
  exit Do02

  case("CONS") ;
  call SetOption( "Lumping Phase",   W2, cPhase )
  call SetOption( "Lumping Cons",    W3, cCons )
  call SetOption( "Lumping Species", W4, cSpecies )
  call Box_Lumping_Add( cPhase, cCons, cSpecies)

  case default
  call Fatal_("BOX.LUMPING = "//trim(W1)//" unknown keyword")

  end select

  enddo Do02

  end select
  enddo Do01

  !// Close Input File
  close(F)

  call Info_("Box_Lumping_Read Done")

  call Box_Lumping_Add_CLump
  call Box_Lumping_Print(6)

end subroutine Box_Lumping_Read

!---

subroutine Box_Lumping_Add( cPhase,cCons, cSpecies)
  use M_Global_Vars, only : vSpc
  implicit none
  character(len=*) :: cPhase, cCons, cSpecies
  !---
  call Info_("Box_Lumping_Add")
  call Lumping_Map_Add( cPhase, cCons, cSpecies, vLumpMap, Size_LumpMap, vSpc) 

  end subroutine Box_Lumping_Add

  !---

  subroutine Box_Lumping_Add_CLump
  use M_StringUtils
  use M_Global_Vars, only : vSpc
  use M_Box_System_Vars
  use M_SetOption
  implicit none

  integer :: iMap, nMap
  integer :: iAq, iSpc
  character(len=40) :: SpcName
  logical :: LPhaseWater
  character(len=80) :: cPhase, cCons, cSpecies
  !---

  call Info_("Box_Lumping_Add_Residual")

  do iAq = 1, nAq

  iSpc = vOrdAq(iAq)
  SpcName = trim(vSpc(iSpc)%NamSp) 

  call Find_Lumping_Map_Species(iSpc, iMap, nMap, vLumpMap, Size_LumpMap)
  if (nMap>1) call Fatal_("Multiple Species Definition For "// trim(SpcName) ) 

  if (iMap>0) then
  LPhaseWater = Lumping_Map_IsPhaseWater(vLumpMap(iMap))
  if ( LPhaseWater ) then
  ! nothing to do Ok
  else
  call Fatal_("Aqueous Species "// trim(SpcName) &
  //" Lumped in Phase "// trim(vLumpMap(iMap)%cPhase)  )
  end if
  else
  call SetOption( "Lumping Phase",   'WATER  ', cPhase   )
  call SetOption( "Lumping Cons",    'C-LUMP ', cCons    )
  call SetOption( "Lumping Species", SpcName,   cSpecies )
  call Box_Lumping_Add( cPhase, cCons, cSpecies)
  end if
  end do

  iCons_CLUMP = String_Index(vLumpMap(1:Size_LumpMap)%cCons , 'C-LUMP')
  iCons_H2O   = String_Index(vLumpMap(1:Size_LumpMap)%cCons , 'H2O')
  call Box_Lumping_Set_CLump_OutPutName

end subroutine Box_Lumping_Add_CLump

!---

subroutine Box_Lumping_Set_CLump_OutPutName
  use M_StringUtils
  integer :: idxCONS_H2O
  if (iCons_H2O==0) then
  CLump_OutputName = "H2O" 
  else
  CLump_OutputName = "C-LUMPPPP"
  end if
end subroutine Box_Lumping_Set_CLump_OutPutName

!---

subroutine Box_Lumping_Print(F)
  use M_Box_Debug_Vars, only : LDebug_Lumping
  implicit none
  integer :: iMap
  integer :: F
  !---
  call Info_("Box_Lumping_Print")

  if ( LDebug_Lumping ) then
    call Lumping_Map_Print_Header(F)
    do iMap = 1, Size_LumpMap
    call Lumping_Map_Print(vLumpMap(iMap), iMap, F)
    end do
    write(*,*) 
  end if

end subroutine Box_Lumping_Print

!---

subroutine Box_Lumping_Output_Coupler(F)
  use M_Box_Debug_Vars, only : LDebug_Lumping
  implicit none
  integer, intent(in) :: F
  !---
  call Box_Lumping_Output_Warning_H2OLump(F)
  call Box_Lumping_Output_Structure_Coupler(F)

end subroutine Box_Lumping_Output_Coupler

!---

subroutine Box_Lumping_Output_Simres(F)
  implicit none
  integer, intent(in) :: F
  !---
  call Box_Lumping_Output_Warning_H2OLump(F)
  call Box_Lumping_Output_Structure_Simres(F)   
  call Box_Lumping_Output_Composition_Simres(F)

end subroutine Box_Lumping_Output_Simres

!---

subroutine Box_Lumping_Output_LumpCompo(F)
  implicit none
  integer, intent(in) :: F
  !---
  call Box_Lumping_Output_Composition_Info(F)

end subroutine Box_Lumping_Output_LumpCompo



subroutine Box_Lumping_Output_Warning_H2OLump(F)
  implicit none
  integer, intent(in) :: F
  if ( trim(CLump_OutputName) == "H2O" ) then

  write(F,'(A)') ' '
  write(F,'(A)') '<<###########################################'
  write(F,'(A)') '<< WARNING : LUMPING MODEL C-LUMP = H2O      '
  write(F,'(A)') '<<###########################################'
  write(F,'(A)') ' '
end if

end subroutine Box_Lumping_Output_Warning_H2OLump
!---

subroutine Box_Lumping_Output_Structure_Simres(F)
  use M_Box_System_Vars
  use M_Global_Vars, only : vSpc, vEle
  use M_Basis_Vars,  only : isW, MWSv
  use M_string
  implicit none
  integer, intent(in) :: F
  !---
  real(dp), parameter :: minalfa = 1.D-6
  !---
  type(T_Lumping_Map)  :: M
  real(dp) :: MolarMass_CONS
  character(len=20) :: Name_CONS 
  real(dp) :: Coef_EAQtoCONS 
  integer  :: iMap, iSpecies
  integer  :: iAq, iBasis, iEle, iCp, iSpc
  real(dp) :: MolarMass_CLUMP
  integer  :: iSpecies_CLUMP
  logical  :: LPhaseWater
  REAl(dp) :: alfa
  integer :: countCLUMP
  logical :: OkView

  !--
  call Info_("Box_Lumping_Output_Structure_Simres")    
  !---    

  !// C-LUMP H2O Equivalence 
  MolarMass_CLUMP = MWSv
  iSpecies_CLUMP  = iSW

  !// MOLWGT CONS
  write(F,'(A)') '<<-------------------------------------------------------------'
  write(F,'(A)') '<< COORES PVT LUMPED MODEL'
  write(F,'(A)') '<<-------------------------------------------------------------'

  write(F,'(A)') ' '
  write(F,'(A)') '<<------------ MOLWGT CONS -----------------'
  write(F,'(A)') ' ' 

  do iMap = 1, Size_LumpMap
  M = vLumpMap(iMap)   

  if (M%iCons .eq. iMap ) then  

    ! Compute Molar Mass
    if ( M%iCons .eq. iCons_CLUMP ) then             
      MolarMass_CONS = MolarMass_CLump *1.D3  
      Name_CONS = CLump_OutputName
    else
      MolarMass_CONS = vSpc(M%iSpecies)% WeitKg *1.D3  
      Name_CONS = M%cCons
    end if          

    ! Define Components
    if ( trim(Name_CONS) == "H2O" ) then
      write(F,'((A10,A10, A),G18.11)' ) 'MOLWGT-W  ', '  ' , ' = ', MolarMass_CONS
    else
      write(F,'((A10,A10,A),G18.11)' )  'MOLWGT    ' , trim(Name_CONS), ' = ', MolarMass_CONS
    end if
  end if

  end do


end subroutine Box_Lumping_Output_Structure_Simres

!---

subroutine Box_Lumping_Output_Structure_Coupler(F)
  use M_Box_System_Vars
  use M_Global_Vars, only : vSpc, vEle
  use M_Basis_Vars,  only : isW, MWSv
  use M_string
  implicit none
  integer, intent(in) :: F
  !---
  real(dp), parameter :: minalfa = 1.D-6
  !---
  type(T_Lumping_Map)  :: M
  real(dp) :: MolarMass_CONS
  character(len=20) :: Name_CONS
  real(dp) :: Coef_EAQtoCONS 
  integer  :: iMap, iSpecies
  integer  :: iAq, iBasis, iEle, iCp, iSpc
  real(dp) :: MolarMass_CLUMP
  integer  :: iSpecies_CLUMP
  logical  :: LPhaseWater
  REAl(dp) :: alfa
  integer :: countCLUMP
  logical :: OkView   
  character(len=20) :: NameSpc
  character(len=3)  :: NameEle

  !--
  call Info_("Box_Lumping_Output_Structure_Coupler")    
  !---    

  !// C-LUMP H2O Equivalence 
  MolarMass_CLUMP = MWSv
  iSpecies_CLUMP  = iSW

  !// MOLWGT CONS
  write(F,'(A)') '<<-------------------------------------------------------------'
  write(F,'(A)') '<< COORES COUPLER LUMPING MODEL'
  write(F,'(A)') '<<-------------------------------------------------------------'

  !// CONStoC
  write(F,'(A)') ' '
  write(F,'(A)') '<<------------- CONStoC -----------------'

  countCLUMP = 0

  do iMap = 1, Size_LumpMap
  M = vLumpMap(iMap)      
  LPhaseWater = Lumping_Map_IsPhaseWater(M) 
  OkView = .false.

  if (LPhaseWater) then

  !// check CLUMP
  if (M%iCons .eq. iCons_CLUMP ) then           
  iSpecies = iSpecies_CLUMP
  countCLUMP = countCLUMP +1
  OkView = (countCLUMP <= 1)
  Name_CONS = CLump_OutputName
  else
  iSpecies = M%iSpecies
  OkView = .true.
  Name_CONS = M%cCons
  end if

  !// write result if OkView
  if (Okview) then
  write(F,'(A)') ' '
  !write(*,*) "Looking for Species", vSpc(iSpecies)%NamSp
  iBasis = 0
  do iAq =1, nAq
  iSpc = VordAq(iAq)
  if ( iSpc == iSpecies ) then
  iBasis = iAq
  exit
  end if
  end do       
  iAq = iBasis	     
  if (iAq==0) call fatal_("Wrong Aqueous Species " // trim ( vSpc(M%iSpecies)%NamSp) )

  do iCp = 1, nCp
  alfa = tAlfAq(iCp, iAq)
  iEle = vCpnBox(iCp)%iEle
  NameEle = vEle(iEle)%NamEl
  call string_lowcase ( NameEle )
  call string_subst ( NameEle, '_', ' ')
  if ( abs(alfa) > minalfa ) then
  write(F,'(3(A10,1X),A, G18.11)' ) "CONStoC     ",  &
  trim(Name_CONS),  trim(NameEle), " = ", alfa 
  end if
  end do
  end if
  end if
  end do

  !// EAQtoCONS
  write(F,'(A)') ' '
  write(F,'(A)') '<<------------- EAQtoCONS -----------------'
  write(F,'(A)') ' '

  do iMap = 1, Size_LumpMap
  M = vLumpMap(iMap)      
  LPhaseWater = Lumping_Map_IsPhaseWater(M) 
  if (LPhaseWater) then
  if (M%iCons .eq. iCons_CLUMP ) then
  Coef_EAQtoCONS = vSpc(M%iSpecies)% WeitKg / MolarMass_CLUMP
  Name_CONS = CLump_OutputName
  else
  Coef_EAQtoCONS = One
  Name_CONS =  M%cCons
  end if
  NameSpc = trim(M%cSpecies)
  call string_lowcase(NameSpc)
  write(F,'(A12, A20, A10 ,A, G18.11)' ) "EAQtoCONS    ", NameSpc, Name_CONS, " = ", Coef_EAQtoCONS
  end if

  end do

  write(F,'(A)') ' '
  write(F,'(A)') '<<------ Specific options for arxim -------'
  write(F,'(A)') ' '
  write(F,'(A)') 'OPT_CHIMIE          = 4'

end subroutine Box_Lumping_Output_Structure_Coupler

!---

subroutine Box_Lumping_Output_Composition_Simres(F)
  use M_Box_Vars, only : vMolF
  use M_Box_System_Vars
  use M_Global_Vars, only : vSpc, vEle
  use M_Basis_Vars,  only : isW, MWSv
  use M_string
  implicit none
  !---
  type(T_Lumping_Map)  :: M

  integer  :: iMap, iCp, iEle
  integer  :: iAq, iBasis, iSpc

  real(dp) :: Coef_EAQtoCONS 
  real(dp) :: MolarMass_CLUMP
  character(len=20) :: Name_CONS
  logical  :: LPhaseWater
  real(dp) :: alfa

  integer :: nKons, iKons, iKonsReal
  character(len=10), allocatable :: vNameKons(:)
  REAl(dp),          allocatable :: vMolKons(:)
  logical,           allocatable :: vRealKons(:)
  real(dp) :: MolTotalKons, MassH2O, MassMolH2O

  character(len=20) :: NameSpc    
  integer, intent(in) :: F

  !--
  call Info_("Box_Lumping_Output_Composition_Simres")    
  !---    

  MolarMass_CLUMP = MWSv

  nKons = 0
  do iMap = 1, Size_LumpMap
  nKons= max(vLumpMap(iMap)%iCons, nKons)
  end do

  allocate(vMolKons(nKons))
  allocate(vRealKons(nKons))
  allocate(vNameKons(nKons))

  vMolKons(1:nKons)  = Zero
  vRealKons(1:nKons) = .false.
  vNameKons(1:nKons) = 'none'

  !// COMPO EAQtoCONS
  write(F,'(A)') ' '
  write(F,'(A)') '<<-------------  Composition CONS -----------------'
  write(F,'(A)') ' '    

  do iMap = 1, Size_LumpMap
  M = vLumpMap(iMap)      
  LPhaseWater = Lumping_Map_IsPhaseWater(M) 
  if (LPhaseWater) then
  if (M%iCons .eq. iCons_CLUMP ) then
  Coef_EAQtoCONS = vSpc(M%iSpecies)% WeitKg / MolarMass_CLUMP
  Name_CONS = CLump_OutputName
  else
  Coef_EAQtoCONS = One
  Name_CONS = M%cCons
  end if
  !write(*,*) "Looking for Species", vSpc(M%iSpecies)%NamSp
  iBasis = 0
  do iAq =1, nAq
  iSpc = VordAq(iAq)
  if ( iSpc == M%iSpecies ) then
  iBasis = iAq
  exit
  end if
  end do       
  iAq = iBasis	     
  if (iAq==0) call fatal_("Wrong Aqueous Species " // trim ( vSpc(M%iSpecies)%NamSp) )

  iKons = M%iCons
  vNameKons(iKons) = Name_CONS
  vMolKons(iKons)  = vMolKons(M%iCons) + Coef_EAQtoCONS*vMolF(iAq) 
  vRealKons(iKons) = .true.

  if (iSpc == iSW) then 
  MassMolH2O = 0.D0
  do iCp = 1, nCp
  alfa = tAlfAq(iCp, iAq)
  iEle = vCpnBox(iCp)%iEle                
  MassMolH2O = MassMolH2O + alfa* vEle(iEle)%WeitKg                   
  end do
  MassH2O = MassMolH2O* vMolF(iAq) 

  end if
  end if

  end do

  MolTotalKons = sum(vMolKons(1:nKons))    

  do iKons = 1, nKons        
  if ( vRealKons(iKons) ) &
  write(F,'(A14, A10 ,A, G18.11)' ) "<< MOLAL-K   ", vNameKons(iKons), " = ", vMolKons(iKons)/MassH2O        
  end do  

  write(F,'(A)') ' ' 

  do iKons = 1, nKons        
  if ( vRealKons(iKons) ) &
  write(F,'(A14, A10 ,A, G18.11)' ) "WK-CELL IG:I:J:K ", vNameKons(iKons), " = ", vMolKons(iKons)/MolTotalKons        
  end do

  write(F,'(A)') ' '    
  write(F,'(A14, A10 ,A, G18.11)' ) "COMP-WATER TRAP = "
  do iKons = 1, nKons        
  if ( vRealKons(iKons) ) then
  iKonsReal = iKonsReal +1
  if(trim(vNameKons(iKons))=="H2O") then
  write(F,'(A, G18.11, A, I2, A, A10)' ) '< ', vMolKons(iKons)/MolTotalKons , '< [ ', iKonsReal, ' ] ', vNameKons(iKons) 
  else
  write(F,'(G18.11, A, I2, A, A10)' ) vMolKons(iKons)/MolTotalKons , '< [ ', iKonsReal, ' ] ', vNameKons(iKons) 
  end if
  end if
  end do

  deallocate(vMolKons)
  deallocate(vNameKons)
  deallocate(vRealKons)

end subroutine Box_Lumping_Output_Composition_Simres

!---

subroutine Box_Lumping_Output_Composition_Info(F)
  use M_Box_Vars, only : vMolF
  use M_Box_System_Vars
  use M_Global_Vars, only : vSpc, vEle
  use M_Basis_Vars,  only : isW, MWSv
  use M_string
  implicit none
  !---
  type(T_Lumping_Map)  :: M

  integer  :: iMap, iCp, iEle
  integer  :: iAq, iBasis, iSpc

  real(dp) :: Coef_EAQtoCONS 
  real(dp) :: MolarMass_CLUMP
  character(len=20) :: Name_CONS
  logical  :: LPhaseWater
  real(dp) :: alfa

  integer :: nKons, iKons, iKonsReal
  character(len=10), allocatable :: vNameKons(:)
  REAl(dp),          allocatable :: vMolKons(:)
  logical,           allocatable :: vRealKons(:)
  real(dp) :: MolTotalKons, MassH2O, MassMolH2O

  character(len=20) :: NameSpc    
  integer, intent(in) :: F

  !--
  call Info_("Box_Lumping_Output_Composition_Simres")    
  !---    

  MolarMass_CLUMP = MWSv

  nKons = 0
  do iMap = 1, Size_LumpMap
  nKons= max(vLumpMap(iMap)%iCons, nKons)
  end do

  allocate(vMolKons(nKons))
  allocate(vRealKons(nKons))
  allocate(vNameKons(nKons))

  vMolKons(1:nKons)  = Zero
  vRealKons(1:nKons) = .false.
  vNameKons(1:nKons) = 'none'

  !// COMPO EAQtoCONS   
  do iMap = 1, Size_LumpMap
  M = vLumpMap(iMap)      
  LPhaseWater = Lumping_Map_IsPhaseWater(M) 
  if (LPhaseWater) then
  if (M%iCons .eq. iCons_CLUMP ) then
  Coef_EAQtoCONS = vSpc(M%iSpecies)% WeitKg / MolarMass_CLUMP
  Name_CONS = CLump_OutputName
  else
  Coef_EAQtoCONS = One
  Name_CONS = M%cCons
  end if
  !write(*,*) "Looking for Species", vSpc(M%iSpecies)%NamSp
  iBasis = 0
  do iAq =1, nAq
  iSpc = VordAq(iAq)
  if ( iSpc == M%iSpecies ) then
  iBasis = iAq
  exit
  end if
  end do       
  iAq = iBasis
  if (iAq==0) call fatal_("Wrong Aqueous Species " // trim ( vSpc(M%iSpecies)%NamSp) )

  iKons = M%iCons
  vNameKons(iKons) = Name_CONS
  vMolKons(iKons)  = vMolKons(M%iCons) + Coef_EAQtoCONS*vMolF(iAq) 
  vRealKons(iKons) = .true.

  if (iSpc == iSW) then 
  MassMolH2O = 0.D0
  do iCp = 1, nCp
  alfa = tAlfAq(iCp, iAq)
  iEle = vCpnBox(iCp)%iEle                
  MassMolH2O = MassMolH2O + alfa* vEle(iEle)%WeitKg                   
  end do
  MassH2O = MassMolH2O* vMolF(iAq) 

  end if
  end if

  end do

  MolTotalKons = sum(vMolKons(1:nKons))    
  write(F,'(A)') ' '
  write(F,'(A)') '<<----------- Lumping Composition ------------------------' 

  do iKons = 1, nKons        
  if ( vRealKons(iKons) ) &
  write(F,'(A,A10 ,A, G18.11)' ) "< Molality ", vNameKons(iKons), " = ", vMolKons(iKons)/MassH2O        
  end do  


  write(F,'(A)') '<<---'
  do iKons = 1, nKons        
  if ( vRealKons(iKons) ) &
  write(F,'(A,A10 ,A, G18.11)' ) "< Mole Fraction ", vNameKons(iKons), " = ", vMolKons(iKons)/MolTotalKons        
  end do
  write(F,'(A)') '<<-------------------- ' 

  deallocate(vMolKons)
  deallocate(vNameKons)
  deallocate(vRealKons)

end subroutine Box_Lumping_Output_Composition_Info

end module M_Box_Lumping
