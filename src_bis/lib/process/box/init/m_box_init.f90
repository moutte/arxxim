module M_Box_Init

  use M_Kinds
  use M_Trace  
  
  use M_Box_Init_Vars
  use M_Box_System_Vars
  
  implicit none
  private

  !// public members
  public :: BOX_Init_System
  public :: BOX_Init_Reset_Basis

  public :: Box_Init_FluidBox
  public :: BOX_Init_FluidInject

  public :: Box_Init_VMolF_FluidBox
  public :: Box_Init_VMolF_FluidInject
  
  public :: Box_Init_VMolK_Rock
  public :: Box_Init_Texture_Rock
   
  public :: Box_Init_Info

contains

  subroutine Box_Init_Reset_Basis
    !===================================================
    ! Initialize Box System Parameters
    !===================================================
    use M_Box_System    
    implicit none
    call Info_("Box_Init_Reset_Basis")    
    call Box_System_Reset_Basis  
    
  end subroutine BOX_Init_Reset_Basis
  
  !---

  subroutine Box_Init_System
    !===================================================
    ! Initialize Box System Parameters
    !===================================================
    use M_Box_System
    use M_System_Vars, only: vCpn
    use M_Box_System_Vars, only: vCpnBox
    use M_System_Vars, only: TdgK,Pbar
    implicit none
    call Info_("Box_Init_System")
    
    call Box_System_Init
    call Box_Init_Vars_New 
    
    call Store_TPCond(TdgK, Pbar)
    call Store_Composition_Species(vMolSpc_FluidSpc)
    call Store_Composition_Component(vCpn, vCpnBox, vMolCpn_FluidSpc)
    
  end subroutine BOX_Init_System

  !---

  subroutine Box_Init_FluidBox(Ok)
    !===============================================================
    ! Init Box Fluid Composition
    !===============================================================
    ! Step 1. Read the system of constraint on box fluid
    ! Step 2. Compute speciation
    ! Step 3. Store the composition 
    !===============================================================
    use M_Global_Vars, only: vEle,vSpc,vMixFas !! ,vTPcond 
    use M_System_Vars, only: TdgK,Pbar,vCpn,vCpnTot
    use M_Dynam_Vars, only: vCpnInj
    use M_Basis_Vars, only: vOrdAq
    use M_System_Tools,only: System_Build_Custom 
    use M_Box_System_Vars, only: vCpnBox
    use M_Box_Debug_Vars, only : LDebug_Init
    use M_Equil 
	use M_Basis, only : Basis_Change
    ! 
    logical,intent(out):: Ok
    integer :: iCp
    !
    call Info_("Box_Init_FluidBox")

    !// Build a Custom System for Bloc SYSTEM.BOX
    Ok = .true.
    call System_Build_Custom( & 
         & "SYSTEM.BOX", & 
         & vEle, vSpc, vMixFas, vCpnBox,  & 
         & TdgK, Pbar, & 
         & vCpn, Ok) 

    if(Ok) then    
       !// Compute Speciation
       call Equil_Calc("BOX")
	   !// Set Component = Elements and Inert  
	   call Basis_Change("DYN", vCpn)       
       !// Store the solution 
       call Store_TPCond(TdgK, Pbar)
       call Store_Composition_Species(vMolSpc_FluidBox)
       call Store_Composition_Component(vCpn, vCpnBox, vMolCpn_FluidBox)
    else
       !// Use System as System.Box
	   call Info_("Using SYSTEM Water for SYSTEM.BOX")
	   vCpn = vCpnBox       
       vMolSpc_FluidBox = vMolSpc_FluidSpc
       vMolCpn_FluidBox = vMolCpn_FluidSpc       
    end if

    !--- Info vCpn Debug
    if (Ldebug_Init) then
       call Info_("Cpn FluidBox")
       do iCp=1,size(vCpn) 
          call Component_Print(6,vEle,vSpc,vCpn(iCp))
       end do
    end if   

  end subroutine Box_Init_FluidBox

  !---

  subroutine Box_Init_FluidInject(Ok)
    !===============================================================
    ! Init Injected Fluid Composition
    !===============================================================
    ! Step 1. Read the system of constraint on injected fluid
    ! Step 2. Compute speciation
    ! Step 3. Store the composition 
    !===============================================================
    use M_Global_Vars, only: vEle,vSpc,vMixFas !//,vTPcond 
    use M_System_Vars, only: TdgK,Pbar,vCpn,vCpnTot
    use M_Dynam_Vars, only: vCpnInj
    use M_Basis_Vars, only: vOrdAq
    use M_System_Tools,only: System_Build_Custom 
    use M_Box_System_Vars, only: vCpnBox
    use M_Box_Debug_Vars, only : LDebug_Init
    use M_Equil 
    use M_Basis, only : Basis_Change
    implicit none
    ! 
    logical,intent(out):: Ok 
    integer :: iCp

    !---
    call Info_("Box_Init_FluidInject")    

    !// Build a Custom System for Bloc SYSTEM.INJECT 
    Ok = .true.
    
    call System_Build_Custom( & 
         & "SYSTEM.INJECT", & 
         & vEle, vSpc, vMixFas, vCpnBox, & 
         & TdgK, Pbar, & 
         & vCpn, Ok) !,SysType) 

    if(Ok) then    
       !// Compute Speciation
       call Equil_Calc("INJ")
	   !// Set Component = Elements and Inert  
	   call Basis_Change("DYN",vCpn)     
       !// Store the solution
       call Store_Composition_Species(vMolSpc_FluidInject)
       call Store_Composition_Component(vCpn, vCpnBox, vMolCpn_FluidInject)
    else
       !// Use System as System.Inject
       call Info_("Using SYSTEM Water for SYSTEM.INJECT")
	   vCpn = vCpnBox       
	   vMolSpc_FluidInject = vMolSpc_FluidSpc
       vMolCpn_FluidInject = vMolCpn_FluidSpc
    end if	
    
    !--- Info vCpn Debug
    if (Ldebug_Init) then
       call Info_("Cpn FluidInject")
       do iCp=1,size(vCpn) 
          call Component_Print(6,vEle,vSpc,vCpn(iCp))
       end do
    end if   

  end subroutine Box_Init_FluidInject

  !---

  subroutine Box_Init_VMolF_FluidInject
    !===============================================================
    ! Normalize Fluid Composition in Fluid Porosity
    !===============================================================
    use M_Box_Thermo_Vars, only : RhoF
    use M_Box_Vars,        only : PhiF, VBox, vMolFInj
    use M_Global_Vars,     only : vSpc
    use M_Basis_Vars,      only : vOrdAq
    use M_Box_Debug_Vars,  only: LDebug_Init
    implicit none
    !--- local variables
    real(dp) :: Factor
    real(dp) :: MassF
    integer :: iAq, iSpc
    !----
    call Info_("Box_Init_VMolF_FluidInject")
    do iAq= 1, nAq
       iSpc= vOrdAq(iAq)
       vMolFInj(iAq)= vMolSpc_FluidInject(iSpc)
    end do
    
    MassF = Zero
    
    !// Compute Mass Factor 
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       MassF = MassF + vMolFInj(iAq) * vSpc(iSpc)%Weitkg
    end do

    Factor = RhoF/MassF !  [(kg.m-3)] / kg ( m-3 )
    vMolFInj = vMolFInj* Factor

    !// Debug Check 
    if (LDebug_Init) then
       MassF = Zero       
       do iAq = 1, nAq
          iSpc = vOrdAq(iAq)
          MassF = MassF + vMolFInj(iAq) * vSpc(iSpc)%WeitKg 
       end do
       write(*,*) "Fluid Inject Density ( kg.m-3 ) =", MassF, RhoF
    end if
    
  end subroutine Box_Init_VMolF_FluidInject

    !---

  subroutine Box_Init_VMolF_FluidBox
    !===============================================================
    ! Normalize Fluid Composition in Fluid Porosity
    !===============================================================
    use M_Box_Thermo_Vars, only : RhoF
    use M_Box_Vars,        only : PhiF, VBox, vMolF
    use M_Global_Vars,     only : vSpc
    use M_Basis_Vars,      only : vOrdAq
    use M_Box_Debug_Vars, only: LDebug_Init
    implicit none
    !--- local variables
    real(dp) :: Factor
    real(dp) :: MassF
    integer :: iAq, iSpc
    !----
    call Info_("Box_Init_VMolF_Fluid")
    do iAq= 1, nAq
       iSpc= vOrdAq(iAq)
       vMolF(iAq)= vMolSpc_FluidBox(iSpc)
    end do
    
    MassF = Zero
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       MassF = MassF + vMolF(iAq) * vSpc(iSpc)%WeitKg 
    end do

    Factor = RhoF*PhiF*VBox / MassF !  [(kg.m-3)*m3] / kg ( adimensional factor )
    vMolF = vMolF* Factor

    MassF = Zero

     !// Debug Check 
    if (LDebug_Init) then
       do iAq = 1, nAq
          iSpc = vOrdAq(iAq)
          MassF = MassF + vMolF(iAq) * vSpc(iSpc)%WeitKg 
       end do

       write(*,*) "Mass Fluid Box ( kg ) =", MassF, RhoF*PhiF*VBox
    end if
     
  end subroutine Box_Init_VMolF_FluidBox

  !---

  subroutine Box_Init_VMolK_Rock
    !===============================================================
    ! Normalize Mineral volume Fractions in Solid Porosity
    ! Transform volumetric composition in mol composition
    !===============================================================
    use M_Box_Thermo_Vars, only : vRhoK
    use M_Box_Vars,        only : PhiF, vMolK, VBox, vPhiK
    use M_Global_Vars,     only : vKinFas, vKinModel, vFas
    use M_T_KinFas
    use M_Box_Debug_Vars, only: LDebug_Init
    implicit none

    !--- local variables
    real(dp) :: PhiS, Factor
    real(dp) :: FracMk, PhiMk, RhoMk, xMk, RhoMkMass
    integer  :: iMk, iFas
    !----
    call Info_("Box_Init_VMolK_Rock")

    vMolK = Zero
    PhiS = Zero
    do iMk = 1, nMk
       PhiMk = vKinFas(iMk)%Dat%PhiM
       PhiS = PhiS + PhiMk
    end do

    Factor = One / PhiS;

    !// update properties
    PhiS = (1-PhiF) 

    if (LDebug_Init) then
      write(*,*) &
      & "Rock-Fluid Volumes : PhiS = ", PhiS, " (1-PhiF) = ", One- PhiF    
    end if
    
    do iMk = 1, nMk

       !// Mineral Volume Fraction
       FracMk = vKinFas(iMk)%Dat%PhiM * Factor

       !// Phase density 
       iFas = vKinFas(iMk)%iFas
       RhoMk = One / vFas(iFas)%VolM3  ! mol.m-3
       RhoMkMass = vFas(iFas)%WeitKg / vFas(iFas)%VolM3 ! kg.m-3

       !// Species molar fraction 
       xMk = One ! pure mineral phase !

       !// vMolMk 
       vMolK(iMk) = VBox * ( FracMk * PhiS ) * RhoMk * xMk
       vRhoK(iMk)  = RhoMkMass
       vPhiK(iMk) = FracMk * PhiS
    end do

  end subroutine Box_Init_VMolK_Rock

  !---
  
  subroutine Box_Init_Texture_Rock()
    !========================================================================================
    !  Init texture properties for rock by using initial box properties
    !========================================================================================
    use M_Numeric_Const,      only: Pi
    use M_T_KinFas,      only: KinFas_InitTexture
    use M_Global_Vars,   only: vFas,vKinFas
    use M_Box_Kinetics_Vars
    use M_Box_Thermo_Vars,  only: vMolarVol
    use M_Box_Vars
    use M_Box_Debug_Vars, only: LDebug_Init
    !
    implicit none
    real(dp):: VolMinim
    integer :: iMk
    real(dp) :: Surf, SferNumber
    real(dp) :: VolMk, PhiMk

    !---

    SurfMinim= DensitMinim *Pi*4.D0       *RadiusMinim**2  !* VBox
    VolMinim=  DensitMinim *Pi*4.D0/3.D0  *RadiusMinim**3  !* VBox

    !// Debug
    if (LDebug_Init) then 
       write (*,'(A,G15.6)') "SurfMinim=", SurfMinim
       write (*,'(A,G15.6)') "VolMinim =", VolMinim
    end if
    
    do iMk = 1,nMk
     
      vMolarVol(iMk) = vFas(vKinFas(iMk)%iFas)%VolM3  !-> molar volumes of phases
      vKinMinim(iMk) = VolMinim / vMolarVol(iMk)      !-> mimimal mole number

      VolMk = VBox*vPhiK(iMk)
      PhiMk = vPhiK(iMk)

      call KinFas_InitTexture( &
      & vFas,         &   ! in
      & vKinFas(iMk), &   ! in
      & VolMk,        &   ! in
      & Surf,         &   ! out
      & SferNumber     ) ! out 

      vPhiM0(iMk)                 = PhiMk
      vSurfM0(iMk)                = MAX ( Surf, SurfMinim )

      vKinFas(iMk)%Dat%PhiM       = PhiMk
      vKinFas(iMk)%Dat%Surf       = vSurfM0(iMk) 
      vKinFas(iMk)%Dat%SferNumber = SferNumber

      !// Debug
      if (LDebug_Init) then 
      write(*,'(A15,2A,2(A,G15.6))') &
      & vKinFas(iMk)%NamKF, &
      & " Sat   =",vKinFas(iMk)%Dat%cSat,&
      & " PhiM  =",vKinFas(iMk)%Dat%PhiM, &
      & " Minim =",vKinMinim(iMk)
      end if
       
    enddo

    !// Debug
    if (LDebug_Init) then
      write(*,'(A,G15.6)') "VBox     ", VBox
      write(*,'(/,A,/)')   "vSurfM0 = Initial Surfaces of Minerals"

      do iMk=1,nMk       
        write(*,'(A15,1X,4(A,1X,G12.5,1X))') &
        & vKinFas(iMk)%NamKF,&
        & "vSurfM0    =", vSurfM0(iMk),&
        & "/Rho       =", vFas(vKinFas(iMk)%iFas)%WeitKg / vMolarVol(iMk),& !mineral density
        & "/vMolK     =", vMolK(iMk),&
        & "/SferNumber=", vKinFas(iMk)%Dat%SferNumber
      end do
    end if

  end subroutine Box_Init_Texture_Rock
 
 !---

  subroutine Box_Init_Texture_Rock_Old()
    !========================================================================================
    !  Init texture properties for rock by using initial box properties
    !========================================================================================
    use M_Numeric_Const,      only: Pi
    use M_T_KinFas,      only: KinFas_InitTexture
    use M_Global_Vars,   only: vFas,vKinFas
    use M_Box_Kinetics_Vars
    use M_Box_Vars
    !
    implicit none
    real(dp):: VolMinim
    integer :: iMk
    real(dp) :: Surf, SferNumber
    real(dp) :: VolMk, PhiMk
    !---
    SurfMinim= DensitMinim *Pi*4.D0       *RadiusMinim**2  !* VBox
    VolMinim=  DensitMinim *Pi*4.D0/3.D0  *RadiusMinim**3  !* VBox

    do iMk = 1,nMk
    
       VolMk = VBox*vPhiK(iMk)
       PhiMk = vPhiK(iMk)
      call KinFas_InitTexture( &
      & vFas,         &   ! in
      & vKinFas(iMk), &   ! in
      & VolMk,        &   ! in
      & Surf,         &   ! out
      & SferNumber     ) ! out 

      vSurfM0(iMk)                = Surf
      vPhiM0(iMk)                 = PhiMk

      vKinFas(iMk)%Dat%PhiM       = PhiMk
      vKinFas(iMk)%Dat%Surf       = Surf
      vKinFas(iMk)%Dat%SferNumber = SferNumber

    enddo

  end subroutine Box_Init_Texture_Rock_Old
 
!---

subroutine Box_Init_Info
    !===============================================
    ! Print the main results of Box_Init
    !===============================================
    use M_T_Component
    use M_Global_Vars,    only: vEle,vSpc
    use M_Box_Debug_Vars, only: LDebug_Init
    use M_Box_System_Vars
    use M_Box_Vars
    implicit none
    !--- 
    integer :: iCp, iAq, iSpc, iMk
    !---
    if (LDebug_Init) then
      call Info_("Box Fluid Compositions")

      write(*,*)" "
      write(*,'(A6,1X,A20,3A16)') &
      & 'iCp', 'vCpnBox(iCp)%NamCp', 'FluidSpc',   'FluidBox',  'FluidInject'
      write(*,'(A6,1X,A20,3A16)') &
      & '-----', '-----------------', '------------', '------------', '-------------'
      do iCp=1,nCp 
        write(*,'(I6,1X,A20,3G16.6)') &
        & iCp, vCpnBox(iCp)%NamCp, &
        & vMolCpn_FluidSpc(iCp), vMolCpn_FluidBox(iCp), vMolCpn_FluidInject(iCp)
      end do

      write(*,*)" "
      write(*,'(A6,A6, 1X,A20, 3A16)') &
      & 'iAq', 'iSpc', 'vSpc(iSpc)%NamSp',  'FluidSpc', 'FluidBox',  'FluidInject'
      write(*,'(A6,A6, 1X,A20, 3A16)') &
      & '-----', '-----', '------------------', '------------', '------------','-------------'
      do iAq=1,nAq 
        iSpc= vOrdAq(iAq)
        write(*,'(I6,I6, 1X,A20,3G16.6)') &
        & iAq, iSpc, vSpc(iSpc)%NamSp, &
        & vMolSpc_FluidSpc(iSpc),vMolSpc_FluidBox(iSpc), vMolSpc_FluidInject(iSpc)
      end do

      write(*,*)" "
      write(*,'(A6,A6, 1X,A20, A16, A16)') &
      & 'iAq', 'iSpc', 'vSpc(iSpc)%NamSp', 'vMolF',  'vMolFInj'
      write(*,'(A6,A6, 1X,A20, A16, A16)') &
      & '-----', '-----', '------------------', '------------', '-------------'
      do iAq=1,nAq 
        iSpc= vOrdAq(iAq)
        write(*,'(I6,I6, 1X,A20,2G16.6)') &
        & iAq, iSpc, vSpc(iSpc)%NamSp, vMolF(iAq), vMolFInj(iAq)
      end do

      write(*,*)" "
      write(*,'(A6,A6, 1X,A20, A16)') &
      & 'iMk', 'iSpc', 'vSpc(iSpc)%NamSp', 'vMolK'
      write(*,'(A6,A6, 1X,A20, A16)') &
      & '-----', '-----', '------------------', '------------'
      do iMk=1,nMk 
        iSpc= vOrdMk(iMk)
        write(*,'(I6,I6, 1X,A20,2G16.6)') &
        & iMk, iSpc, vSpc(iSpc)%NamSp, vMolK(iMk)
      end do
    end if

  end subroutine Box_Init_Info

  !--------- Small Utils ------------------------------------------------------

  subroutine Store_TPCond(T, P)
    use M_Box_Vars
    implicit none
    real(dp), intent(in) :: T, P
    !--
    TdgK= T
    PBar= P
    TdgKBox= TdgK
    PBarBox= Pbar
    
  end subroutine Store_TPCond

  !--

  subroutine Store_Composition_Species(vMolSpc_Fluid)
    !===============================================================
    ! Store fluid composition and index of species
    !===============================================================
    use M_T_Component
    use M_Equil 
    !--
    implicit none
    real(dp), intent(out) :: vMolSpc_Fluid(:)
    !--
    call Box_AqueousPhase_Get_VMolSpc(vMolSpc_Fluid)
    
  end subroutine Store_Composition_Species

  !---
  
  subroutine Store_Composition_Component(vCpnVar, vCpnMaster, vMolCpnVar)
    !===============================================================
    ! Store vCpnVar Composition by using vCpnMaster index
    ! vCpnVar must be a subset of vCpnMaster
    !===============================================================
    use M_T_Component
    implicit none
    type(T_Component), intent(in) :: vCpnVar(:)
    type(T_Component), intent(in) :: vCpnMaster(:)
    real(dp), intent(out) :: vMolCpnVar(:)

    !--- local variables
    integer :: iCp, jCp
    integer :: iEle

    !---     
    vMolCpnVar = Zero
    do iCp = 1, nCp
      iEle= vCpnVar(iCp)%iEle
      jCp= Component_EleIndex(iEle,vCpnMaster)
      if (jCp>0) then
      vMolCpnVar(jCp)= vCpnVar(iCp)%Mole
      else
      call Fatal_("Component not found in vCpnMaster")
      end if
    end do

  end subroutine Store_Composition_Component

  !---

  function Component_EleIndex(iEle,V) result(Index)
    !===============================================================
    ! Find the index of the component associated with Element iEle
    !===============================================================
    use M_T_Component
    implicit none
    integer,intent(in):: iEle
    type(T_Component),intent(in):: V(:)
    integer:: Index

    !-- local variable
    integer::I

    Index=0
    do I=1,size(V)
       if( iEle == V(I)%iEle ) then
          Index = I
          exit
       end if
    end do

  end function Component_EleIndex

  !---

  subroutine Box_AqueousPhase_Get_vMolSpc(vMolSpc)
    !------------------------------------------------------
    ! Store the amount of aqueous species only 
    !------------------------------------------------------
    use M_Global_Vars, only : vSpc
    implicit none
    real(dp), intent(out)::vMolSpc(:)
    integer :: iSpc, nSpc
    !---
    nSpc = size(vSpc)
    if (.not. (size(vMolSpc).eq.nSpc) ) then
      call Fatal_("Wrong number of Species")
    end if

    do iSpc=1, nSpc
      if (vSpc(iSpc)%Typ .eq. "AQU") then        
      vMolSpc(iSpc) = vSpc(iSpc)%Dat%Mole
      else
      vMolSpc(iSpc) = Zero
      end if
    end do

  end subroutine Box_AqueousPhase_Get_vMolSpc

end module M_Box_Init
