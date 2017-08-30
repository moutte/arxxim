module M_Box_Kinetics

  use M_Kinds
  use M_Trace
  use M_Info_Value

  use M_Box_Kinetics_Vars

  implicit none
  private

  public :: BOX_Compute_Kinetics
  public :: BOX_Compute_MineralStatus
  public :: BOX_Compute_KinActiv
  public :: BOX_Compute_KinSurf
  
contains

  !---

  subroutine BOX_Compute_MineralStatus
    use M_Box_KinRate
    use M_Box_Vars
    use M_Global_Vars, only : vSpc, vFas, vKinFas, vKinModel
    use M_Box_Thermo_Vars, only : vOmegaK
    implicit none
    integer :: iMk
    
    do iMk = 1, nMk
      call Box_KinRate_SatState( &
      & vKinFas(iMk),   & !IN:  mineral: for cSatur, QsKseuil                        
      & vMolK(iMk),       & !IN:  Nr Moles of mineral, used to check nMol<=MolMinim
      & vKinMinim(iMk), & !IN:
      & vOmegaK(iMk),      & !IN:
      & QskIota,        & !IN:
      & vStatusK(iMk))   !OUT: saturation state        
    end do
    
  end subroutine BOX_Compute_MineralStatus

  !---

  subroutine BOX_Compute_KinActiv    
    implicit none
    integer :: iMk
    
    do iMk = 1, nMk
       vLKinActiv(iMk) = ( trim(vStatusK(iMk)) == "PRECIPI" .or. trim(vStatusK(iMk)) == "DISSOLU" ) 
    end do
    
  end subroutine BOX_Compute_KinActiv

  !---
  
  subroutine BOX_Compute_KinSurf   
    use M_Box_Vars
    use M_Global_Vars,       only : vFas, vKinFas
    use M_T_KinFas
    implicit none
    
    !-- local variables    
    integer :: iMk
    
    ! Surface   
    real(dp) :: SR
    real(dp) :: dSR_dXm
    real(dp) :: dSR_dLnXm
    
    do iMk = 1, nMk
      call KinFas_CalcSurf( &
      & vFas,               & !IN:
      & .true.,             & !IN:
      & vMolK(iMk),           & !IN:     Nr Moles of mineral                 
      & vKinFas(iMk),       & !IN/OUT: Mineral (update Radius or SferNumber)
      & SR,                 & !OUT:    SR = Reactive Surface                             
      & dSR_dLnXm,     & !OUT:    dSR | dLnXm
      & dSR_dXm )        !OUT:    dSR | dXm             
      vSrmK(iMk) = SR
    end do    
    
  end subroutine BOX_Compute_KinSurf
  
  !---

  subroutine BOX_Compute_Kinetics( vLnXf, vXm, vVm, dVm_dLnXf, dVm_dXm,  dVm_dLnXm)
    !==================================================
    ! Purpose   : compute Vm and dVm_dvLnX
    !==================================================
    use M_KinRate
    use M_T_KinModel
    use M_T_KinFas
    
    use M_Box_Vars
    use M_Global_Vars,       only : vSpc, vFas, vKinFas, vKinModel
    use M_Box_Solver_Vars,   only : Implicit_Surface, Implicit_ActivFactor
    use M_Box_Thermo_Vars,   only : vLnBuf, vLnGam, vLnAct
    use M_Box_RefState_Vars, only : vLnK_Mk
    use M_Box_System_Vars,   only : tNuMk, iAqH2O, vOrdMk
    use M_Box_KinRate
    use M_Numeric_Tools,only: vAbs

    implicit none
    real(dp), intent(in) :: vLnXf(:)
    real(dp), intent(in) :: vXm(:)
    real(dp), intent(out):: vVM(:)
    real(dp), intent(out):: dVM_dLnXf(:,:)
    real(dp), intent(out):: dVM_dLnXm(:,:)
    real(dp), intent(out):: dVM_dXm(:,:)

    !-- local variables
    logical ::  LSMOOTH_SURFACE     
    integer :: iMk  

    ! Surface
    real(dp) :: SR
    real(dp) :: XSR, dXSR_dSR
    real(dp) :: dSR_dXm
    real(dp) :: dSR_dLnXm
    
    ! Mineral index
    real(dp) :: QsK
    real(dp),dimension(1:nAq) :: dQsK_dLnXf

    ! Saturation State
    character:: cSatur
    logical :: LActiv

    ! Vm(Qsk) Term
    real(dp) :: VmQsK
    real(dp),dimension(1:nAq):: dVmQsK_dLnXf
    
    ! Vm(Act) Term
    real(dp) :: VmAct
    real(dp),dimension(1:nAq):: dVmAct_dLnXf
    
    !--------------------------------------------------
    
    call Info_("BOX_Compute_Kinetics")
    
    SR        = Zero
    dSR_dXm   = Zero
    dSR_dLnXm = Zero
    
    QsK        = Zero
    dQsK_dLnXf = Zero
    
    VmAct        = Zero
    dVmAct_dLnXf = Zero

    VmQsK        = Zero
    dVmQsK_dLnXf = Zero
    
    cSatur = ""
    
    call BOX_Info_Kinetics_Explicit

    call BOX_Info_Kinetics_Header  
    
    do iMk=1,nMk

    !call Info_Value("MINERAL ", vKinFas(iMk)%Name)
    !call Info_Value("INDEX   ", iMk)

    if ( vLKinActiv(iMk) ) then

      !//---------------- ACTIVE Mineral 

      !call Info_Value("ACTIVE MINERAL")

      !========= Mineral Surface ============
      SR=          Zero
      dSR_dLnXm=   Zero
      dSR_dXm=     Zero

      if (Implicit_Surface) then 
        call KinFas_CalcSurf( &
        & vFas,               & !IN:
        & .true.,             & !IN:
        & vXm(iMk),           & !IN:     Nr Moles of mineral                 
        & vKinFas(iMk),       & !IN/OUT: Mineral (update Radius or SferNumber)
        & SR,                 & !OUT:    SR = Reactive Surface                             
        & dSR_dLnXm,     & !OUT:    dSR | dLnXm
        & dSR_dXm )        !OUT:    dSR | dXm

        vKinFas(iMk)%Dat%Surf = SR !! bug test !!
      else
        SR= vSrmK(iMk)
      end if

      !========= compute QsK =================
      QsK=       Zero
      dQsK_dLnXf= Zero

      call KinRate_CalcQsK( &
      & nCp,          & !IN: nCi Parameter
      & 0,            & !IN: nCx Parameter
      & vLnK_Mk(iMk), & !IN: LnK Mineral
      & tNuMk(iMk,1:nPr), & !IN: Mineral Reaction Stokio
      & vLnXf,        & !IN:  Ln( Aqu species mole number ) 
      & vLnGam,       & !IN:  Ln( Aqu species Activity Coeff Gamma)                             
      & vLnBuf,       & !IN:  Ln( Activity of Buffered Species )
      & vLnAct(iAqH2O),  & !IN:  Ln( Activity of Solvent)                                  
      & QsK,          & !OUT  Qsk = mineral saturation index
      & dQsK_dLnXf)     !OUT  dQsk | dLnXf 

      !// cut extreme values of QsK
      QsK = MIN(QsK,QsK_Max)
      QsK = MAX(QsK,QsK_Min)
          
      vKinFas(iMk)%Dat%QsK = QsK !! bug test !!

      !=========== From QsK deduce cSatur, for branching to dissol/precip ====
      cSatur = 'none'

      call KinRate_SatState( &
      & vKinFas(iMk),   & !IN:  mineral: for cSatur, QsKseuil                        
      & vXm(iMk),       & !IN:  Nr Moles of mineral, used to check nMol<=MolMinim
      & Zero,           & !IN:  Check pure mineral presence
      & QsK,            & !IN:
      & QskIota,        & !IN:
      & cSatur)           !OUT: saturation state  

      vKinFas(iMk)%Dat%cSat = cSatur !! bug test !!      
      LActiv = ( trim(cSatur)=="PRECIPI" .or. trim(cSatur)=="DISSOLU" ) 

      !=========== Compute QsK Factor ===========
      VmQsK= Zero
      dVmQsK_dLnXf= Zero

      if ( LActiv ) then
      call KinRate_CalcQsKFactor( & 
      & cSatur,      &        !IN
      & vKinModel(vKinFas(iMk)%iKin),&  !IN: kinetic model of mineral
      & QsK,         &        !IN
      & dQsK_dLnXf,   &       !IN
      & VmQsK,       &        !OUT: the QsK function in the rate law            
      & dVmQsK_dLnXf(1:nAq) ) !OUT: dVmQsk | dLnXf
      end if
       
      !=========== Compute Activity Factor ===========
      VmAct= One
      dVmAct_dLnXf= Zero

      call KinModel_CalcCoeffs(vKinModel(vKinFas(iMk)%iKin),TimeFactor, TdgK)

      if ( LActiv ) then
        !write(*,*) "KinRate_CalcActivFactor"
        call Box_KinRate_CalcActivFactor(& 
        & cSatur, &                        !IN:  Saturation State
        & vKinModel(vKinFas(iMk)%iKin), &  !IN:  kinetic model of mineral
        & vLnAct, &                        !IN:  Aqueous species Ln(Activity)
        & VmAct,  &                        !OUT: the Activity function in the rate law
        & dVmAct_dLnXf(1:nAq)  )           !OUT: dVmAct | dLnXf
      end if

      LSMOOTH_SURFACE = .false.
      if ( LSMOOTH_SURFACE ) then

        !========== Surface Modifiers ============================
        call Smooth_Function ( SR, XSR, dXSR_dSR, 0.5*SurfMinim, SurfMinim )

        !========== Compose Vm Function and Derivatives ==========
        vVm(iMk)= XSR * VmAct * VmQsK

        dVM_dLnXf(iMk,1:nAq)= &
        &    XSR * VmAct * dVmQsK_dLnXf(1:nAq) &
        &  + XSR * dVmAct_dLnXf(1:nAq) * VmQsK 

        dVM_dXm(iMk,iMk)= dXSR_dSR * dSR_dXm* VmAct * VmQsK
        dVm_dLnXm(iMk,iMk) = Zero

      else

        !========== Compose Vm Function and Derivatives ==========
        vVm(iMk)= SR * VmAct * VmQsK

        dVM_dLnXf(iMk,1:nAq)= &
        &    SR * VmAct * dVmQsK_dLnXf(1:nAq) &
        &  + SR * dVmAct_dLnXf(1:nAq) * VmQsK 

        dVM_dXm(iMk,iMk)= dSR_dXm* VmAct * VmQsK
        dVm_dLnXm(iMk,iMk) = Zero
      
      end if
       
    else
       
      !//---------------- NOT-ACTIVE Mineral 

      !call Info_Value("NOT-ACTIVE MINERAL")

      vVm(iMk)             = Zero
      dVM_dLnXf(iMk,1:nAq) = Zero
      dVM_dXm(iMk,iMk)     = Zero
      dVm_dLnXm(iMk, iMk)  = Zero

    end if

      call BOX_Info_Kinetics(iMk, vSpc(vOrdMk(iMk))%NamSp, cSatur, LActiv, &
      & vXm(iMk), QsK, SR, VmAct, VmQsK, vVm(iMk))

    end do

  end subroutine BOX_Compute_Kinetics
  
  !---

  subroutine Smooth_Function ( x, f, df_dx, a, b )
    real(dp), intent(in)  :: x, a, b
    real(dp), intent(out) :: f, df_dx
    real(dp) :: delta 
    
    !//-----------------------
    delta = b-a

    if ( x < a ) then
      f     = 0.D0
      df_dx = 0.D0
    elseif ( x < b ) then
      f     = 0.5*( (x-a) / delta )**2
      df_dx = delta * ( (x-a) / delta )
    else
      f     = x 
      df_dx = 1.D0
    end if

  end subroutine Smooth_Function
  
  !---
  
  subroutine Box_Info_Kinetics_Header()
    !==================================================
    ! Purpose   : Info Debug Kinetics Header
    !==================================================
    use M_Box_Debug_Vars, only : LDebug_KinRate
    implicit none
    
    if (LDebug_KinRate) then 
      write(*,'(A6,3(1X,1A12), 6A14)') &
      & '  iMk ', 'iMkSpc%Name ', 'cSatur ', 'LActiv  ', &
      & '  Xm  ', 'QsK ',&
      & '  SRm ', 'VmAct' ,  'VmQsK ' ,&
      & '  Vm '           

      write(*,'(A6,3(1X,1A12), 6A14)') &
      & '------',  '------------', '------------',&
      & '-----------', '-----------', '-----------', &
      & '-----------', '-----------', '-----------', &
      & '-----------'
    end if
     
  end subroutine Box_Info_Kinetics_Header
  
  !---

  subroutine BOX_Info_Kinetics( iMk, NameMk, cSatur, LActiv, XmMk, QsK, SR, VmAct, VmQsK , VmMk)
    !=====================================================
    ! Purpose   : Info Debug Kinetics Line for Mineral Mk
    !=====================================================
    use M_Box_Debug_Vars, only : LDebug_KinRate
    use M_Box_Vars
    
    implicit none
    
    integer :: iMk
    character(len=*) :: NameMk
    logical :: LActiv
    character(len=12) :: FlagActiv
    real(dp), intent(in) :: XmMk,  QsK, SR,  VmAct, VmQsK, VmMk
    character(len=*) :: cSatur
    !
    if (LDebug_KinRate) then 
      FlagActiv = "NOT-ACTIVE"
      if ( LActiv ) FlagActiv = "ACTIVE"
      !
      write(*,'(I6,3(1X,1A12),6G14.5)') &
      & iMk, NameMk, cSatur, FlagActiv, &
      & XmMk, QsK,&
      & SR, VmAct, VmQsK ,&
      & VmMk 

    end if
   
  end subroutine BOX_Info_Kinetics
  
  !---
  
  subroutine BOX_Info_Kinetics_Explicit
    use M_Global_Vars,       only : vSpc
    use M_Box_Debug_Vars,    only : LDebug_KinRate
    use M_Box_System_Vars
    use M_Box_Vars
    implicit none
    integer :: iMk  
    character(len=12) :: FlagKinActiv

    call Info_("BOX_Info_Kinetics_Explicit ")

    if (LDebug_KinRate) then 
      write(*,'(A6,3(1X,1A12), 3A14)') &
      & '  iMk ', 'iMkSpc%Name ', 'vStatusK', 'vKinActiv', &
      & '  MolK ', 'SrmK', 'KinRateK'

      write(*,'(A6,2(1X,1A12), 3A14)') &
      & '------', '------------', '------------',& 
      & '-----------', '-----------', '-----------'

      do iMk = 1, nMk
      
        FlagKinActiv = "NOT-ACTIVE"
        if ( vLKinActiv(iMk) ) FlagKinActiv = "ACTIVE"
        !
        write(*,'(I6,3(1X,1A12),3G14.5)') &
        & iMk, vSpc(vOrdMk(iMk))%NamSp, vStatusK(iMk), FlagKinActiv, &
        & vMolK(iMk), vSrmK(iMk), vKinRateK(iMk)

      end do
    end if
  
  end subroutine BOX_Info_Kinetics_Explicit
       
end module M_Box_Kinetics
