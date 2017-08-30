module M_Box_Prepare

  !----------------------------------------------------------------
  ! Compute terms needed by the equations 
  ! Allow to control the level of explicitation of the model
  !----------------------------------------------------------------
  ! - TimeLoop  Explicit terms 
  ! - GammaLoop Explicit terms 
  ! - Implicit Terms for Residual ( values )
  ! - Implicit Terms for Jacobian ( values + derivatives )
  !----------------------------------------------------------------
  
  use M_Box_Residual_Vars
  
  public :: BOX_Prepare_TimeStep
  public :: BOX_Prepare_GammaStep
  public :: BOX_Prepare_Residual
  public :: BOX_Prepare_Jacobian

contains

  subroutine BOX_Prepare_Guess
    !-----------------------------------------------------------------
    ! Compute Solution Guess for Newton Loop
    !-----------------------------------------------------------------
    use M_Box_Vars
    use M_Box_Residual_Vars
    use M_Box_Debug_Vars, only : LDebug_Newton
    use M_Box_Load_Vars
    use M_Box_Utils
    use M_Box_Info
    implicit none
    
    call Info_("Box_Prepare_Guess")
    
    !// Residual Variables Initial Guess for Newton 
    vXf(1:nAq)   = vMolF(1:nAq)
    vLnXf(1:nAq) = Limited_Log(vXf(1:nAq))
    vXm(1:nMk)   = vMolK(1:nMk)
    
    if (LDebug_Newton) call Box_Info_Amount(vXf, vXm , 'Prepare_Guess', 6)
       
  end subroutine BOX_Prepare_Guess

  !---
  
  subroutine BOX_Prepare_TimeStep
    !-----------------------------------------------------------------
    ! Compute terms for the TimeLoop Sequence
    !-----------------------------------------------------------------
    use M_Box_RefState
    use M_Box_Volume
    use M_Box_Inject
    use M_Box_Outflow
    use M_Box_Kinetics
    
    implicit none
    
    call Info_("BOX_Prepare_TimeStep")

    ! Compute Fluid Volumes
    call BOX_Compute_Volume
    
    ! Compute Injection Rate
    call BOX_Compute_Inject
    call BOX_Compute_Outflow
    
    ! Compute Reference State Potential
    call Box_Compute_RefStatePotential
    call Box_Compute_LogK
    
    ! Explicit Kinetics Parameters
    call BOX_Compute_MineralStatus
    call BOX_Compute_KinActiv
    call BOX_Compute_KinSurf
    
  end subroutine BOX_Prepare_TimeStep

  !---
  
  subroutine BOX_Prepare_GammaStep(vX)
    !-----------------------------------------------------------------
    ! Compute terms for the GammaLoop Sequence
    !-----------------------------------------------------------------
    use M_Box_Entries
    use M_Box_Gamma
    use M_Box_KinThermo
    use M_Box_Kinetics

    implicit none
    real(dp),dimension(:),intent(in) :: vX
    
    call Info_("BOX_Prepare_GammaStep")

    ! Compute Natural Variables from Entries
    call BOX_Compute_Variables(vX) 

    ! Compute Thermodynamics 
    call BOX_Compute_Gamma 
    call BOX_Compute_WaterProperties
    
    ! Compute Thermodynamics 
    call BOX_Compute_Omega
    call Box_Compute_MineralStatus
    
  end subroutine BOX_Prepare_GammaStep

  !---
  
  subroutine BOX_Prepare_Residual(vX)
    !-----------------------------------------------------------------
    ! Compute implicit terms for Residual ( values )
    !-----------------------------------------------------------------
    use M_Box_Entries
    use M_Box_Gamma
    use M_Box_KinThermo
    use M_Box_Kinetics
    
    implicit none
    real(dp),dimension(:),intent(in) :: vX
    
    call Info_("BOX_Prepare_Residual")

    ! Compute Natural Variables From Entries
    call BOX_Compute_Variables(vX) 

    ! Compute Thermodynamics Water
    call Box_Compute_Activity(vXf, vLnXf, vLnai, dLnai_dLnXf)

    ! Compute Thermodynamics Water-Rock
    call Box_Compute_QsK( vLnXf, vLnai, dLnai_dLnXf, &      ! IN
         &      vLnQsK, dLnQsK_dLnXf, vQsK, dQsK_dLnXf )    ! OUT 
    
    ! Compute Vm and dVm
    call BOX_Compute_Kinetics(vLnXf, vXm, vVm, dVm_dLnXf, dVm_dXm, dVm_dLnXm)
    
  end subroutine BOX_Prepare_Residual

  !---
  
  subroutine BOX_Prepare_Jacobian(vX)
    !-----------------------------------------------------------------
    ! Compute implicit terms for Jacobian ( values + derivatives )
    !-----------------------------------------------------------------
    use M_Box_Entries
    use M_Box_Gamma
    use M_Box_KinThermo
    use M_Box_Kinetics
    use M_Box_Kinrate
    implicit none
    real(dp),dimension(:),intent(in) :: vX
    
    call Info_("BOX_Prepare_Jacobian")

    ! Compute Natural Variables From Entries
    call BOX_Compute_Variables(vX) 
    
    ! Compute Thermodynamics Water
    call Box_Compute_Activity(vXf, vLnXf, vLnai, dLnai_dLnXf)

    ! Compute Thermodynamics Water-Rock
    call Box_Compute_QsK( vLnXf, vLnai, dLnai_dLnXf, &      ! IN
         &      vLnQsK, dLnQsK_dLnXf, vQsK, dQsK_dLnXf )    ! OUT 

    ! Compute Vm and dVm
    call BOX_Compute_Kinetics(vLnXf, vXm, vVm, dVm_dLnXf, dVm_dXm, dVm_dLnXm)
    
  end subroutine BOX_Prepare_Jacobian
  
end module M_Box_Prepare
