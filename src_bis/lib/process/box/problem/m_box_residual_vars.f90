module M_Box_Residual_Vars

  !==================================================================================
  ! BOX VARIABLES FOR RESIDUAL AND JACOBIAN
  !==================================================================================

  use M_Kinds
  use M_Trace
  use M_Box_Vars, only : nMk, nAq, nAs, nCp
  
  implicit none

  public
  
  !// Problem Size
  integer, private :: nX 

  !// Natural Variables
  real(dp),allocatable:: vTotf(:)
  real(dp),allocatable:: vXf(:)
  real(dp),allocatable:: vLnXf(:)
  real(dp),allocatable:: vXm(:)
  
  !// Equations
  real(dp),allocatable:: vEqBal(:)
  real(dp),allocatable:: vEqEquil(:)
  real(dp),allocatable:: vEqKin(:)  
  
  !// System Value, Equations and Jacobian
  real(dp),allocatable:: vXSystem(:)
  real(dp),allocatable:: vEqSystem(:)
  real(dp),allocatable:: tJacSystem(:,:)
  
  !// Thermodynamics
  real(dp), allocatable:: vLnai(:)         
  real(dp), allocatable:: dLnai_dLnXf(:,:)  ! Derivatives

  real(dp), allocatable:: vLnQsK(:)        
  real(dp), allocatable:: dLnQsK_dLnXf(:,:) ! Derivatives
  
  real(dp), allocatable:: vQsK(:)         
  real(dp), allocatable:: dQsK_dLnXf(:,:)   ! Derivatives

  !// Kinetics
  real(dp), allocatable:: vVM(:)
  real(dp), allocatable:: dVM_dLnXf(:,:)    ! Derivatives
  real(dp), allocatable:: dVM_dLnXm(:,:)    ! Derivatives
  real(dp), allocatable:: dVM_dXm(:,:)      ! Derivatives
  
contains

  !--

  subroutine Box_Residual_Vars_New
    implicit none
    
    call Info_("Box_Residual_Vars_New")
    
    !// Problem Size
    nX= nAq + nMk
    
    !// Physical Variables
    allocate(vTotf(nCp))
    allocate(vXf(nAq))
    allocate(vLnXf(nAq))
    allocate(vXM(nMk))

    !// Equations
    allocate(vEqBal(nCp))
    allocate(vEqEquil(nAs))
    allocate(vEqKin(nMk))
    
    !// System Value, Equations and Jacobian
    allocate(vXSystem(nX))
    allocate(vEqSystem(nX))
    allocate(tJacSystem(nX, nX))
    
    !// Thermodynamics
    allocate(vLnai(nAq))         
    allocate(dLnai_dLnXf(nAq,nAq))  

    allocate(vLnQsK(nMk))        
    allocate(dLnQsK_dLnXf(nMk,nAq)) 
  
    allocate(vQsK(nMk))         
    allocate(dQsK_dLnXf(nMk,nAq))   

    !// Kinetics
    allocate(vVm(nMk))
    allocate(dVm_dLnXf(nMk,nAq))    
    allocate(dVm_dLnXm(nMk,nMk))    
    allocate(dVm_dXm(nMk,nMk))      

    call Box_Residual_Vars_Zero
    
  end subroutine Box_Residual_Vars_New

  !---

  subroutine Box_Residual_Vars_Zero
    implicit none
    
    call Info_("Box_Resiual_Vars_Zero")

    !// Problem Size
    nX= nX 
    
    !// Physical Variables
    vTotf(:)= Zero
    vXf(:)= Zero
    vLnXf(:)= Zero
    vXm(:)= Zero
    
    !// Equations
    vEqBal(:)= Zero
    vEqEquil(:)= Zero
    vEqKin(:)= Zero  
    vEqSystem(:)= Zero

    !// Jacobian
    vXSystem(:)= Zero
    vEqSystem(:)=Zero
    tJacSystem(:,:)=Zero

     !// Thermodynamics
    vLnai(:)= Zero
    dLnai_dLnXf(:,:)= Zero  

    vLnQsK(:)= Zero        
    dLnQsK_dLnXf(:,:)= Zero 
  
    vQsK(:)= Zero         
    dQsK_dLnXf(:,:)= Zero   

    !// Kinetics
    vVm(:)= zero
    dVm_dLnXf(:,:)= Zero    
    dVm_dLnXm(:,:)= Zero    
    dVm_dXm(:,:)= Zero      
    
  end subroutine Box_Residual_Vars_Zero

  !---
  
  subroutine Box_Residual_Vars_Delete
    implicit none
    
    call Info_("Box_Residual_Vars_Delete")
    
    call Box_Residual_Vars_Zero
    
    !// Problem Size
    nX = 0

    !// Physical Variables
    deallocate(vTotf)
    deallocate(vXf)
    deallocate(vLnXf)
    deallocate(vXM)

    !// Equations
    deallocate(vEqBal)
    deallocate(vEqEquil)
    deallocate(vEqKin)
    
    !// System Value, Equations and Jacobian
    deallocate(vXSystem)
    deallocate(vEqSystem)
    deallocate(tJacSystem)
    
    !// Thermodynamics
    deallocate(vLnai)         
    deallocate(dLnai_dLnXf)

    deallocate(vLnQsK)        
    deallocate(dLnQsK_dLnXf)
  
    deallocate(vQsK)         
    deallocate(dQsK_dLnXf)

    !// Kinetics
    deallocate(vVm)
    deallocate(dVm_dLnXf)
    deallocate(dVm_dLnXm)
    deallocate(dVm_dXm)
    
  end subroutine Box_Residual_Vars_Delete


end module M_Box_Residual_Vars
