module M_Box_Kinetics_Vars

  !==================================================================================
  ! KINETIC REACTION AND MINERAL TEXTURE VARIABLES
  !==================================================================================

  use M_Kinds
  use M_Trace

  use M_T_KinModel,only: T_KinModel
  use M_Box_Param_Vars

  implicit none
  public

  !// reaction_limiter mode
  real(dp):: Reaction_limiter = 1.D0      !kinetics limiter 
  real(dp),parameter:: QsK_Max= 1.D+30    !over  saturation limiter 
  real(dp),parameter:: QsK_Min= 1.D-30    !under saturation limiter 

  !// texture
  character(len=7):: sModelSurf= "SPHERE" !"CRUNCH" !"SPHERE",
  real(dp), allocatable:: vSurfM0(:)      !1:nMk, surface of phase at ref'time
  real(dp), allocatable:: vPhiM0(:)       !1:nMk, vol'fraction of phase at ref'time
  
  ! default texture parameters for secondary minerals
  real(dp):: RadiusMinim= 1.0D-9  !radius (meter), 1.0D-9= 1 nm
  real(dp):: DensitMinim= 1.D+15  !number of grains / m^3 fluid volume
  real(dp):: SurfMinim            !=1.0D-6 ! lower limit for surface (per volume) 

  !// mineral reaction
  real(dp) :: QskIota

  !// mineral status
  integer:: nMkA0                       !nMkA !,= count(vLKinActiv), number of "active" kinetic phases
  logical,  allocatable:: vLKinActiv(:) !Active Minerals
  logical,  allocatable:: vKinPrm(:)    !Primary Minerals Logical

  !// kinetics rate Vm
  real(dp), allocatable:: vSrmK(:)     !1:nMk, Mineral Explicit Surface
  real(dp), allocatable:: vVmAct(:)    !1:nMk, kinetic constant * (activators/inhibitors)
  real(dp), allocatable:: vVmQsK(:)    !1:nMk, affinity factor (e.g. QsK-1)
  real(dp), allocatable:: vKinRateK(:) !1:nMk, value of vVm at begin time step (end former time step)
  character(len=10), allocatable:: vStatusK(:) !1:nMk, mineral status
  
contains

  subroutine Box_Kinetics_Vars_New

    implicit none

    call Info_("Box_Kinetics_Vars_New")

    allocate(vSurfM0(1:nMk))  
    allocate(vPhiM0(1:nMk))  

    allocate(vVmAct(1:nMk))  
    allocate(vSrmK(1:nMk))
    allocate(vVmQsK(1:nMk))   
    allocate(vKinRateK(1:nMk))   
    allocate(vStatusK(1:nMk))   

    allocate(vLKinActiv(1:nMk)) 
    allocate(vKinPrm   (1:nMk)) 
    
    call Box_Kinetics_Vars_Zero

  end subroutine Box_Kinetics_Vars_New

  !--

  subroutine Box_Kinetics_Vars_Zero

    implicit none

    call Info_("Box_Kinetics_Vars_Zero")

    vSurfM0(:)= Zero
    vPhiM0(:)= Zero

    vVmAct(:)= Zero
    vSrmK(:)= Zero
    vVmQsK(:)= Zero
    vKinRateK(:)= Zero
    vStatusK(:)= 'none'

    nMkA0= 0
    vLKinActiv(:)= .false.
    vKinPrm(:)= .false.

  end subroutine Box_Kinetics_Vars_Zero

  !--

  subroutine Box_Kinetics_Vars_Delete

    implicit none

    call Info_("Box_Kinetics_Vars_Delete")

    call Box_Kinetics_Vars_Zero

    deallocate(vSurfM0)  
    deallocate(vPhiM0)  

    deallocate(vSrmK)
    deallocate(vVmAct)   
    deallocate(vVmQsK)
    deallocate(vKinRateK)   
    deallocate(vStatusK)   

    deallocate(vLKinActiv) 
    deallocate(vKinPrm) 

  end subroutine Box_Kinetics_Vars_Delete

end module M_Box_Kinetics_Vars
