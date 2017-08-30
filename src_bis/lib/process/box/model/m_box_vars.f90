module M_Box_Vars

  !=====================================================
  ! BOX PROBLEM parameterS AND VARIABLES
  !=====================================================

  use M_Kinds
  use M_Trace

  use M_Box_Param_Vars

  implicit none
  public 

  !================= parameterS 
  public :: nAq
  public :: nMk
  public :: nCp
  public :: nPr
  public :: nAs    

  !================= VARIABLES 

  !// pvt conditions
  real(dp):: VBox     ! volume of the "box", m3
  real(dp):: TdgKBox  ! conditions Temperature User
  real(dp):: PbarBox  ! conditions Pressure User
  real(dp):: TdgK     ! conditions Temperature Reference
  real(dp):: Pbar     ! conditions Pressure Reference

  !// time
  real(dp):: TimeFactor= 1.0_dp  ! current Time Factor
  real(dp):: Time      ! current Time
  real(dp):: dTime     ! current dTime

  !// amounts
  real(dp):: EpsMineral    ! seuil annulation des mineraux vMolK
  real(dp), allocatable :: vMolF(:)    !1:nAq, nr moles of aqu. SPECIES in the box
  real(dp), allocatable :: vMolK(:)    !1:nMk, nr moles of kin'species in the box

  !// volumes
  real(dp)              :: PhiF        ! fluid volume fraction
  real(dp), allocatable :: vPhiK(:)    ! mineral volume fractions
  real(dp)              :: VolF        ! volume of Fluid, i.e. porosity*Vbox (and the gas ?)

  !// source injection
  real(dp):: Finj    ! flow rate input,  m3.s-1
  real(dp),allocatable ::  vMolFInj(:) ! input concentration, mol.s-1/m3 of fluid
  real(dp),allocatable ::  vCInj(:)    ! input concentration, mol.s-1/m3 of fluid
  real(dp),allocatable ::  vQInj(:)    ! mass rate input, mol.s-1

  !// output flow
  real(dp):: Fout    ! darcy flow rate output, m3.s-1
  real(dp):: Vout    ! real  flow rate output, m3.s-1 

  !================= SPECIAL

  real(dp), allocatable :: vKinMinim(:) !1:nMk, minimal nr moles of kin'species in the box
  integer,  allocatable :: vKinPrm(:)
  real(dp):: PhiInert                   ! volume fraction occupied by inactive minerals 

contains

  subroutine Box_Vars_New

    implicit none
    call Info_("Box_Newton_Vars_New")

    !================= VARIABLES 

    !// amounts
    allocate( vMolF(nAq) )
    allocate( vMolK(nMk) )

    !// volumes
    allocate( vPhiK(nMk) )
    
    !// source injection
    allocate( vMolFInj(nAq) )
    allocate( vCInj(nCp) )
    allocate( vQInj(nCp) )

    !================= SPECIAL

    allocate( vKinMinim(nMk) ) 
    allocate( vKinPrm(nMk) )
     
    Call Box_Vars_Zero

  end subroutine Box_Vars_New

  !---

  subroutine Box_Vars_Zero
    implicit none
    call Info_("Box_Newton_Vars_Zero")
    
    !================= VARIABLES 
    
    !// pvt conditions
    VBox= Zero
    TdgK= Zero
    Pbar= Zero

    !// time
    Time= Zero
    dTime= Zero

    !// amounts
    vMolF= Zero
    vMolK= Zero

    !// volumes
    PhiF= Zero
    vPhiK=Zero
    VolF= Zero

    !// source injection
    Finj= Zero
    vMolFInj= Zero 
    vCInj= Zero 
    vQInj= Zero 

    !// output flow
    Fout=Zero
    Vout=Zero

    !================= SPECIAL

    vKinMinim= Zero
    vKinPrm= 0
    PhiInert= Zero
  
  end subroutine Box_Vars_Zero

  !---

  subroutine Box_Vars_Delete

    implicit none
    call Info_("Box_Newton_Vars_Delete")
    !--
    call Box_Vars_Zero

    deallocate( vMolF)
    deallocate( vMolK)

    deallocate(vPhiK)

    deallocate(vMolFInj)
    deallocate(vCInj)
    deallocate(vQInj)
 
    deallocate(vKinMinim)
    deallocate(vKinPrm)

  end subroutine Box_Vars_Delete

  !---

end module M_Box_Vars
