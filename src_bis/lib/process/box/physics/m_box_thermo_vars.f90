module M_Box_Thermo_Vars

  !==================================================================================
  ! THERMODYNAMICAL VARIABLES
  !==================================================================================

  use M_Kinds
  use M_Trace
  
  use M_Box_Param_Vars

  implicit none
  public

  !// molar weigth
  real(dp), allocatable:: vWeitSp(:)  !1:nAq, vSpc(vOrdAq(1:nAq))%WeitKg
  real(dp), allocatable:: vWeitCp(:)  !1:nCp, vEle(vOrdCp(1:nCp))%WeitKg

  !// density and molar volume
  real(dp):: RhoF                      !fluid density, kg.m^-3
  real(dp), allocatable:: vRhoK(:)     !1:nMk mineral density, kg.m^-3
  real(dp), allocatable:: vMolarVol(:) !1:nMk, Molar volume (=vFas(vKinFas(1:nMk)%iFas)%V)

  !// chemical potential and activity
  real(dp), allocatable:: vLnMolal(:)  !Ln(Molalityy) aqu'species in box
  real(dp), allocatable:: vLnAct(:)    !Ln(Activity) aqu'species in box
  real(dp), allocatable:: vLnGam(:)    !Ln(Activ'Coeff) aqu'species in box 
  real(dp), allocatable:: vLnBuf(:)    !Ln(Activ'Coeff) of buffers au_species in the box 
  
  !// minerals
  real(dp), allocatable:: vOmegaK(:)  !QsK mineral saturation index
 

  real(dp):: pH
  real(dp):: pE
  real(dp):: Ionic
  
contains 

  subroutine Box_Thermo_Vars_New
    use M_Box_Vars
    implicit none
    
    call Info_("Box_Thermo_Vars_New")
    
    allocate(vWeitSp(nAq))
    allocate(vWeitCp(nCp))

    allocate(vRhoK(nMk))
    allocate(vMolarVol(nMk))
    
    allocate(vLnMolal(nAq))
    allocate(vLnAct(nAq))
    allocate(vLnGam(nAq))
    allocate(vLnBuf(nAq))

    allocate(vOmegaK(nMk)) 

    call Box_Thermo_Vars_Zero
    
  end subroutine Box_Thermo_Vars_New

!---

  subroutine Box_Thermo_Vars_Zero
    use M_Box_Vars
    implicit none
    
    call Info_("Box_Thermo_Vars_Zero")
   
    RhoF= Zero
    vRhoK(:)= Zero
    
    vLnMolal(:)= Zero
    vLnAct(:)= Zero
    vLnGam(:)= Zero
    vLnBuf(:)= Zero
    
    vOmegaK = Zero

    pH= Zero
    pE= Zero
    Ionic= Zero
    
  end subroutine Box_Thermo_Vars_Zero

  !---
  
  subroutine Box_Thermo_Vars_Delete
    implicit none

    call Info_("Box_Thermo_Vars_Del")
    
    deallocate(vWeitSp)
    deallocate(vWeitCp)

    deallocate(vRhoK)
    deallocate(vMolarVol)
    
    deallocate(vLnMolal)
    deallocate(vLnAct)
    deallocate(vLnGam)
    deallocate(vLnBuf)

    deallocate(vOmegaK)
            
  end subroutine Box_Thermo_Vars_Delete

end module M_Box_Thermo_Vars
