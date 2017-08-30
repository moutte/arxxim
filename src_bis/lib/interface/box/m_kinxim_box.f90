!**********************************************************************
!* File name:     M_KINXIM.f90
!*                procedures principales KINXIM 
!*
!* Author:        Anthony Michel (anthony.michel@ifp.fr)
!*
!* Created:       2006
!*
!* Modification:  J.Moutte 2008
!*
!*
!**********************************************************************

module M_KINXIM_BOX

  use M_Trace
  
  implicit none
  private
  
  public :: KINXIM_Initstatique
  public :: KINXIM_Speinit
  public :: KINXIM_InitDynamique
  public :: KINXIM_CalculDyn 
  public :: KINXIM_End

contains

  !---

  subroutine KINXIM_Initstatique
    !--------------------------------------------------
    ! Purpose   : Initialisation statique
    !             - lecture database
    !             - lecture directives et parametres user
    !             - allocation des tableaux 
    !--------------------------------------------------
    use M_Global_Build
    use M_Global_Tools
    use M_Files
    !
    implicit none
    logical:: OkCmd, Ok
    call Info_("KINXIM_Initstatique")
    !// Start Program
    write(*,'(A)') "ARXIM START"
    
    !// Open Trace Files
    call Trace_Init
    !
    !// Read The Arxim File 
    call Files_PathRead(Ok)
    if (.not. (Ok)) call Fatal_("Error in Reading Input File")
    call Files_BuildInput
    
    !// Read The Global Variables
    call Global_Zero
    call Global_Build 
    
  end subroutine KINXIM_Initstatique

  !---

  subroutine KINXIM_Speinit
    !--------------------------------------------------
    ! Purpose   : Calcul de speciation initiale
    !--------------------------------------------------
    
    use M_Box_Calc
    call Info_("KINXIM_Speinit")
    call Box_Calc_Init
    
  end subroutine KINXIM_Speinit

  !---

  subroutine KINXIM_InitDynamique
    !--------------------------------------------------
    ! Purpose   : Initialisation dynamique
    !             - preparation du calcul dynamique
    !--------------------------------------------------
    use M_Box_Calc
    call Info_("KINXIM_InitDynamique")
    call Box_Calc_Start
  end subroutine KINXIM_InitDynamique

  !---

  subroutine  KINXIM_CalculDyn
    !--------------------------------------------------
    ! Purpose   : Calcul Dynamique
    !             - calcul dynamique entre T et T+DT              
    !--------------------------------------------------
    use M_Box_Calc
    call Info_("KINXIM_CalculDyn")
    call Box_Calc_Restart
    call Box_Calc_Prepare
    call Box_Calc_Compute
    call Box_Calc_Finalize
    
end subroutine KINXIM_CalculDyn

  !---

  subroutine KINXIM_End
    !--------------------------------------------------
    ! Purpose   : 
    !             - Fermeture unites fichiers
    !             - Desallocation memoire 
    !--------------------------------------------------
    use M_Box_Calc
    call Info_("KINXIM_End")
    call Box_Calc_End
  end subroutine KINXIM_End

  !---

  subroutine InitStockVar
    !--------------------------------------------------
    ! Purpose   : Init the stock-variables
    !--------------------------------------------------
    call Info_("InitStockVar")
  end subroutine InitStockVar

end module M_KINXIM_BOX

