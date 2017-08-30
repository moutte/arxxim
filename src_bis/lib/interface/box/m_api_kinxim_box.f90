!! *******************************************************************
!! * File name:    M_API_KINXIM
!! *
!! * Purpose:      API module de geochimie KINXIM
!! *
!! * Author:       Anthony Michel (anthony.michel@ifp.fr)
!! *
!! * Created:      2006
!! *
!! ********************************************************************

module M_API_KINXIM_BOX

  use M_Trace
  use M_Kinds,only: dp
  implicit none
  private

  character (len=10) :: APITAG = 'KIMXIM-BOX'

  ! API STANDARD
  public :: Init_API_KINXIM_INTERACT
  public :: Init_API_KINXIM
  public :: Compute_API_KINXIM
  public :: Output_API_KINXIM
  public :: End_API_KINXIM

  ! CALCUL DYNAMIQUE NORMAL
  public :: Compute_API_KINXIM_noargs

  ! EXTRACTION VARIABLES POUR OUTPUT
  public :: Extract_KINXIM
  
contains

  !---

  subroutine Init_API_KINXIM_INTERACT
  !lt !DEC$ ATTRIBUTES DLLEXPORT :: Init_API_KINXIM_INTERACT
    !--------------------------------------------------
    ! Purpose : 
    ! - LECTURE INPUT FILE
    ! - INIT BOX
    !--------------------------------------------------    
    use M_KINXIM_BOX
    use M_EXCHANGE_KINXIM_BOX
    use M_BOX_COUPLER_VARS
    !--- 
    implicit none
    real(dp) :: Tfinal_save
    call Info_("Init_API_KINXIM")
    write(*,*) "INIT API [", APITAG, "]" 
    
    ! lecture database et directives statiques + allocation memoire
    !call KINXIM_InitStatique
    
    ! calcul de speciation initiale si necessaire
    call KINXIM_Speinit

    ! lecture database et directives dynamiques + allocation memoire
    call KINXIM_InitDynamique
    
    ! calcul T = 0
    !call KINXIM_Get_TMore(Tfinal_save)
    !call KINXIM_Set_TMore(0.D0)
    !call KINXIM_CalculDyn   
    !call KINXIM_Set_TMore(TFinal_save)

    LCouplerActive = .false.

  end subroutine Init_API_KINXIM_INTERACT

  subroutine Init_API_KINXIM
  !lt !DEC$ ATTRIBUTES DLLEXPORT :: Init_API_KINXIM
    !--------------------------------------------------
    ! Purpose : 
    ! - LECTURE INPUT FILE
    ! - INIT BOX
    !--------------------------------------------------    
	use M_KINXIM_BOX
    use M_BOX_COUPLER_VARS
    !--- 
    implicit none
    call Info_("Init_API_KINXIM")
    write(*,*) "INIT API [", APITAG, "]" 
    ! lecture database et directives statiques + allocation memoire
    call KINXIM_InitStatique
    
    ! calcul de speciation initiale si necessaire
    call KINXIM_Speinit

    ! lecture database et directives dynamiques + allocation memoire
    call KINXIM_InitDynamique
    
    ! mode couple
    LCouplerActive = .true.

  end subroutine Init_API_KINXIM

  !---
  
  subroutine End_API_KINXIM
    !--------------------------------------------------
    ! Purpose   : 
    ! - FIN DU CALCUL
    ! - SORTIES KINXIM
    !--------------------------------------------------
!lt !DEC$ ATTRIBUTES DLLEXPORT :: End_API_KINXIM
	 use M_KINXIM_BOX
    use M_OUTPUT_KINXIM_BOX
    implicit none
    !--- 
    call Info_("End_API_KINXIM")
    call KINXIM_end
    
    write(*,*) "END API [", APITAG, "]"
    write(*,*) "PERFECT"
    
  end subroutine End_API_KINXIM
  
  !---
  
  subroutine Output_API_KINXIM
!lt !DEC$ ATTRIBUTES DLLEXPORT :: Output_API_KINXIM
    !--------------------------------------------------
    ! Purpose   : FIN ET SORTIES KINXIM
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2006
    !--------------------------------------------------
    use M_OUTPUT_KINXIM_BOX
    implicit none
    !--- 
    call Info_("Output_API_KINXIM")
    call KINXIM_sorties_COORES
   
  end subroutine Output_API_KINXIM
  
  !---

  subroutine Compute_API_KINXIM_noargs
!lt !DEC$ ATTRIBUTES DLLEXPORT :: Compute_API_KINXIM_noargs
    !--------------------------------------------------
    ! Purpose   : CALCUL DYNAMIQUE NORMAL
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2006
    !--------------------------------------------------
	
    use M_KINXIM_BOX
    implicit none   
    !---
    call Info_("Compute_API_KINXIM_noargs")
    write(*,*) "COMPUTE API [", APITAG, "] NOARGS"
    
    ! calcul dynamique effectif
    
    call KINXIM_CalculDyn   

  end subroutine Compute_API_KINXIM_noargs

  !---

  subroutine Compute_API_KINXIM( &
       TINIT, DT, &
       NE, PHIM, &
       PRESSURE, TEMPER, &
       QCin, Vout, &
       VOLUME, &
       dNE, dPHIM, &
       NEout,&
       dtzero, dtend, niter, nbtimestep, &
       nbEAQ, nbC, nbMIN, &
       & LTAGCELLOUT, ierror )  
!lt !DEC$ ATTRIBUTES DLLEXPORT :: Compute_API_KINXIM
    !--------------------------------------------------
    ! Purpose   : CALCUL DYNAMIQUE AVEC RESET VAR
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2006
    !--------------------------------------------------
    use M_mathconst
    use M_KINXIM_BOX
    use M_EXCHANGE_KINXIM_BOX
    implicit none   

    ! IN 
    integer, intent(in) :: nbEAQ
    integer, intent(in) :: nbC
    integer, intent(in) :: nbMIN
    logical :: LTAGCELLOUT
    real (dp), intent(in)  :: TINIT
    real (dp), intent(in)  :: DT
    real (dp), intent(in)  :: VOLUME
    real (dp), intent(in)  :: Vout    ! Debit volumique sortant
    real (dp), intent(in), dimension(nbC)    :: QCin    ! Debit molaire entrant / element
    real (dp), intent(in), dimension(nbEAQ)  :: NE      ! Especes aqueuses     
    real (dp), intent(in), dimension(nbMIN)  :: PHIM    ! porosite Mineraux
    real (dp), intent(in) :: PRESSURE, TEMPER           ! pression et temperature
    real (dp), intent(in) :: dtzero

    ! OUT
    real (dp), dimension(nbEAQ), intent(out) :: dNE    ! variation des especes aqueuses
    real (dp), dimension(nbMIN), intent(out) :: dPHIM  ! variation des Mineraux
    real (dp), dimension(nbEAQ), intent(out) :: NEout  ! variation des Mineraux
    integer, intent(out) :: niter
    integer, intent(out) :: nbtimestep
    real (dp), intent(out) :: dtend
    integer, intent(out) :: ierror

    !--- local variables
    real (dp), dimension(nbEAQ)  :: NE_DT      ! Especes aqueuses
    real (dp), dimension(nbMIN)  :: PHIM_DT    ! porosite Mineraux
    real (dp) :: dtmini
    integer :: nbtimestep_init  !, i
    
    !--- instructions

    call Info_("Compute_API_KINXIM")

    !---- RESET + TRANSformatION NECESSAIRES (M3/LITRE, s en y etc ...)
    
    ! pas de temps
    call KINXIM_set_t ( TINIT   )
    call KINXIM_set_tmore ( DT  ) 
    dtmini = min ( dtzero, DT  ) ! pour ne pas depasser DT
    dtend = 0
    niter = 0
    ierror = 0	
    !  
    call KINXIM_set_dt ( dtmini )
        
!!$    call KINXIM_set_NE   ( ( NE  / VOLUME  )* m3_in_liter ) ! boite de 1.d-3 m3 
!!$    call KINXIM_set_PHIM ( PHIM )
!!$    call KINXIM_set_QCin ( ( Qcin  / VOLUME ) * m3_in_liter / year_in_second  )
!!$    call KINXIM_set_Vout ( ( Vout  / VOLUME ) * m3_in_liter / year_in_second  )

    call KINXIM_set_NE   ( ( NE  / VOLUME  ) ) ! boite de 1.m3 
    call KINXIM_set_PHIM ( PHIM )
    call KINXIM_set_PRESSURE_TEMPERATURE(PRESSURE, TEMPER)
    call KINXIM_set_QCin ( ( Qcin  / VOLUME )    )
    call KINXIM_set_Vout ( ( Vout  / VOLUME )    )
    
    call KINXIM_set_TAGCELLOUT ( LTAGCELLOUT ) 

    !write(333,*) '--------------------------'
    !write(333,*) ' DT ', DT
    !write(333,*) ' NE  / VOLUME',  NE  / VOLUME
    !write(333,*) ' PHIM', PHIM
    !write(333,*) ' Qcin / VOLUME', Qcin / VOLUME
    !write(333,*) ' Vout / VOLUME', Vout / VOLUME

    !---- CALCUL  -----------------------------------------------------------
    nbtimestep_init = 0

    !call KINXIM_get_nbtimestep ( nbtimestep_init )
    call KINXIM_reinit_stockvar
        
    call KINXIM_CalculDyn ! calcul dynamique chimie

    call KINXIM_get_nbtimestep ( nbtimestep  ) 
    
    !---- RECUPERATION DES VARIABLES ----------------------------------------

    call KINXIM_get_nbtimestep ( nbtimestep  ) 
    call KINXIM_get_ERROR  ( ierror)
    call KINXIM_get_NE     ( NE_DT )
    call KINXIM_get_PHIM   ( PHIM_DT )
    call KINXIM_get_dt     ( dtend )
    call KINXIM_get_NEout  ( NEout )
    call KINXIM_get_ntotaliter  ( niter )
    nbtimestep = nbtimestep - nbtimestep_init 

    !---- DifFERENCE et TRANSformatION PAR HOMOTHETIE INVERSE ---------------
    
!!$    NEout = NEout* VOLUME * liter_in_m3
!!$    dNE(1:nbEAQ) = NE_DT(1:nbEAQ)* VOLUME * liter_in_m3 - NE(1:nbEAQ) 
!!$    dPHIM(1:nbMIN) = PHIM_DT(1:nbMIN) - PHIM(1:nbMIN)

    NEout = NEout* VOLUME 
    dNE(1:nbEAQ) = NE_DT(1:nbEAQ)* VOLUME  - NE(1:nbEAQ) 
    dPHIM(1:nbMIN) = PHIM_DT(1:nbMIN) - PHIM(1:nbMIN)
    
   !! call pause_
    
  end subroutine Compute_API_KINXIM

  !---

  subroutine  Extract_KINXIM(logGammaE, logQsKM, Ionic, RadM, Charge)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: Extract_KINXIM
    ! extraction de variables supplementaires
    ! pour les sorties 
    
    use M_EXCHANGE_KINXIM_BOX
    implicit none
   
    
    real(dp), dimension(:) :: logGammaE
    real(dp), dimension(:) :: logQsKM
    real(dp) :: Ionic
    real(dp), dimension(:) :: RadM
    real(dp) :: Charge
    !---
     call Info_("Extract_KINXIM")
     call KINXIM_get_logGAMMAEAQ(logGammaE)
    call KINXIM_get_logQsKMIN(logQsKM)
    call KINXIM_get_Ionic(Ionic)
    call KINXIM_get_RadM(RadM)
    call KINXIM_get_Charge(Charge)
    !---------------

  end subroutine  Extract_KINXIM


 
end module M_API_KINXIM_BOX
