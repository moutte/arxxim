module M_exchange_KINXIM_BOX

  use M_string
  use M_Kinds
  implicit none
  private

  !-------- STATIQUES  !-- GET ( KINXIM -> COORES)

  !-- parametres
  public :: KINXIM_get_nbC
  public :: KINXIM_get_nbEAQ
  public :: KINXIM_get_nbMIN

  !-- noms
  public :: KINXIM_get_nomC
  public :: KINXIM_get_nomEAQ
  public :: KINXIM_get_nomMIN

  !-- masse-molaire
  public :: KINXIM_get_massmolC

  ! species properties
  public :: KINXIM_get_a0E
  public :: KINXIM_get_ZE
  public :: KINXIM_get_KE

  ! rock properties
  public :: KINXIM_get_KM
  public :: KINXIM_get_Radm
  public :: KINXIM_get_Srm

  !-- matrices de transformation

  public :: KINXIM_get_EAQtoC
  public :: KINXIM_get_MINtoC

  !-- densite molaire
  public :: KINXIM_get_RHOM
  public :: KINXIM_set_RHOM

  !-------- VARIABLES !-- GET et SET ( KINXIMMEDXE <-> COORES)

  !-- variables

  public :: KINXIM_get_NC
  public :: KINXIM_get_NE
  public :: KINXIM_get_PHIM

  public :: KINXIM_set_NE
  public :: KINXIM_set_Vout
  public :: KINXIM_set_Qcin
  public :: KINXIM_set_PHIM
  public :: KINXIM_set_PRESSURE_TEMPERATURE

  !-- gestion du temps
  public :: KINXIM_get_t
  public :: KINXIM_get_tmore
  public :: KINXIM_get_dt

  public :: KINXIM_set_t
  public :: KINXIM_set_dt
  public :: KINXIM_set_tmore
  public :: KINXIM_get_NEout

  !-- controle vecteur de calcul
  public :: KINXIM_get_X

  !------------------------------------------
  ! variables secondaires pour output
  !------------------------------------------
  public :: KINXIM_get_logGAMMAEAQ  ! activite ( var secondaire affichage )
  public :: KINXIM_get_logQsKMIN    ! activite ( var secondaire affichage )
  public :: KINXIM_get_Ionic        ! activite ( var secondaire affichage )
  public :: KINXIM_get_Charge       ! Charge   ( var secondaire affichage )
  public :: KINXIM_get_contextMIN   ! Context Mineral ( Precip,Dissol,Test)
  public :: KINXIM_get_SALINITY     ! Salinity computed with equiv Cl
  ! -----------------------------------------

  public :: KINXIM_get_ERROR
  public :: KINXIM_get_niter
  public :: KINXIM_get_ntotaliter
  public :: KINXIM_get_nbtimestep

  public :: KINXIM_reinit_stockvar  ! reinit stockvar

  !------------------------------------------

  public :: KINXIM_set_TAGCELLOUT   ! set dynamic-ouput on/off
  public :: KINXIM_set_TAGCELLHIST  ! set dynamic-log on/off

  !--------------------------------------------

contains

  !----------------------------------------------------------------------------
  !***********  ELEMENTS ET STOECHIOMETRIE
  !----------------------------------------------------------------------------

  !*************------------------------------

  subroutine  KINXIM_get_nbC(xnbC)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nbC
    use M_common_KINXIM_BOX
    implicit none
    integer :: xnbC
    !--
    xnbC = nCp

  end subroutine KINXIM_get_nbC

  !---

  subroutine  KINXIM_get_ZE(xZE)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_ZE
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xZE
    integer :: iAq, iSpc
    !---
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       xZE(iAq) = vSpc(iSpc)%Z
    end do
  end subroutine KINXIM_get_ZE

  !---

   subroutine  KINXIM_get_Srm(xSrm )
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_Srm
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xSrm
    real(dp) :: Pi
    !---
    Pi = atan(1.d0)*4.D0
    xSrm = 0.D0
    
  end subroutine KINXIM_get_Srm

  !---

  subroutine  KINXIM_get_Radm(xRadm )
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_Radm
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xRadm
    !---
    xRadm = 0.D0

  end subroutine KINXIM_get_Radm

  !---

  subroutine  KINXIM_get_KM(xKM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_KM
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xKM
    integer :: xnbMIN
    real(dp) :: lnK
    !---
    xKM = 1.D0
    
  end subroutine KINXIM_get_KM

  !---

  subroutine  KINXIM_get_KE(xKE)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_KE
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xKE
    integer :: xnbEAQ
    !---
    xKE = 1.D0
    
  end subroutine KINXIM_get_KE

  !---

  subroutine  KINXIM_get_EAQtoC(xEAQtoC, xnbEAQ, xnbC)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_EAQtoC

    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:,:) :: xEAQtoC
    integer :: xnbC
    integer :: xnbEAQ
    integer :: iAq, iC
    !--- 
    ! nCp = xnbC
    ! nAq = xnbEAQ
    !--
    do iAq=1, nAq
       do iC = 1, nCp
          xEAQtoC(iAq, iC) = tAlfAq(iC, iAq)
       end do
    end do

  end subroutine KINXIM_get_EAQtoC


  !*************------------------------------

  subroutine  KINXIM_get_MINtoC(xMINtoC, xnbMIN, xnbC)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_MINtoC

    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:,:) :: xMINtoC
    integer :: xnbMIN
    integer :: xnbC
    integer :: iC, iMk
    !--- 
    ! nCp = xnbC
    ! nMk = xnbMIN
    !--
    do iMk = 1, nMk
       do iC = 1, nCp
          xMINtoC(iMk, iC) = tAlfMk(iC, iMk)
       end do
    end do

  end subroutine KINXIM_get_MINtoC

  !----------------------------------------------------------------------------
  !***********  DEBITS
  !----------------------------------------------------------------------------

  subroutine KINXIM_set_Vout(xVout)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_Vout
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xVout
    !---

    LCouplerOutflow = .true.
    Vout = xVout
    !write(*,*) 'I SET Vout = ', Vout , xVout

  end subroutine KINXIM_set_Vout

  !*************------------------------------

  subroutine  KINXIM_set_Qcin(xQcin)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_Qcin
    use M_echomat
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xQcin
    integer :: xnbC
    logical :: LDEBUG_KINXIM_INPUT = .false.
    character(len=3), allocatable :: xnomC(:)
    !---
    !write(*,*) "I Set Qcin"
    
    LCouplerInject = .true. ! active coupler source mode
    vQInj(1:nCp) = -xQcin(1:nCp)
    
    !----// Debug
    if (LDEBUG_KINXIM_INPUT) then
       write(*,*) '------- Input Source Terms'
       write(*,*) "xQcin", xQcin
       write(*,*) "nCp=", nCp
       !-
       call echovec_real( xQcin, nCp, 'QTotInj', TRTAG='T', unit=6)
       call echovec_real( vQInj, nCp, 'vQInj',   TRTAG='T', unit=6)
    end if
  end subroutine KINXIM_set_Qcin


  !----------------------------------------------------------------------------
  !***********  WATER
  !----------------------------------------------------------------------------

  subroutine  KINXIM_get_nbEAQ(xnbEAQ)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nbEAQ
    use M_common_KINXIM_BOX
    implicit none
    integer :: xnbEAQ
    !---
    xnbEAQ =  nAq

  end subroutine KINXIM_get_nbEAQ

  !---

  subroutine  KINXIM_get_NC(xNC, xnbEAQ, xnbC) ! NC dans la phase aqueuse !
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_NC
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xNC
    integer :: xnbC, xnbEAQ
    real(dp), dimension(xnbEAQ) :: xNE
    real(dp), dimension(xnbEAQ,xnbC) :: xEAQtoC
    !---
    ! stoechiometrie + calcul a partir de xnE
    call KINXIM_get_NE(xNE)
    call KINXIM_get_EAQtoC(xEAQtoC, xnbEAQ, xnbC )
    xNC(1:xnbC) = matmul(transpose(xEAQtoC(1:xnbEAQ, 1:xnbC)),xNE(1:xnbEAQ))

  end subroutine KINXIM_get_NC

  !*************------------------------------

  subroutine KINXIM_get_nomEAQ(xnom)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nomEAQ
    use M_common_KINXIM_BOX
    implicit none
    character(len=*), dimension(:) :: xnom
    integer :: xnbEAQ, iAq, iSpc
    !---
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       xnom(iAq) = vSpc(iSpc)%NamSp
       call string_lowcase(  xnom(iAq) )
    end do

    !// debug
    !write(*,*) '------'
    !do iAq=1, nAq
    !   write(*,'(A,I4,A,A)') "nomEAQ [ ", iAq," ] = ", xnom(iAq)
    !end do
    
  end subroutine KINXIM_get_nomEAQ

  !*************------------------------------

  subroutine KINXIM_get_nomC(xnom)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nomC
    use M_common_KINXIM_BOX
    implicit none
    character(len=*), dimension(:) :: xnom
    integer :: iC
    !---
    do iC = 1, nCp
       xnom(iC) = vCpn(iC)%NamCp
       call string_lowcase( xnom(iC) )
       call string_subst ( xnom(iC), '_', ' ')
    end do
    
    !// debug
    !write(*,*) '------'
    !do iC=1, nCp
    !   write(*,'(A,I4,A,A)') "nomC   [ ", iC," ] = ", xnom(iC)
    !end do
       
  end subroutine KINXIM_get_nomC

  !*************------------------------------

  subroutine KINXIM_get_massmolC(xmassmolC)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_massmolC
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xmassmolC
    integer :: iC
    !---
    do iC = 1, nCp
       xmassmolC(iC) = vEle(vCpn(iC)%iEle)%WeitKg
    end do
  end subroutine KINXIM_get_massmolC

  !*************------------------------------

  subroutine  KINXIM_set_NE(xNE)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_NE
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xNE
    integer :: iAq
    !---
    do iAq = 1, nAq
       vMolF(iAq) = xNE(iAq)
    end do

  end subroutine KINXIM_set_NE

  !---

  subroutine  KINXIM_get_NE(xNE)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_NE
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xNE
    integer :: iAq
    !---    
    do iAq = 1, nAq
        xNE(iAq) = vMolF(iAq)
    end do

  end subroutine KINXIM_get_NE

  !--------

  subroutine  KINXIM_get_a0E(xa0E)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_a0E
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xa0E
    integer :: iAq, iSpc
    !---
    
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       xa0E = vSpc(iSpc)%AquSize
    end do

  end subroutine KINXIM_get_a0E

  !--------

  subroutine  KINXIM_get_logGAMMAEAQ(xlogGAMMA)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_logGAMMAEAQ
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xlogGAMMA
    integer :: xnbEAQ, iAq
    real(dp) :: log10
    !---
    log10 = log(10.D0)

    do iAq = 1, nAq
       xlogGAMMA(iAq) = vLnGam(iAq)/log10
    end do


  end subroutine KINXIM_get_logGAMMAEAQ


  !----------------------------------------------------------------------------
  !***********  MINERAUX
  !----------------------------------------------------------------------------

  subroutine  KINXIM_get_nbMIN(xnbMIN)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nbMIN
    use M_common_KINXIM_BOX
    implicit none
    integer :: xnbMIN
    !---
    xnbMIN = nMk 

  end subroutine KINXIM_get_nbMIN

  !*************------------------------------

  subroutine KINXIM_get_nomMIN(xnom)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nomMIN
    use M_common_KINXIM_BOX
    implicit none
    character(len=*), dimension(:) :: xnom
    integer :: iSpc, iMk
    !--
    do iMk = 1,nMk
       iSpc = vOrdMk(iMk)
       xnom(iMk) = vKinFas(iMk)%NamKF
       !xnom(iMk) = vSpc(iSpc)%Name
       call string_lowcase( xnom(iMk) )
    end do

    !// debug
    ! write(*,*) '------'
    !do iMk=1, nMk
    !   write(*,'(A,I4,A,A)') "nomMIN [ ", iMk," ] = ", xnom(iMk)
    !end do

  end subroutine KINXIM_get_nomMIN

   !*************------------------------------


  subroutine  KINXIM_get_PHIM(xPHIM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_PHIM
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xPHIM
    integer :: iMk, iSpc
    real(dp) :: MassPhase, RhoPhase, PhiPhase
    !--
    !write(*,*) "I GET PHIM !!"
    do iMk = 1, nMk
       iSpc = vOrdMk(iMk)
       MassPhase = vMolK(iMk)*vSpc(iSpc)%WeitKg
       RhoPhase  = vRhoK(iMk)
       !write(*,*) "--------- iMk =", iMk
       !write(*,*) vRhoK(iMk)
       !write(*,*) "vMolK=", vMolK(iMk)
       !write(*,*) Vbox
       PhiPhase  = MassPhase/RhoPhase/VBox
       xPHIM(iMk) = PhiPhase  
       !write(*,*) "PhiPhase =", iMk, PhiPhase
       !write(*,*) "--------- "
    end do

  end subroutine KINXIM_get_PHIM


  !*************------------------------------

  subroutine  KINXIM_get_logQsKMIN(xlogQKMIN)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_logQsKMIN
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xlogQKMIN
    integer :: iMk
    real(dp) :: Log10
    !---
    do iMk = 1, nMk
       xlogQKMIN(iMk) = log10(vOmegaK(iMk))
    end do
  end subroutine KINXIM_get_logQsKMIN

  !*************------------------------------

  subroutine  KINXIM_get_Ionic(xIonic)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_Ionic
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xIonic
    !---
    xIonic = Ionic

  end subroutine KINXIM_get_Ionic

  !*************------------------------------

  subroutine  KINXIM_get_RHOM(xRHOM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_RHOM
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xRHOM
    integer :: iMk, iSpc
    !---
    do iMk = 1, nMk
       iSpc = vOrdMk(iMk)
       xRHOM(iMk) = vRhoK(iMk)/vSpc(iSpc)%WeitKg
    end do

  end subroutine KINXIM_get_RHOM

  !---

  subroutine  KINXIM_set_RHOM(xRHOM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_RHOM
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xRHOM
    !---
    write(*,*) "TOdo : KINXIM_set_RhoM."
    write(*,*) "QUESTION : Is it possible to modify mineral molar volumes in Arxim v3 ?"

  end subroutine

  !---

  subroutine  KINXIM_set_PHIM(xPHIM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_PHIM
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xPHIM
    integer :: iMk, iSpc
    real(dp) :: PhiPhase, MassPhase 
    !---
    !write(*,*) "I SET PHIM !!"
    do iMk = 1, nMk
       iSpc = vOrdMk(iMk)
       PhiPhase  = xPHIM(iMk)
       MassPhase = PhiPhase*VBox*vRhoK(iMk)
       vMolK(iMk) = MassPhase/vSpc(iSpc)%WeitKg
       !write(*,*) "--------- iMk =", iMk
       !write(*,*) vRhoK(iMk)
       !write(*,*) PhiPhase, MassPhase
       !write(*,*) "vMolK=", vMolK(iMk)
    end do
    
  end subroutine KINXIM_set_PHIM

  !---

   !*************------------------------------


  subroutine  KINXIM_set_PRESSURE_TEMPERATURE(xPRESSURE, xTEMPER)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_PRESSURE_TEMPERATURE
    use M_common_KINXIM_BOX
    implicit none
    real(dp), intent(in) :: xPRESSURE, XTEMPER

    ! set (T,P) conditions

    !write(*,*) 'KINXIM I set (T,P)=', xPRESSURE, xTEMPER
    TdgKBox = xTEMPER
    PbarBox = xPRESSURE*1.D-5

  end subroutine KINXIM_set_PRESSURE_TEMPERATURE

  !---

  subroutine KINXIM_get_contextMIN(xcontextM)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_contextMIN
    use M_common_KINXIM_BOX
    implicit none
    character, dimension(:) :: xcontextM
    integer :: iMk
    !--
    do iMk = 1, nMk
       select case(trim(vStatusK(iMk)))
       case("DISSOLU") ; xcontextM(iMk)="D"             
       case("PRECIPI") ; xcontextM(iMk)="P"
       case("MINIMAL") ; xcontextM(iMk)="M"
       case("INERT")   ; xcontextM(iMk)="I"
       case default
          xcontextM(iMk)="X"
          write(*,*) "Error context unknown :", vStatusK(iMk)
          stop "wrong context detected"          
       end select
    end do
  end subroutine KINXIM_get_contextMIN

  !*************------------------------------

  subroutine KINXIM_get_X(xX)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_X
    use M_common_KINXIM_BOX
    implicit none
    real(dp), dimension(:) :: xX
    integer :: nbX
    integer :: iX, iMk, iAq
    !---
    nbX = nAq + nMk
    
    do iAq = 1,nAq
       iX = iAq
       xX(iX) =  vMolF(iAq)
    end do
    
    do iMk = 1,nMk
      iX = nAq + iMk
      xX(iMk) =  vMolK(iMk)
    end do
    
  end subroutine KINXIM_get_X


  !*************------------------------------

  subroutine KINXIM_get_t(xt)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_t
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xt
    !----
    xt = Time 

  end subroutine KINXIM_get_t
  
  
  subroutine KINXIM_get_tmore(xt)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_t
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xt
    !----
    xt = Tfinal - Time  

  end subroutine KINXIM_get_tmore

  !*************------------------------------

  subroutine KINXIM_get_dt(xdt)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_dt
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xdt
    !----
    !write(*,*) "I Get DT", dTime
    xdt = dTime

  end subroutine KINXIM_get_dt

!---

  subroutine KINXIM_set_tmore(xtmore)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_tmore
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xtmore
    !---
    Tfinal = Time + xtmore
   ! write(*,*) "I Set Tfinal", Tfinal

  end subroutine KINXIM_set_tmore

  !*************------------------------------

  subroutine KINXIM_set_t(xt)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_t
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xt 
    !---
    Time = xt
    
  end subroutine KINXIM_set_t

  !*************------------------------------

  subroutine KINXIM_set_dt(xdt)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_dt
    use M_common_KINXIM_BOX
    !---
    implicit none
    real(dp) :: xdt
    !---
    !write(*,*) "I Set DT", xdt
    
    DTime = xdt
    
  end subroutine KINXIM_set_dt

  !*************------------------------------

  subroutine KINXIM_get_ERROR(xierror)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_ERROR
    use M_common_KINXIM_BOX
    implicit none
    integer  :: xierror
    !---
    if (TimeLoop_Succes) then
       xierror = IerrorChemistry
    else
       xierror = 1
    end if
    
  end subroutine KINXIM_get_ERROR

  !*************------------------------------

  subroutine KINXIM_get_ntotaliter(xntotaliter)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_ntotaliter
    use M_common_KINXIM_BOX
    implicit none
    integer  :: xntotaliter
    !---
    xntotaliter = Total_NewtonIter

  end subroutine KINXIM_get_ntotaliter

  !*************------------------------------

  subroutine KINXIM_get_niter(xniter)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_niter
    use M_common_KINXIM_BOX
    implicit none
    integer  :: xniter
    !---
    xniter = NewtonIter

  end subroutine KINXIM_get_niter

  !*************------------------------------

  subroutine KINXIM_get_nbtimestep(xnbtimestep)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_nbtimestep
    use M_common_KINXIM_BOX
    implicit none
    integer  :: xnbtimestep
    !---
    xnbtimestep = iTimeStep

  end subroutine KINXIM_get_nbtimestep

 !*************------------------------------

  subroutine KINXIM_set_TAGCELLOUT(lval)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_TAGCELLOUT
    use M_common_KINXIM_BOX
    implicit none
    logical :: lval
    !---
    !!DebugCoores = lval !!! a verifier

  end subroutine KINXIM_set_TAGCELLOUT

    !*************------------------------------

  subroutine KINXIM_set_TAGCELLHIST(lval)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_set_TAGCELLHIST
    use M_common_KINXIM_BOX
    implicit none
    logical :: lval
    !---
    !!iDebug = 0
    !!    if (lval) iDebug = 1

  end subroutine KINXIM_set_TAGCELLHIST


   !*************------------------------------

   !! recuperation du nombre de moles moyen en sortie (mol)

  subroutine  KINXIM_get_NEout(xNEout)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_NEout
    use M_common_KINXIM_BOX
    use M_stockvar_KINXIM_BOX, only : size_X, computemean_stockvar, LSTOCK, print_stockvar
    use M_echomat

    implicit none

    ! -- arguments
    real(dp), dimension(:) :: xNEout
    integer :: xnbeaq

    !--- local variables
    real(dp), dimension(:), allocatable :: Xout
    real(dp) :: DT, dtnu, bilanERR, xVout !, tol
    logical :: LDEBUG_KINXIMNEOUT = .false.

    !--- instructions
    !write(666,*) "LSTOCK= ", LSTOCK
    !write(666,*) "size_X = ", size_X
    
    allocate( Xout (size_X) )
    !!call print_stockvar("stocky.var")
    call computemean_stockvar(Xout, xVout, DT, dtnu, bilanERR)

    !-- extraction des coordonnees especes
    xnbeaq = nAq
    xNEout(1:nAq) = Xout (1:nAq)

    !-- debug
    if ( LDEBUG_KINXIMNEOUT) then
       write(*,*) "-------- BILAN Flux KINXIM out ----------------------------"
       write(*,*) "DT       = ", DT
       write(*,*) "dt moyen = ", dtnu
       write(*,*) "Vout     = ", xVout
       write(*,*) "bilanERR = ", bilanERR
       call echovec_real(XNEout, xnbeaq, "XEout", trtag='T')
       call echovec_real(XNEout * xVout, xnbeaq, "FluxEout", trtag='T')
       call echovec_real(XNEout * xVout * DT, xnbeaq, "Total Moles sorties", trtag='T')
       write(*,*) "---------------------------------------------------------"
     end if

    !-- deallocation
    deallocate( Xout )

  end subroutine KINXIM_get_NEout

  !*************------------------------------

  subroutine  KINXIM_reinit_stockvar
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_reinit_stockvar
    use M_stockvar_KINXIM_BOX
    implicit none
    !-- instructions

    call reset_stockvar

  end subroutine KINXIM_reinit_stockvar

   !-------

  subroutine KINXIM_get_SALINITY(xsalinity)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_SALINITY
    use M_common_KINXIM_BOX
    implicit none
    real(dp) xsalinity
    character(len=10) :: NAME(10)
    real(dp) :: VAlenCY(10)
    real(dp) :: MwtNa = 23
    real(dp) :: MwtCl = 35.5
    integer :: xnbC, xnbEAQ
    real(dp), allocatable :: nC_(:)
    character(len=10), allocatable :: nomC_(:)
    integer :: i, j
    !------------------------------
    ! data BASE
    !------------------------------
    NAME(1)  = 'na'
    NAME(2)  = 'ca'
    NAME(3)  = 'sr'
    NAME(4)  = 'ba'
    NAME(5)  = 'mg'
    NAME(6)  = 'k'
    NAME(7)  = 'li'
    NAME(8)  = 'mn'
    NAME(9)  = 'fe'
    NAME(10) = 'so4'
    !-------------
    VAlenCY(1) = 1.
    VAlenCY(2) = 2.
    VAlenCY(3) = 2
    VAlenCY(4) = 2.
    VAlenCY(5) = 2.
    VAlenCY(6) = 1.
    VAlenCY(7) = 1.
    VAlenCY(8) = 2.
    VAlenCY(9) = 0.
    VAlenCY(10)= 0.
    !------------------------------
    call KINXIM_get_nbC(xnbC)
    call KINXIM_get_nbEAQ(xnbEAQ)
    allocate(nomC_(1:xnbC))
    allocate(nC_(1:xnbC))

    call KINXIM_get_nomC(nomc_(1:xnbC))
    call KINXIM_get_NC(nC_(1:xnbC), xnbEAQ, xnbC)
    xSalinity = 0.D0

    do i=1, 10
       do j=1, xnbC
          if ( NAME(i).eq.nomC_(j) ) then
             !write(*,*) 'matching ', i,',', j, ' | ', NAME(i),'=', nomC_(j)
             xSalinity = xSalinity + nC_(j)*VAlenCY(i)
          end if
       end do
    end do
    !write(*,*)  xSalinity
    xSalinity = xSalinity * (MwtNa+MwtCl)
    deallocate(nomC_, nC_)

  end subroutine KINXIM_get_SALINITY

  !----

  subroutine  KINXIM_get_Charge(xCharge)
!lt !DEC$ ATTRIBUTES DLLEXPORT :: KINXIM_get_Charge
    use M_common_KINXIM_BOX
    implicit none
    real(dp) :: xCharge
    !real(dp) :: Ztmp
    integer :: i !, js
    real(dp), allocatable :: ChargeE(:)
    real(dp), allocatable :: NE(:)
    integer :: nbEAQ
    !---
    call KINXIM_get_nbEAQ(nbEAQ)

    allocate(NE(nbEAQ))
    allocate(ChargeE(nbEAQ))
    call KINXIM_get_NE(NE)
    call KINXIM_get_ZE(ChargeE)

    xCharge = 0.D0
    do i=1,nbEAQ
       xCharge =  xCharge + NE(i)*ChargeE(i)
    end do

  end subroutine KINXIM_get_Charge

  !----------------------------------------------------------------------------
  !---                  EXTENSIONS POSSIBLES GAZ et SOLUTIONS SOLIDES       ---
  !---                        A REDEFINIR                                   ---
  !----------------------------------------------------------------------------

  end module M_exchange_KINXIM_BOX
