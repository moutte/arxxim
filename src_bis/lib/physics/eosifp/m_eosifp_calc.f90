module M_EosIfp_Calc
!* ======================================================================
!  COORES Simulator      ifP 
!  Translation from Coores Subroutines to a stand alone module by
!  A.Michel, ifP,  2005 - 2008
! -----------------------------------------------------------------------
!  Thermodynamical Properties of EOS Fluid 
!  Models : Cubic Law with binary interactions coefficients
!  Parameters available : PR, SRK
!-----
!  Vectorial Computation 
!  Each point can have a different state (T,P,(Zk)) 
!-----
!  Number of Points     : LVTH
!  Number of Components : NCONS
!-----
!  Results : 
!  1. Covolume :  Z = PV/RT 
!  2. Fugacity Coefficient Logarithm : ln(Phi)
!* ======================================================================
  use M_Kinds,only: dp
  implicit none
  private

  !// private UTIL parameterS
  real(dp), parameter :: ZERO = 0.D0
  real(dp), parameter :: UN   = 1.D0 
  real(dp), parameter :: TIERS = 1.D0 / 3.D0
  real(dp), parameter :: PI = 3.141592654D0
  real(dp), parameter :: DEPIS3 = 2.D0* PI / 3.D0
  !---
  real(dp), parameter :: R = 8.31439D0
  !---
  logical, parameter :: SCAL = .false.
  logical, parameter :: VECT = .true.

  !// public functionS

  public :: EosIfp_VFUGA
  public :: EosIfp_VCUBZ

contains 

!//----------------------------------------------------------------------

subroutine EosIfp_VCUBZ ( &
&   LVTH, NCONS & 
& , T, P, X, Z  &
& , COEFG, COEF1, COEF2, PBI, DLTA1, DLTA2 &
& , PB, PA, GB, AKJGA, CZ1, SAIJXJ, RPB, RBBRTD &
& , CZBAS2,CZALS3,CZ2S3,DISCRI,TRAINV,Z1,Z2,Z3) 
!* ======================================================================
!C/doC/UVCUBZ
!C/1 SIMUSCOPP Simulator.      ifP & Beicip-Franlab.
!C/2 --------------------------------------------- Version 2.0 April 1996
!C/4      RESOLUTION OF THE CUBIC Z-FACTOR EQUATION AND VAPOR/LIQUID
!C/4      IDENTifICATION
!C/4       _____________________________________________________________
!C/4      !  CUBIC EQUATION OF STATE                                    !
!C/4      !  -----------------------                                    !
!C/4      !     P  =  R.T / (V-B) - A(T) / ( (V+DLTA1.B) + (V+DLTA1.B) )!
!C/4      !                                                             !
!C/4      !     DEFINE :                                                !
!C/4      !                                                             !
!C/4      !               GA = A.P / (R.T)**2                           !
!C/4      !               GB = B.P / (R.T)                              !
!C/4      !               Z  = P.V / (R.T)                              !
!C/4      !                                                             !
!C/4      !     THE COMPRESSIBILITY FACTOR Z CAN BE EXPRESSED AS :      !
!C/4      !                                                             !
!C/4      !     Z**3  +  CZ2.Z**2  +  CZ1.Z  +  CZ0  =  0               !
!C/4      !                                                             !
!C/4      !           SD  = DLTA1 + DLTA2                               !
!C/4      !           PD  = DLTA1 * DLTA2                               !
!C/4      !           CZ2 = (SD - 1).BB - 1                             !
!C/4      !           CZ1 = AA - SD.BB + (PD-SD).BB**2                  !
!C/4      !           CZ0 = -(BB.AA + PD.BB.BB.(BB + 1))                !
!C/4      !                                                             !
!C/4      !      FOR MIXTURES, THE parameterS A AND B ARE DEFINED USING !
!C/4      !     THE FOLLOWING                                           !
!C/4      !                                                             !
!C/4      !  MIXING RULE:                                               !
!C/4      !  -----------                                                !
!C/4      !                                                             !
!C/4      !     A = SIGMA_I SIGMA_J ( X_I.X_J.A_IJ )                    !
!C/4      !     B = SIGMA_I ( X_I.B_I )                                 !
!C/4      !                                                             !
!C/4      !                                                             !
!C/4      !                                                             !
!C/4      !  P, V, T : PRESSURE, VOLUME, TEMPERATURE                    !
!C/4      !  A_IJ  = C1_IJ + C2_IJ.SQRT(T) + CG_IJ*T                    !
!C/4      !_____________________________________________________________!
!C/4 ====================================================================
!C/5   PENG D.Y. & ROBINSON D.B. ; 1976 ; "A NEW TWO-CONSTANT EQUATION OF
!C/5         STATE" ; IND. ENG. CHEM. FUNDAM., VOL. 15, NO. 1, PP 59-64
!C/5 --------------------------------------------------------------------
!C/6   COEF1(NCONS,NCONS)  :
!C/6       C1ij = (1-DLTAij).SQRT(PAi.PAj).(1+BTAi)(1+BTAj)
!C/6   COEFG(NCONS,NCONS)  :
!C/6       CGij = (1-DLTAij).SQRT(PAi.PAj).BTAi.BTAj/SQRT(TCi.TCij
!C/6   COEF2(NCONS,NCONS)  :
!C/6       C2ij = -(1-DLTAij).SQRT(PAi.PAj).BTAi.BTAj/SQRT(TCi.TCj)
!C/6              [(SQRT(TCi).(1+BTAi)/BTAi + SQRT(TCj).(1+BTAj)/BTAj ]
!C/6   PBI(NCONS)        : = OMEGA_B*(R*TC)/PC
!C/6   DLTA1,DLTA2        : EQUATION OF STATE parameterS
!C/6 ---------------------
!C/6   P(LVTH) : PRESSURE
!C/6   T(LVTH) : TEMPERATURE
!C/6   X(LVTH,NCONS) : MOLE FRACTION OF THE CONSIDERED PHASE
!C/7 --------------------------------------------------------------------
!C/7   Z(LVTH): PHYSICAL SOLUTION OF THE CUBIC EQUATION
!C/7 ---------------------------------------------------
!C/7   PB(LVTH)                   = SIGMA_i (Bi.Xi)
!C/7   PA(LVTH)                   = SIGMA_i (SIGMA_j (A_ij.Xj).Xi)
!C/7   SAIJXJ(LVTH,NCONS)        = 2*SIGMA_j(Aij.Xj)
!C/7   GB(LVTH)                   = PB.P/(R.T)
!C/7   AKJGA(LVTH)                = PA.P/(R.T)**2
!C/7   CZ1(LVTH)                  = AA - SD.BB + (PD-SD).BB**2
!C/7   RPB(LVTH)    : 1/PB
!C/7   RBBRTD(LVTH) : 1/(R.T.B.B.(DLTA2-DLTA1)
!C/7 ------------------------------------------
!C/8   UCMGPL, UCMGT2
!C/9   ...
!C/11 -------------------------------------------------------------------
!C/11  CZBAS2(LVTH) : CZ0 AND ...  BETA/2
!C/11  CZALS3(LVTH) : -CZ2 AND ... ALPHA/3
!C/11  CZ2S3(LVTH)  : -CZ2/3
!C/11  DISCRI(LVTH) : -DISCRIMINANT ...AND MAX (DISCRIMINANT , 0)
!C/11  TRAINV(LVTH) : logical
!C/11  Z1(LVTH)     : SOLUTIONS OF THE CUBIC EQUATION
!C/11  Z2(LVTH)     :
!C/11  Z3(LVTH)     :
!C/FINdoC/ UVCUBZ /======================================================

  !----------- moduleS ENTIERS TAILLES ------------------------

  implicit none
  integer, intent(in) :: LVTH, NCONS
  real(dp), intent(in) :: COEFG(NCONS,NCONS),COEF1(NCONS,NCONS)
  real(dp), intent(in) :: COEF2(NCONS,NCONS),PBI(NCONS)
  !---
  real(dp), intent(in) :: T(LVTH), P(LVTH), X(LVTH,NCONS)
  real(dp), intent(out) :: Z(LVTH)
  !---
  real(dp) :: PB(LVTH), PA(LVTH), GB(LVTH), AKJGA(LVTH), CZ1(LVTH), SAIJXJ(LVTH,NCONS)
  real(dp) :: RPB(LVTH), RBBRTD(LVTH)
  logical :: TRAINV(LVTH)
  real(dp) :: CZBAS2(LVTH), CZALS3(LVTH), CZ2S3(LVTH), DISCRI(LVTH)
  real(dp) :: Z1(LVTH), Z2(LVTH), Z3(LVTH)

  !----------- FONCTION EXTERNE  ------------------------------------
  !// Integree dans la module f90
  ! integer ULSUM

  !----------- LOCAL VARIABLES
  integer :: I, J, K
  real(dp) :: RE, SQRDIS, RnoneG
  real(dp) :: ALPHA, BETA
  real(dp) :: CUBAL
  real(dp) :: PDLTA, SDLTA, DLTA1, DLTA2 
  real(dp) :: PHI
  real(dp) :: SCALif
  real(dp) :: SI
  real(dp) :: RR2T2, RD2MD1
  real(dp) :: ZMGBL, ZPD2GB , RZPDGB, ZD12BL, ZD12SB 
  real(dp) :: AIJXJ
  real(dp) :: PAIJ11 

  !!$write(*,*) COEF1(1:NCONS, 1:NCONS)
  !!$write(*,*) COEF2(1:NCONS, 1:NCONS)
  !!$write(*,*) COEFG(1:NCONS, 1:NCONS)
  !!$write(*,*) X(1:LVTH, 1:NCONS)
  !!$write(*,*) T(1:LVTH)
  !!$write(*,*) P(1:LVTH)
  !!$write(*,*) PBI(1:NCONS)
  !!$write(*,*) DLTA1,DLTA2

  !C/FONC/       ULSUM
  !*-----------------------------------------------------------------------
  !*     ----------------  COEFFICIENTS  PA, PB, SAIJXJ, GA, GB DE LA PHASE
  do I = 1, LVTH 
    PAIJ11  = COEF1(1,1) + SQRT(T(I))*COEF2(1,1) + COEFG(1,1)*T(I)
    AKJGA(I) = PAIJ11 * X(I,1)
  end do

  !--------------------------------------------------------
  do I = 1, LVTH
    PB(I) = PBI(1) * X(I,1)
    PA(I) = AKJGA(I) * X(I,1)
  end do


  !--------------------------------------------------------
  do  J = 2, NCONS
    do  I = 1, LVTH
      AIJXJ = (COEF1(1,J)+SQRT(T(I))*COEF2(1,J) &
        &            +  COEFG(1,J)*T(I))*X(I,J)
      PA(I) = PA(I) +  AIJXJ * X(I,1)
      AKJGA(I) = AKJGA(I) + AIJXJ
    end do
  end do

  do I = 1, LVTH
    SAIJXJ(I,1) = AKJGA(I) + AKJGA(I)
  end do

  !--------------------------------------------------------

  do  K = 2, NCONS

    do I = 1, LVTH
      AKJGA(I) = (COEF1(K,1) + SQRT(T(I))*COEF2(K,1) &
        &   + COEFG(K,1)*T(I)) * X(I,1)
      PA(I) = PA(I) +  AKJGA(I) * X(I,K)
    end do

    do J = 2, NCONS
      do I = 1, LVTH
        AIJXJ = (COEF1(K,J) + SQRT(T(I))*COEF2(K,J) &
          &  + COEFG(K,J)*T(I)) * X(I,J)
        PA(I) = PA(I) +  AIJXJ * X(I,K)
        AKJGA(I) = AKJGA(I) + AIJXJ
      end do
    end do

    do I = 1, LVTH
      SAIJXJ(I,K) = AKJGA(I) + AKJGA(I)
      PB(I) = PB(I) +  PBI(K) * X(I,K)
    end do

  end do
  do I = 1, LVTH
    RR2T2 = UN / (R * R * T(I) * T(I))
    AKJGA(I) = PA(I) * P(I) * RR2T2
    GB(I) = PB(I) * P(I) / (R * T(I))
  end do

  !-------------------------- *************  -----------------------------
  !--------------------------    CALCUL      -----------------------------
  !-------------------------- *************  -----------------------------

  !*     --------------------------------------- COEFFICIENTS DE LA CUBIQUE
  do I = 1, LVTH
    SDLTA = DLTA1 + DLTA2
    PDLTA = DLTA1 * DLTA2
    CZBAS2(I) = -GB(I) * ( PDLTA * GB(I) * (1.+GB(I)) + AKJGA(I) )
    CZ1(I) = AKJGA(I) -  GB(I) *( SDLTA - (PDLTA - SDLTA) * GB(I) )
    CZALS3(I) = ( SDLTA - 1. ) * GB(I) - 1.
  end do

  !*     -------- CHANGEMENT DE VARIABLE, POUR RAMENER A LA FORME CANONIQUE
  !*              Z  =  X  -  CZALS3 / 3 ,
  !*              D'OU     X**3  +  ALPHA . X  +  BETA  =  0
  do I = 1, LVTH
    CZ2S3(I) = CZALS3(I) * TIERS
    ALPHA = CZ1(I) -  CZALS3(I) * CZ2S3(I)
    CZALS3(I) = ALPHA * TIERS
    BETA = CZBAS2(I) -  CZ2S3(I) * ( CZ2S3(I) * CZ2S3(I)  + ALPHA )
    CZBAS2(I) = BETA * .5
    CUBAL = CZALS3(I) * CZALS3(I) * CZALS3(I)
    !*     -------------- DISCRI = - (LE DISCRIMINANT),  POUR  "UCMGPL" ...
    DISCRI(I) = - CUBAL -  CZBAS2(I) * CZBAS2(I)
  end do

  !*     ----------------------------  EVALUATION DES RACINES DE LA CUBIQUE
  !*                 INITIALISATION DE  Z2,  Z3  SI UNE SEULE RACINE REELLE
  RnoneG = -1.
  do I = 1, LVTH
    Z2(I) = RnoneG
    Z3(I) = RnoneG
  end do

  !*     ------- TRAINV : LOGIQUE, CRITERE POUR MAILLES A 3 RACINES REELLES
  call UCMGPL ( (/.true./), SCAL, (/.false./), SCAL, LVTH, DISCRI, TRAINV )
  !*          ======
  !*         ---  INITIALIS. DE Z1 A L'UNIQUE RACINE REELLE EVENTUELLE
  !*
  do I = 1, LVTH
  !!$  !!write(*,*) "DISCRI=", DISCRI(I), TRAINV(I)
    DISCRI(I) = MAX ( -DISCRI(I), ZERO )
  end do

  !*
  do I = 1, LVTH
    SQRDIS = SQRT ( DISCRI(I) )
    RE = - CZBAS2(I) + SQRDIS
    SI = - CZBAS2(I) - SQRDIS
    RE = SIGN ( ABS ( RE ) ** TIERS, RE )
    SI = SIGN ( ABS ( SI ) ** TIERS, SI )
    Z1(I) = RE + SI - CZ2S3(I)
  end do
  !*
  if ( ULSUM ( LVTH, TRAINV, 1 ) .GT. 0 ) then
    do I = 1, LVTH
      !*        --------- TROIS RACINES REELLES, EVENTUELLEMENT UNE doUBLE...
      if ( TRAINV(I) ) then
        SCALif = SQRT ( -CZALS3(I) )
        PHI = ACOS (  CZBAS2(I) / ( SCALif * CZALS3(I) )  )
        PHI = TIERS * PHI
        SCALif = 2. * SCALif
        Z1(I) = SCALif * COS ( PHI ) -  CZ2S3(I)
        PHI = PHI + DEPIS3
        Z2(I) = SCALif * COS ( PHI ) -  CZ2S3(I)
        PHI = PHI + DEPIS3
        Z3(I) = SCALif * COS ( PHI ) -  CZ2S3(I)
      end if
    end do
  end if

  !*     -------------------------------------- CHOIX DE LA RACINE PHYSIQUE
  do  I = 1 , LVTH
    !*      --- ON ELIMINE LA RACINE INTERMEDIAIRE (INSTABLE), RESTE Z1 < Z2
    Z(I)  = Z1(I)
    Z1(I) = MIN( Z(I), Z2(I), Z3(I) )
    Z2(I) = MAX( Z(I), Z2(I), Z3(I) )
    Z(I)  = Z1(I)
    !*      -------------------------------- ON ELIMINE LA RACINE Z1 SI < GB
    !*                                        IL RESTE GB< Z1 <= Z2
    TRAINV(I) = Z1(I) .GT. GB(I)
  end do

  call UCMGT2 ( Z, .true. , Z2, .true. , LVTH, TRAINV, Z1 )

  RD2MD1     = UN / (DLTA2 - DLTA1)
  do I = 1 , LVTH
    RPB(I)     = UN / PB(I)
  end do

  do I = 1 , LVTH
    RBBRTD(I)  = RPB(I) * RPB(I) * RD2MD1 / ( R * T(I) )
  !!$  write(*,*) "RBBRTD=", RBBRTD(I), RPB(I), RD2MD1, R, T(I)
  end do

  !*    --------------------- CALCUL DE L'ENERGIE LIBRE DE Z1(I) --> -Z(I)
  !*                          CALCUL DE L'ENERGIE LIBRE DE Z2(I) --> -Z3(I)
  do I = 1 , LVTH
    ZMGBL     = log( Z1(I) - GB(I) )
    ZPD2GB    = Z1(I) + DLTA2 * GB(I)
    RZPDGB    = UN / ( Z1(I) +  DLTA1 * GB(I) )
    ZD12BL    = LOG ( ZPD2GB * RZPDGB )
    ZD12SB    = RBBRTD(I) * ZD12BL * PA(I) * PB(I)
    Z(I)      = ZMGBL + ZD12SB - Z1(I)
    !*
    ZMGBL     = log( Z2(I) - GB(I) )
    ZPD2GB    = Z2(I) + DLTA2 * GB(I)
    RZPDGB    = UN / ( Z2(I) +  DLTA1 * GB(I) )
    ZD12BL    = LOG ( ZPD2GB * RZPDGB )
    ZD12SB    = RBBRTD(I) * ZD12BL * PA(I) * PB(I)
    Z3(I)     = ZMGBL + ZD12SB - Z2(I)
    TRAINV(I) = Z3(I) .GT. Z(I)
  end do

  call UCMGT2 ( Z2, .true. , Z1, .true. , LVTH, TRAINV, Z )

  return

end subroutine EosIfp_VCUBZ

!//----------------------------------------------------------------------

subroutine EosIfp_VFUGA ( &
&   LVTH, NCONS & 
& , Z, FUG &
& , PBI, DLTA1, DLTA2 &
& , PB, PA, GB, AKJGA, CZ1, SAIJXJ, RPB, RBBRTD &
& , ZPD2GB, RZPDGB, ZD12BL, ZD12SB, RZMGB, ZM1SB )
!* ======================================================================
!C/doC/UVFUGA
!C/1 SIMUSCOPP Simulator.      ifP & Beicip-Franlab.
!C/2 --------------------------------------------- Version 2.0 April 1996
!C/4   FUGACITY COEFFICIENT CALCULATION
!C/4 ====================================================================
!C/5   PENG D.Y. & ROBINSON D.B. ; 1976 ; "A NEW TWO-CONSTANT EQUATION OF
!C/5         STATE" ; IND. ENG. CHEM. FUNDAM., VOL. 15, NO. 1, PP 59-64
!C/6 --------------------------------------------------------------------
!C/6   Z(LVTH)            :
!C/6 -------------------------
!C/6   PBI(NCONS)          = OMEGA_B*(R*TC)/PC
!C/6   DLTA1,DLTA2         : EQUATION OF STATE parameterS
!C/6 -------------------------
!C/6   PB(LVTH)            = SIGMA_I (BI.XI)
!C/6   PA(LVTH)            = SIGMA_I (SIGMA_J (A_IJ.XJ).XI)
!C/6   SAIJXJ(LVTH,NCONS)  = 2*SIGMA_J(AIJ.XJ)
!C/6   GB(LVTH)            = PB.P/(R.T)
!C/6   AKJGA(LVTH)         = PA.P/(R.T)**2
!C/6   CZ1(LVTH)           = AA - SD.BB + (PD-SD).BB**2
!C/6   RPB(LVTH)           = 1/PB
!C/6   RBBRTD(LVTH)        = 1/(R.T.B.B.(DLTA2-DLTA1)
!C/7 --------------------------------------------------------------------
!C/7   FUG(LVTH,NCONS)     = LOG (FUGACITY COEFFICIENT)
!C/7 -----------------------
!C/7   ZPD2GB(LVTH)        = Z + DLTA2.GB
!C/7   RZPDGB(LVTH)        = 1/(Z + DLTA1.GB)
!C/7   ZD12BL(LVTH)        = log( (Z + DLTA2.B)/(Z + DLTA1.B) )
!C/7   ZD12SB(LVTH)        = RBBRTD.log( (Z + DLTA2.B)/(Z + DLTA1.B) )
!C/7   RZMGB(LVTH)         = 1 /( Z - GB )
!C/7   ZM1SB(LVTH)         = ( Z - 1 )/PB
!C/7 --------------------------------------------------------------------
!C/8
!C/9
!C/FINdoC/ UVFUGA /======================================================

!----------- moduleS ENTIERS TAILLES ------------------------
implicit none

!----------- ARGUMENTS  ------------------------------------
integer :: LVTH, NCONS
real(dp), intent(in) :: Z(LVTH) ,PB(LVTH),PA(LVTH),GB(LVTH),PBI(NCONS) &
  &               ,AKJGA(LVTH),CZ1(LVTH),SAIJXJ(LVTH,NCONS),RPB(LVTH),RBBRTD(LVTH)

real(dp), intent(out) :: FUG(LVTH,NCONS)  ,ZPD2GB(LVTH),RZPDGB(LVTH) &
  &                      ,ZD12BL(LVTH),ZD12SB(LVTH) &
  &                      ,RZMGB(LVTH),ZM1SB(LVTH)

real(dp) :: DLTA1, DLTA2
!-----------

integer :: I,J, K
!*-----------------------------------------------------------------------
do I = 1, LVTH
!!$  write(*,*) "-------------- INPUT = ", i
!!$  write(*,*) Z(I)
!!$  write(*,*) GB(I)
!!$  write(*,*) DLTA2
!!$  write(*,*) DLTA1
!!$  write(*,*) RBBRTD(I)
!!$  write(*,*) RPB(I)

  RZMGB(I)  = UN /( Z(I) - GB(I) )
  ZPD2GB(I) = Z(I) + DLTA2 * GB(I)
  RZPDGB(I) = UN / ( Z(I) +  DLTA1 * GB(I) )
  ZD12BL(I) = LOG ( ZPD2GB(I) * RZPDGB(I) )
  ZD12SB(I) = RBBRTD(I) * ZD12BL(I)
  ZM1SB(I)  = ( Z(I) - UN ) * RPB(I)
  
!!$  write(*,*) "--------------  OUTPUT = ", i
!!$  write(*,*) RZMGB(I)
!!$  write(*,*) ZPD2GB(I)
!!$  write(*,*) RZPDGB(I)
!!$  write(*,*) ZD12BL(I)
!!$  write(*,*) ZD12SB(I)
!!$  write(*,*) ZM1SB(I)
end do
!*     ---------------------------------- LOG. DU COEFFICIENT DE FUGACITE
do K = 1, NCONS
  do I = 1, LVTH
    FUG(I,K) = PBI(K) * ZM1SB(I)  &
    &        + log(RZMGB(I))      &
    &        - ZD12SB(I) * ( SAIJXJ(I,K) * PB(I) - PA(I) * PBI(K) )
  end do
end do

return
end subroutine EosIfp_VFUGA

!//---------------------------------------------------------------------
!      local utils subroutines  
!//---------------------------------------------------------------------

subroutine UCMGT2 ( VAL1, VTSF1, VAL2, VTSF2, LONG, CRITR, RES )

!C/doC/    UCMGT2 /======================================================
!C/1 SIMUSCOPP Simulator.      ifP & Beicip-Franlab.
!C/2 --------------------------------------------- Version 2.0 April 1996
!C/4      UTILITAIRE VECTORIEL,  "UCMGT", 2E VERSION  :      PLUS RAPIDE
!C/4      -                       -----   -
!C/4      SANS LE CALCUL DU NOMBRE DE VALEURS A VRAI DU VECTEUR CRITERE.
!C/5   CRAY RESEARCH, INC. ; "FORTRAN (CFT) REFERENCE MANUAL, SR-0009"
!C/6 ------------------------------------------------------------------
!C/6  VAL1(LONG), VAL2() : VECTEUR(S)/SCALAIRE(S) DE DEPART ( CF. VTSF- )
!C/6  VTSF1/2 : VAL1/2() = VECTEUR LONGUEUR LONG (.T.) OU SCALAIRE (.F.)
!C/6  CRITR(LONG) : VECTEUR DE CHOIX
!C/7  RES(1..LONG) : VECTEUR resultAT, VAL1(I) SI CRITR(I), VAL2(I) SINON
!C/FINdoC/ UCMGT2 /======================================================
implicit none

real(dp) :: VAL1(*), VAL2(*), RES(*)
logical :: VTSF1, VTSF2
logical :: CRITR(*)
integer :: LONG
!------
integer :: I
!*-----------------------------------------------------------------------
!*     ---> programMATION SCALAIRE
!*          ----------------------
if ( VTSF1 .and. VTSF2 ) then
  do I = 1, LONG
    if ( CRITR(I) ) then
      RES(I) = VAL1(I)
    else
      RES(I) = VAL2(I)
    end if
  end do
elseif ( VTSF1 ) then
  do I = 1, LONG
    if ( CRITR(I) ) then
      RES(I) = VAL1(I)
    else
      RES(I) = VAL2(1)
    end if
  end do
elseif ( VTSF2 ) then
  do I = 1, LONG
    if ( CRITR(I) ) then
      RES(I) = VAL1(1)
    else
      RES(I) = VAL2(I)
    end if
  end do
else
  do I = 1, LONG
    if ( CRITR(I) ) then
      RES(I) = VAL1(1)
    else
      RES(I) = VAL2(1)
    end if
  end do
end if

return
end subroutine UCMGT2

!//----------------------------------------------

subroutine UCMGPL( VAL1, VTSF1, VAL2, VTSF2, LONG, CRITR, RES )

!C/doC/    UCMGPL /======================================================
!C/1 SIMUSCOPP Simulator.      ifP & Beicip-Franlab.
!C/2 --------------------------------------------- Version 2.0 April 1996
!C/4      UTILITAIRE VECTORIEL, PROCHE DE LA FONCTION CVMGP ( CRAY )
!C/4      -                                           - ---
!C/5   CRAY RESEARCH INC. ; "FORTRAN (CFT) REFERENCE MANUAL, SR0009".
!C/6 ------------------------------------------------------------------
!C/6  VAL1(), VAL2() : VECTEUR(S)/SCALAIRE(S) DE DEPART, (CF.VTSF1,VTSF2)
!C/6  VTSFI : LOGIQUE, VALI = VECTEUR(1..LONG) SI .T., = SCALAIRE SI .F.
!C/6  CRITR(LONG) : VECTEUR DE CHOIX, REEL
!C/7  RES(1..LONG) : VECT. resultAT, VAL1(I) SI CRITR(I)>0, VAL2(I) SINON
!C/8   CVMGP SI CRAY
!C/FINdoC/ UCMGPL /======================================================
  implicit none
  
  logical:: VAL1(*),VAL2( * ),RES( * )
  real(dp) :: CRITR( * )
  logical :: VTSF1, VTSF2
  integer :: LONG
  
  logical :: VAL1I, VAL2I
  integer :: I
  !*-----------------------------------------------------------------------
  !*
  !*     ---> programMATION SCALAIRE
  !*          ----------------------
  if ( VTSF1 .and. VTSF2 ) then
    do I = 1, LONG
      if ( CRITR(I) .GE. .0 ) then
        RES(I) = VAL1(I)
      else
        RES(I) = VAL2(I)
      end if
    end do
  elseif ( VTSF1 ) then
    VAL2I = VAL2(1)
    do  I = 1, LONG
      if ( CRITR(I) .GE. .0 ) then
        RES(I) = VAL1(I)
      else
        RES(I) = VAL2I
      end if
    end do
  elseif ( VTSF2 ) then
    VAL1I = VAL1(1)
    do I = 1, LONG
      if ( CRITR(I) .GE. .0 ) then
        RES(I) = VAL1I
      else
        RES(I) = VAL2(I)
      end if
    end do
  else
    VAL1I = VAL1(1)
    VAL2I = VAL2(1)
    do I = 1, LONG
      if ( CRITR(I) .GE. .0 ) then
        RES(I) = VAL1I
      else
        RES(I) = VAL2I
      end if
    end do
  end if

  return

end subroutine UCMGPL

!//----------------------------------------------

integer function ULSUM ( LONG, VECT, ISAUT )
!C/doC/    ULSUM  /======================================================
!C/1 SIMUSCOPP Simulator.      ifP & Beicip-Franlab.
!C/2 --------------------------------------------- Version 2.0 April 1996
!C/4      UTILITAIRE VECTORIEL, PROCHE DE LA FONCTION ILSUM  (SCILIB) :
!C/4      -                                            ----
!C/4      RECHERCHE SUR UN VECTEUR LOGIQUE DU NOMBRE DE COMPOSANTE(S)
!C/4      AYANT LA VALEUR  "VRAI".
!C/5   CRAY RESEARCH, INC. ; "LIBRARY REFERENCE MANUAL, SR-0014" ;
!C/5                                     (doCUMENTATION  "ILSUM", SCILIB)
!C/6 ------------------------------------------------------------------
!C/6  VECT(1,1..LONG) : VECTEUR A EXPLORER
!C/6  ISAUT : PAS EVENTUEL DE LA RECHERCHE, TRADUIT COMME PREMIERE
!C/6       :                                  dimension DECLAREE DE  VECT
!C/7  ULSUM : NOMBRE DE COMPOSANTES A  "VRAI"
!C/FINdoC/ ULSUM  /======================================================
!
  implicit none
  
  integer :: ISAUT, LONG
  logical :: VECT(ISAUT,*)
  !---
  integer :: I
  !*-----------------------------------------------------------------------
  !*
  !*     ---> programMATION SCALAIRE
  !*          ----------------------
  ULSUM = 0
  do I = 1, LONG
    if ( VECT(1,I) )   ULSUM = ULSUM + 1
  end do

  return
end function ULSUM

end module M_EosIfp_Calc
