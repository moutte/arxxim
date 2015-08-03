MODULE M_Solmodel_Calc_Debye_Hueckel
  !------------------------------------------------------------
  ! Debye-Hueckel models of activity for Aqueous Species
  !
  ! DH1    : DebyeHueckel
  ! DH2    : DebyeHueckel + Bdot
  !
  ! DH1EQ3 : DebyeHueckel + EQ3 for neutral species
  ! DH2EQ3 : DebyeHueckel + Bdot + EQ3 for neutral species
  !
  !------------------------------------------------------------

  USE M_Kinds
  USE M_Numeric_Const
  USE M_Trace,ONLY: fTrc,iDebug,T_
  IMPLICIT NONE
  PRIVATE

  PUBLIC:: Solmodel_Calc_Debye_Hueckel

CONTAINS

  SUBROUTINE Solmodel_Calc_Debye_Hueckel( &
  & TdgK, Pbar, vZSp, vMolal,          & !IN
  & dhA, dhB, vA0, bDot, iModel, MWSv, & !IN
  & vLnGam, LnActSv, OsmoSv)

    INTEGER, INTENT(IN) :: vZSp(:)
    REAL(dp),INTENT(IN) :: vMolal(:)
    REAL(dp),INTENT(IN) :: TdgK, Pbar
    REAL(dp),INTENT(IN) :: vA0(:)
    REAL(dp),INTENT(IN) :: dhA
    REAL(dp),INTENT(IN) :: dhB
    REAL(dp),INTENT(IN) :: bDot
    REAL(dp),INTENT(IN) :: MWSv
    INTEGER, INTENT(IN) :: iModel

    REAL(dp),INTENT(OUT):: vLnGam(:)
    REAL(dp),INTENT(OUT):: LnActSv
    REAL(dp),INTENT(OUT):: OsmoSv

    REAL(dp):: Sum_Solute
    REAL(dp):: LnGam_Neutral
    REAL(dp):: A0,X
    REAL(dp):: bDot_Coeff
    REAL(dp):: IonStr, SqrtIoSt
    INTEGER :: vZ2(SIZE(vZSp))
    INTEGER :: nSolut,iAq

    nSolut= SIZE(vMolal)
    Sum_Solute= SUM(vMolal(:))
    vZ2(:)= vZSp(:)*vZSp(:)
    !
    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= DOT_PRODUCT(vZ2(:),vMolal(:)) /Two
    SqrtIoSt= SQRT(IonStr)
    !
    !--------------------------------------------------------- Solute --
    !//======== bDot Models Specific Parameters
    bDot_Coeff= Zero
    !

    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12

    SELECT CASE(iModel)

    CASE(4,5) !("DH2","DH2EQ3")
      bDot_Coeff= bDot

    CASE(10) !("HKF81") ! provisional
      bDot_Coeff= 0.03D0 !bDot
      ! value taken by Dolejs & Wagner,GCA, 2008

    END SELECT
    !
    !//======= Eq3 Models Specific Parameters
    LnGam_Neutral= Zero
    SELECT CASE(iModel)

    CASE(3,5) !("DH1EQ3", "DH2EQ3")
       ! retrieved from EQ36, v7.0, user's manual (Wolery, 1992)
       ! activ'coeff' of CO2(aq) in a (otherwise pure) NaCl solution
       ! as a function of ionic strength
       ! uses an expression proposed by Drummond, 1981
       LnGam_Neutral= (-1.0312D0 + 0.0012806D0*TdgK +255.9D0/TdgK)*IonStr &
       &            - (0.4445D0  - 0.001606D0 *TdgK              )*IonStr/(One+IonStr)
       !from Wolery'92, p41:
       !this should be applied ONLY to non-polar species (O2(aq), H2(aq), N2(aq), ...),
       !for which a salting out effect is expected

    CASE(10) !("HKF81")
      LnGam_Neutral= 0.03D0 *IonStr *Ln10
      ! value taken by Dolejs & Wagner,GCA, 2008

    END SELECT

    DO iAq=1,nSolut
      IF (vZSp(iAq)/=0) THEN
        !//========================================== Charged Species ==
        SELECT CASE(iModel)

        CASE(2,3,4,5) !("DH1","DH1EQ3","DH2","DH2EQ3")
          X= - dhA *vZ2(iAq) *SqrtIoSt/(One+vA0(iAq)*dhB*SqrtIoSt) &
          &  + bDot_Coeff*IonStr
          vLnGam(iAq)= X *Ln10

        CASE(10) !("HKF81")
          A0= 3.72D0
          ! Y=
          X= - dhA*vZ2(iAq)*SqrtIoSt/(One+A0*dhB*SqrtIoSt) &
          &  + LOG(One +MWsv *Sum_Solute) &
          &  + bDot_Coeff*IonStr
          vLnGam(iAq)= X *Ln10

        END SELECT
      ELSE
        !//========================================== Neutral species ==
        vLnGam(iAq)= LnGam_Neutral
      END IF
    END DO
    !--------------------------------------------------------/ Solute --

    !-------------------------------------------------------- solvent --
    !DebyeHueckel + BDot, cf EQ3NR-Doc, p.41
    LnActSv= Ln10 *Two/3.0D0 *dhA *SqrtIost**3 *Sigma(4.0D0*dhB*SqrtIoSt) &
    &      - Ln10            *bDot*IonStr*IonStr &
    &      - Sum_Solute
    !vLnAct(isW)= Zero
    LnActSv=  LnActSv *MWSv
    OsmoSv= - LnActSv /MWsv /Sum_Solute !osmotic coeff'
    !-------------------------------------------------------/ solvent --

  END SUBROUTINE Solmodel_Calc_Debye_Hueckel

  FUNCTION Sigma(x)
    ! for calcul. Water Acivity,
    ! cf equ.23-40 in LewisRandall'61, used also in EQ3NR
    REAL(dp):: Sigma
    REAL(dp),INTENT(IN):: x

    Sigma= 3.0D0 *(One + x - One/(One+x) - Two*LOG(One+x))/x/x/x

  END FUNCTION Sigma

ENDMODULE M_Solmodel_Calc_Debye_Hueckel
