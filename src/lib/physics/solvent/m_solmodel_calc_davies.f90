MODULE M_Solmodel_Calc_Davies
!-----------------------------------------------------------------------
! Davies models of activity for Aqueous Species
!
! DAV_1 :  Davies 1
!   after Grenthe et al, Modelling In Aquatic Chemistry, OECD,Chap.IX,p.330
!   (same as in Bethke,'96, p110)
!   Davies, with Gamma(Neutral)=1.0
!
! DAV_2 :  Davies 2
!   after Grenthe et al, Modelling In Aquatic Chemistry, OECD,p.330
!   Davies, with Gamma(Neutral) a FUNCTION of IonStrength
!
! SAMSON : Samson et al. 1999
!   Modeling chemical activity effects in strong ionic solutions.
!   Computational Materials Science 15: 285-294
!
!-----------------------------------------------------------------------

  USE M_Numeric_Const
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,iDebug,T_

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: Solmodel_Calc_Davies

CONTAINS

  SUBROUTINE Solmodel_Calc_Davies( &
  & vZSp, vMolal, &
  & dhA, iModel, MWSv, &
  & vLnGam, LnActSv, OsmoSv) !,ERR)

    INTEGER, INTENT(IN)::  vZSp(:)
    REAL(dp),INTENT(IN)::  vMolal(:)
    REAL(dp),INTENT(IN)::  dhA
    REAL(dp),INTENT(IN)::  MWSv
    INTEGER, INTENT(IN)::  iModel

    REAL(dp),INTENT(OUT):: vLnGam(:)
    REAL(dp),INTENT(OUT):: LnActSv
    REAL(dp),INTENT(OUT):: OsmoSv

    REAL(dp):: LnGam_Charged, LnGam_Neutral
    REAL(dp):: Sum_Solute
    REAL(dp):: IonStr, SqrtIoSt
    INTEGER :: vZ2(SIZE(vZSp))
    INTEGER :: nSolut,iAq
    !-------------------------------------------------------------------

    nSolut= SIZE(vMolal)
    vZ2(:)= vZSp(:)*vZSp(:)

    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= DOT_PRODUCT(vZ2(:),vMolal(:)) /Two

    !-- "TRUNCATED DAVIES"
    !-- IF(IonStr>0.3) IonStr- 0.3D0

    SqrtIoSt= SQRT(IonStr)

    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12

    !---------------------------------------------------------< Solute--
    SELECT CASE(iModel)
    !
    CASE(6,7) !"DAV_1","DAV_2")
      LnGam_Charged = -dhA &
      &               *( SqrtIoSt/(One + SqrtIoSt) -0.3D0*IonStr ) &
      &               *Ln10
      ! there's also a version of Davies equation
      ! where coeff 0.2 is used instead of 0.3 !?!
    !
    CASE(9) !("SAMSON")
      LnGam_Charged = -dhA &
      &               *( SqrtIoSt/(One + SqrtIoSt) &
      &                 -(0.2D0 -4.17D-5*IonStr) *IonStr ) &
      &               *Ln10
    !
    END SELECT
    
    SELECT CASE(iModel)
    !
    CASE(6,9) !("DAV_1","SAMSON")
      LnGam_Neutral = Zero
      !-> official Davies equation
    !
    CASE(7) !"DAV_2")
      LnGam_Neutral = 0.064D0*IonStr*Ln10
      !-> non-official equation (where is it from ??)
    !
    END SELECT
    
    DO iAq=1,nSolut
      IF (vZSp(iAq)/=0) THEN
        !------------------Charged Species--
        ! e.g. Bethke, eq.7.4, p110
        vLnGam(iAq)= vZ2(iAq) *LnGam_Charged
      ELSE
        !------------------Neutral Species--
        vLnGam(iAq)= LnGam_Neutral
      ENDIF
    ENDDO
    !---------------------------------------------------------</Solute--
    
    !-------------------------------------------------------- Solvent --
    Sum_Solute = SUM(vMolal(:))
    !___________LewisRandall,23-39,p347; EQ3NR-Doc, p.39
    LnActSv= Ln10 *Two/3.0D0 *dhA  *SqrtIost**3 *Sigma(SqrtIoSt) &
    &      - Ln10            *dhA  *IonStr*IonStr &
    &      - Sum_Solute
    LnActSv=  LnActSv * MWSv
    OsmoSv= - LnActSv /MWsv / Sum_Solute
    !--------------------------------------------------------/Solvent --

  END SUBROUTINE Solmodel_Calc_Davies

  REAL(dp) FUNCTION Sigma(x)
  ! for calcul. Water Acivity,
  ! cf equ.23-40 in LewisRandall'61, used also in EQ3NR
    REAL(dp),INTENT(IN):: x
    
    Sigma=  3.0D0*(One + x - One/(One+x) - Two*LOG(One+x))/x/x/x
  
  END FUNCTION Sigma

ENDMODULE M_Solmodel_Calc_Davies
