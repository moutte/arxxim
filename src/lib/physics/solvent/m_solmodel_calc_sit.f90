MODULE M_Solmodel_Calc_SIT
  !---------------------------------------------------------------
  ! SIT model of activity for Aqueous Species
  !          after Grenthe et al, ModellingInAquaticChemistry,OECD,Chap.IX,p.330
  !----------------------------------------------------------------

  USE M_Numeric_Const
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,iDebug,T_
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC:: Solmodel_Calc_SIT

CONTAINS

  SUBROUTINE Solmodel_Calc_SIT( &
  & vZSp, vMolal, &
  & dhA, Model, MWSv, &
  & vLnGam, LnActSv, OsmoSv)
  
    INTEGER, INTENT(IN)     :: vZSp(:)
    REAL(dp),INTENT(IN)     :: vMolal(:)
    REAL(dp),INTENT(IN)     :: dhA
    REAL(dp),INTENT(IN)     :: MWSv
    CHARACTER(*),INTENT(IN) :: Model
    REAL(dp),INTENT(OUT)    :: vLnGam(:)
    REAL(dp),INTENT(OUT)    :: LnActSv
    REAL(dp),INTENT(OUT)    :: OsmoSv
    !
    REAL(dp):: LnGam_Charged, LnGam_Neutral
    REAL(dp):: Sum_Solute
    REAL(dp):: IonStr, SqrtIoSt, X
    INTEGER :: vZ2(SIZE(vZSp))
    INTEGER :: nSolut,iAq
    
    nSolut= SIZE(vMolal)
    vZ2(:)= vZSp(:)*vZSp(:)
    !    
    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= DOT_PRODUCT(vZ2(:),vMolal(:)) /Two
    !
    SqrtIoSt= SQRT(IonStr)
    !
    !--------------------------------------------------------- Solute --
    LnGam_Charged = -dhA *(SqrtIoSt/(One + 1.5D0*SqrtIoSt))*Ln10
    !
    LnGam_Neutral = Zero
    !
    DO iAq=1,nSolut
      IF (vZSp(iAq)/=0) THEN
        !--- Charged Species --!
        !e.g. Bethke, eq.7.4, p110
        vLnGam(iAq)= vZ2(iAq) *LnGam_Charged
      ELSE
        !--- Neutral Species --!
        vLnGam(iAq)= LnGam_Neutral
      ENDIF
    ENDDO
    !--------------------------------------------------------/ Solute --
    !
    !-------------------------------------------------------- Solvent --
    Sum_Solute = SUM(vMolal(:))
    !-- LewisRandall,23-39,p347; EQ3NR-Doc, p.39
    LnActSv= Ln10 *Two/3.0D0 *dhA  *SqrtIost**3 *Sigma(SqrtIoSt) &
    &      - Ln10            *dhA  *IonStr*IonStr &
    &      - Sum_Solute
    !--- provisoire
    !LnActSv= Ln10 *Two/3.0D0 *dhA  *SqrtIost**3 *Sigma(SqrtIoSt) !&
    !&      - Sum_Solute
    X= Ln10*dhA /IonStr/3.375D0 *Sigma(1.5D0*SqrtIoSt) - One
    !LnActSv= Zero
    !---/
    OsmoSv= - X
    LnActSv=  X *MWSv *Sum_Solute
    !-------------------------------------------------------/ Solvent --!

  END SUBROUTINE Solmodel_Calc_SIT

  REAL(dp) FUNCTION Sigma(x)
  !-- for computing Water Acivity,
  !-- -> equ.23-40 in LewisRandall'61, p.347, used also in EQ3NR
    REAL(dp),INTENT(IN):: x
    
    REAL(dp):: XX
    
    XX= One+x
    Sigma=  3.0D0*(XX - One/XX - Two*LOG(XX)) /X**3
    
  END FUNCTION Sigma
  
ENDMODULE M_Solmodel_Calc_SIT
