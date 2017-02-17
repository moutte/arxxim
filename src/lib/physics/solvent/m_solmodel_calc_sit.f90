module M_Solmodel_Calc_SIT
!-----------------------------------------------------------------------
! SIT model of activity for Aqueous Species
! after Grenthe et al, Modelling In Aquatic Chemistry,OECD,Chap.IX,p.330
!-----------------------------------------------------------------------

  use M_Numeric_Const
  use M_Kinds
  use M_Trace,only: fTrc,iDebug,T_
  
  implicit none
  private

  public:: Solmodel_Calc_SIT

contains

  subroutine Solmodel_Calc_SIT( &
  & vZSp, vMolal, &
  & dhA, Model, MWSv, &
  & vLnGam, LnActSv, OsmoSv)
  
    integer, intent(in)     :: vZSp(:)
    real(dp),intent(in)     :: vMolal(:)
    real(dp),intent(in)     :: dhA
    real(dp),intent(in)     :: MWSv
    character(*),intent(in) :: Model
    real(dp),intent(out)    :: vLnGam(:)
    real(dp),intent(out)    :: LnActSv
    real(dp),intent(out)    :: OsmoSv
    !
    real(dp):: LnGam_Charged, LnGam_Neutral
    real(dp):: Sum_Solute
    real(dp):: IonStr, SqrtIoSt, X
    integer :: vZ2(size(vZSp))
    integer :: nSolut,iAq
    
    nSolut= size(vMolal)
    vZ2(:)= vZSp(:)*vZSp(:)
    !    
    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= dot_product(vZ2(:),vMolal(:)) /Two
    !
    SqrtIoSt= SQRT(IonStr)
    !
    !------------------------------------------------------------ Solute
    LnGam_Charged = -dhA *(SqrtIoSt/(One + 1.5D0*SqrtIoSt))*Ln10
    !
    LnGam_Neutral = Zero
    !
    do iAq=1,nSolut
      if (vZSp(iAq)/=0) then
        !--- Charged Species --!
        !e.g. Bethke, eq.7.4, p110
        vLnGam(iAq)= vZ2(iAq) *LnGam_Charged
      else
        !--- Neutral Species --!
        vLnGam(iAq)= LnGam_Neutral
      end if
    end do
    !-----------------------------------------------------------/ Solute
    !
    !----------------------------------------------------------- Solvent
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
    !----------------------------------------------------------/ Solvent

  end subroutine Solmodel_Calc_SIT

  real(dp) function Sigma(x)
  !-- for computing Water Acivity,
  !-- -> equ.23-40 in LewisRandall'61, p.347, used also in EQ3NR
    real(dp),intent(in):: x
    
    real(dp):: XX
    
    XX= One+x
    Sigma=  3.0D0*(XX - One/XX - Two*log(XX)) /X**3
    
  end function Sigma
  
end module M_Solmodel_Calc_SIT
