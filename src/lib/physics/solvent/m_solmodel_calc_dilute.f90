module M_Solmodel_Calc_Dilute
  !---------------------------------------------------------------
  ! models of activity for dilute aqueous solutions
  !
  ! IDEAL : a_i = mi , gamme_i = 1
  !         a_w = xw , gamma_w = 1
  !
  !---------------------------------------------------------------

  use M_Numeric_Const
  use M_Kinds
  use M_Trace,    only: fTrc,iDebug,T_
  
  implicit none
  
  private

  public:: Solmodel_Calc_Dilute

contains

  subroutine Solmodel_Calc_Dilute( & 
  & vLnGam, LnActSv, OsmoSv) !,ERR)
    
    !real(dp),intent(in)     :: vMolal(:) != solute species molalities
    !real(dp),intent(in)     :: MolWeitSv
    real(dp),intent(out)    :: vLnGam(:)
    real(dp),intent(out)    :: LnActSv
    real(dp),intent(out)    :: OsmoSv
    
    !//================ Solute Species
    vLnGam(:)= Zero

    !//================ Solvent
    !solvent activity= mole fraction = Nw/(Nw+SumNs) = 1 /(1+SUM(Ns/Nw) = 1/(1+ MW.SUM(m_i)
    !! LnActSv= -log( One + MolWeitSv *sum(vMolal(:)) )
    !use approximation for dilute solutions: Solvent>>Solute
    LnActSv= Zero
    !
    OsmoSv=  One
    !that's an approximation, valid for a DILUTE ideal solution
    !lnActW= ln(1/(1+ MW*SUM(m_i)) ~= - MW*SUM(m_i)
    !-> OsmoSv = -lnActW /MW /SUM(m_i) ~= +1

  end subroutine Solmodel_Calc_Dilute
  
end module M_Solmodel_Calc_Dilute
