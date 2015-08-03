MODULE M_Solmodel_Calc_Dilute
  !---------------------------------------------------------------
  ! models of activity for dilute aqueous solutions
  !
  ! IDEAL : a_i = mi , gamme_i = 1
  !         a_w = xw , gamma_w = 1
  !
  !---------------------------------------------------------------

  USE M_Numeric_Const
  USE M_Kinds
  USE M_Trace,    ONLY: fTrc,iDebug,T_
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC:: Solmodel_Calc_Dilute

CONTAINS

  SUBROUTINE Solmodel_Calc_Dilute( & 
  & vLnGam, LnActSv, OsmoSv) !,ERR)
    
    !REAL(dp),INTENT(IN)     :: vMolal(:) != solute species molalities
    !REAL(dp),INTENT(IN)     :: MolWeitSv
    REAL(dp),INTENT(OUT)    :: vLnGam(:)
    REAL(dp),INTENT(OUT)    :: LnActSv
    REAL(dp),INTENT(OUT)    :: OsmoSv
    
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

  END SUBROUTINE Solmodel_Calc_Dilute
  
ENDMODULE M_Solmodel_Calc_Dilute
