!! *******************************************************************
!! * File name:    M_common_KINXIM
!! *
!! * Purpose:      ACCES AUX VARIABLES UTILES POUR M_EXCHANGE_KINXIM
!! *
!! * Author:       Anthony Michel (anthony.michel@ifp.fr)
!! *
!! * Created:      2006
!! * Updated:      2007  J.Moutte, A.Michel
!! *
!! ********************************************************************

module M_common_KINXIM_BOX

  use M_Global_Vars,   only : vSpc, vEle, vKinFas
  use M_System_Vars,   only : vCpn

  use M_Box_Param_Vars,  only : nAq, nMk, nCp
  
  use M_Box_System_Vars, only : &
       & tAlfAq, tAlfMk,        &
       & vOrdMk, vOrdAq
  
  use M_Box_Vars, only :   & 
       & Time, DTime,      &
       & VBox, TdgKBox, PbarBox, &
       & vMolK, vMolF, &
       & vQinj, Vout      
  
  use M_Box_Newton_Vars,   only : NewtonIter, Total_NewtonIter
  use M_Box_TimeLoop_Vars, only : Tfinal, iTimeStep, TimeLoop_Succes 
  use M_Box_Thermo_Vars,   only : vLnGam, vRhoK, vOmegaK, Ionic
  use M_Box_Kinetics_Vars, only : vStatusK
  use M_Box_Coupler_Vars,  only : LCouplerInject, LCouplerOutflow, IerrorChemistry

  implicit none
  private
  
  !// Species Global
  public :: vSpc
  public :: vEle
  public :: vCpn
  public :: vKinFas

   !// Sizes 
  public :: nAq
  public :: nMk
  public :: nCp
  
  !// Conditions
  public :: Vbox
  public :: TdgKBox
  public :: PbarBox
 
  !// State
  public :: VmolF
  public :: VmolK
 
  !// Species Index
  public :: vOrdAq
  public :: vOrdMk

  !// Species Formula
  public :: tAlfAq
  public :: tAlfMk
  
  !// Time
  public :: Tfinal
  public :: Time
  public :: DTime
  
  !// Numerics
  public :: NewtonIter, Total_NewtonIter
  public :: iTimeStep, TimeLoop_Succes 
  
  !// Thermo Vars
  public :: vLnGam
  public :: vRhoK
  public :: vOmegaK
  public :: Ionic
   
  !// Kinetic Vars
  public :: vStatusK

  !// Source Term
  public :: vQinj

  !// Flow Out
  public :: Vout
  
  public :: LCouplerInject
  public :: LCouplerOutflow
  public :: IerrorChemistry


end module M_common_KINXIM_BOX
