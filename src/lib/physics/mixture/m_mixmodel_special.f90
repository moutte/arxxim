module M_MixModel_Special
  use M_Kinds
  use M_Dtb_Const, only: R_jk
  use M_Trace,only: Stop_
  implicit none
  !
  private
  !
  public:: MixModel_Special_Activities

contains

subroutine MixModel_Special_Activities( &
& sKey,       &
& TdgK, Pbar, &
& vX,         &
& vLPole,     &
& vLActIdl,   &
& vLGam)
  character(len=*),intent(in) :: sKey
  real(dp),intent(in) :: TdgK, Pbar
  real(dp),intent(in) :: vX(:)
  logical, intent(in) :: vLPole(:)
  real(dp),intent(out):: vLActIdl(:)
  real(dp),intent(out):: vLGam(:)
  !
  integer:: N
  real(dp),allocatable:: vActId(:)
  !
  N= size(vX)
  allocate(vActId(N))
  !
  select case(trim(sKey))

  case("SPECIAL_2")        ;  call Special_2(TdgK,vX,vActId,vLGam)
  case("SPECIAL_3")        ;  call Special_3(TdgK,vX,vActId,vLGam)

  case("FLUID_KERJAC")     ;  call Fluid_KerJac(TdgK,Pbar,vX,vActId,vLGam)
  case("FLUID_DUAN92")     ;  call Fluid_Duan92(TdgK,Pbar,vX,vActId,vLGam)

  case("FELSPAR_GH84")     ;  call Felspar_Gh84(TdgK,vX,vActId,vLGam)
  case("FELSPAR_HP03A")    ;  call Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)
  case("FELSPAR_HP03B")    ;  call Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)
  case("FELSPAR_HP03C")    ;  call Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)

  case("CORDIERITE_HP05")  ;  call Cordierite_HP05(TdgK,vX,vActId,vLGam)
  case("CORDIERITE_HPWEB") ;  call Cordierite_HPWeb(TdgK,vX,vActId,vLGam)

  case("BIOTITE_HP05")     ;  call Biotite_HP05(TdgK,vX,vActId,vLGam)
  case("BIOTITE_HP98")     ;  call Biotite_HP98(TdgK,vX,vActId,vLGam)
  case("BIOTITE_HPWEB")    ;  call Biotite_HPWeb(TdgK,vX,vActId,vLGam)
  case("MUSCOVITE_HP05")   ;  call Muscovite_HPWeb(TdgK,Pbar,vX,vActId,vLGam)
  case("MICA_VIDAL")       ;  call Mica_Vidal(TdgK,Pbar,vX,vActId,vLGam)

  case("OPX_HP05")         ;  call Opx_HP05(TdgK,vX,vActId,vLGam)
  case("CHLORITE_VIDAL")   ;  call Chlorite_Vidal(TdgK,Pbar,vX,vActId,vLGam)
  case("CHLORITE_VIDAL2")  ;  call Chlorite_Vidal_Ideal(TdgK,Pbar,vX,vActId,vLGam)
  case("CHLORITE_FE")      ;  call Chlorite_Fe(TdgK,vX,vActId,vLGam)
  case("CHLORITE_HP05")    ;  call Chlorite_HPWeb(TdgK,vX,vActId,vLGam)

  !--- AVCHENKO MODELS
  case("BIOTITE_AVCH")     ;  call Biotite_Avchenko(TdgK,vX,vActId,vLGam)
  case("CORDIERITE_AVCH")  ;  call Cordierite_Avchenko(TdgK,vX,vActId,vLGam)
  case("OPX_AVCH")         ;  call Opx_Avchenko(TdgK,vX,vActId,vLGam)
  case("TALC_AVCH")        ;  call Talc_Avchenko(TdgK,vX,vActId,vLGam)
  case("EPIDOTE_AVCH")     ;  call Epidote_Avchenko(TdgK,vX,vActId,vLGam)
  case("ALKSPAR_AVCH")     ;  call Alkspar_Avchenko(TdgK,Pbar,vX,vActId,vLGam)
  case("PLAGIO_AVCH")      ;  call Plagio_Avchenko(TdgK,vX,vActId,vLGam)
  !---/AVCHENKO MODELS

  case default
    call Stop_(trim(sKey)//" is NOT A SPECIAL MODEL")

  end select
  !
  where(vLPole(:))  ;  vLActIdl(:)= log(vActId(:))
  elsewhere         ;  vLActIdl(:)= Zero
  end where
  !
  deallocate(vActId)
  !
  return
end subroutine MixModel_Special_Activities

subroutine Special_2( &
& TdgK,               &
& X,                  &
& vActId,vLGam)
  !MIXTURE.MODEL SPECIAL_01
  !  MODEL SPECIAL
  !  POLE
  !    AA_M1
  !    AB_M1
  !  end
  !end
  !---------------------------------------------------------------------
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !---------------------------------------------------------------------
  real(dp):: W12
  !---------------------------------------------------------------------
  W12= 3.D4
  !
  vActId(1)= X(1)
  vActId(2)= X(2)
  !
  vLGam(1)= X(2)*(One-X(1))*W12
  vLGam(2)= X(1)*(One-X(2))*W12
  vLGam(:)= vLGam(:) /R_jk/TdgK
  !
end subroutine Special_2

subroutine Special_3( &
& TdgK, &
& X,    &
& vActId,vLGam)
  !MIXTURE.MODEL SPECIAL_03
  !  MODEL SPECIAL
  !  POLE
  !    MIN1
  !    MIN2
  !    MIN3
  !  end
  !end
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  real(dp):: W12,W23,W31
  !---------------------------------------------------------------------
  W12= 3.D4  ;  W23= 3.D4  ;  W31= 3.D4
  !
  vActId(1)= X(1)
  vActId(2)= X(2)
  vActId(3)= X(3)
  !
  !--------------------------------------------- regular mixing model --
  vLGam(1)= (One-X(1)) *( X(2)*W12 + X(3)*W31 ) - X(2)*X(3)*W23
  vLGam(2)= (One-X(2)) *( X(3)*W23 + X(1)*W12 ) - X(3)*X(1)*W31
  vLGam(3)= (One-X(3)) *( X(1)*W31 + X(2)*W23 ) - X(1)*X(2)*W12
  vLGam(:)= vLGam(:) /R_jk/TdgK
  !
end subroutine Special_3

subroutine Fluid_KerJac(TdgK,Pbar,vX,vActId,vLGam)
  use M_Eos_KerJac
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  real(dp):: XCO2
  real(dp):: LnPhiH2O_Pur,LnPhiCO2_Pur
  real(dp):: LnPhiH2O_Mix,LnPhiCO2_Mix
  real(dp):: LnActH2O,    LnActCO2
  real(dp):: V_H2O_Pur,   V_CO2_Pur
  real(dp):: V_m3
  logical :: Update_Pure
  !---------------------------------------------------------------------
  vActId(:)= vX(:)
  vLGam(:)=  Zero !vLGam(:) /R_jk/TdgK
  !
  XCO2= vX(2)
  !
  Update_Pure= .true.
  !
  call EoS_H2O_CO2_KerJac( & !
  & TdgK,Pbar,XCO2,              & !in
  & Update_Pure,                 & !in
  & LnPhiH2O_Pur,LnPhiCO2_Pur,   & !inout
  & V_H2O_Pur,V_CO2_Pur,         & !inout !volumes in cm3 !!!
  & LnPhiH2O_Mix,LnPhiCO2_Mix,   & !out
  & LnActH2O,    LnActCO2,       & !out
  & V_m3)
  !
  vLGam(1)= LnActH2O - log(vX(1))
  vLGam(2)= LnActCO2 - log(vX(2))
  !
  return
end subroutine Fluid_KerJac

subroutine Fluid_Duan92(TdgK,Pbar,vX,vActId,vLGam)
  use M_Mixmodel_Duan
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  real(dp):: LnPhiH2O_Pur,LnPhiCO2_Pur
  real(dp):: V_H2O_Pur,   V_CO2_Pur
  real(dp):: LnActH2O,    LnActCO2
  real(dp):: V_m3
  logical :: Update_Pure
  !---------------------------------------------------------------------
  vActId(:)= vX(:)
  vLGam(:)=  Zero !vLGam(:) /R_jk/TdgK
  !
  Update_Pure= .true.
  LnPhiH2O_Pur= 0.D0
  LnPhiCO2_Pur= 0.D0
  !
  call EoS_H2O_CO2_Duan92(     & !
  & TdgK,Pbar,vX,              & !in
  & Update_Pure,               & !in
  & LnPhiH2O_Pur,LnPhiCO2_Pur, & !inout
  & V_H2O_Pur,   V_CO2_Pur,    & ! inout
  & V_m3,LnActH2O,LnActCO2)      ! OUT
  !
  vLGam(1)= LnActH2O -log(vX(1))
  vLGam(2)= LnActCO2 -log(vX(2))
  !
  return
end subroutine Fluid_Duan92

subroutine Plagio_Wrk( &
& TdgK, &
& X,    &
& vActId,vLGam)
!--
!-- G-FSP- Ghiorso, 1984
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL FELSPAR_WORK
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE
  !    ANORTHITE
  !  end
  !end

  vActId(1)= X(1) *( One -X(2)**2)       ! albite
  vActId(2)= X(2) *((One +X(2))/2.0D0)**2 ! anorthite

  !-- G-FSP- Ghiorso, 1984
  vLGam(1)= X(2)**2 *( 28211.0974D0 -39485.4998D0*X(1) )
  vLGam(2)= X(1)**2 *( 8468.3475D0  +39485.4998D0*X(2) )
  vLGam(3)= Zero
  !
  vLGam(:)= vLGam(:) /R_jk/TdgK

  return
end subroutine Plagio_Wrk

subroutine Felspar_Gh84( &
& TdgK, &
& X,    &
& vActId,vLGam)
!--
!-- G-FSP- Ghiorso, 1984
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL FELSPAR_GH84
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE
  !    ANORTHITE
  !    SANIDINE
  !  end
  !end
  !
  vActId(1)= X(1) *(Two -X(1))
  vActId(2)= X(2) *(One +X(2)) /4.0D0
  vActId(3)= X(3)
  !
  !-- G-FSP- Ghiorso, 1984
  vLGam(1)= X(2)**2 *( 28211.0974D0 -39485.4998D0*X(1) )
  vLGam(2)= X(1)**2 *( 8468.3475D0  +39485.4998D0*X(2) )
  vLGam(3)= Zero
  !
  vLGam(:)= vLGam(:) /R_jk/TdgK
  !
  return
end subroutine Felspar_Gh84

subroutine Felspar_HP03( &
& TdgK, Pbar, &
& X, &
& vActId,vLG)
!--
!-- (Holland and Powell (2003)) -- David Dolejs, 13-March-05
!--
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL FELSPAR_HP03A
  !  MODEL SPECIAL
  !  POLE
  !    MICROCLINE
  !    ALBITE
  !    ANORTHITE-C1
  !  end
  !end
  !! ..or..
  !MIXTURE.MODEL FELSPAR_HP03B
  !  MODEL SPECIAL
  !  POLE
  !    SANIDINE
  !    ALBITE
  !    ANORTHITE-C1
  !  end
  !end
  !! ..or..
  !MIXTURE.MODEL FELSPAR_HP03C
  !  MODEL SPECIAL
  !  POLE
  !    SANIDINE
  !    HIGH-ALBITE
  !    ANORTHITE-C1
  !  end
  !end
  !
  real(dp):: V1,V2,V3,F1,F2,F3,P
  !
  P=  Pbar *1.0D-3 ! P in kbar !!
  !
  V1= One
  V2= 0.643D0
  V3= One
  !
  F1= X(1) *V1 /( X(1)*V1 +X(2)*V2 +X(3)*V3 )
  F2= X(2) *V2 /( X(1)*V1 +X(2)*V2 +X(3)*V3 )
  F3= X(3) *V3 /( X(1)*V1 +X(2)*V2 +X(3)*V3 )
  !
  vLG(1)= (One-F1)*F2 *Two *V1/(V1+V2) *(25.1D3 -10.8D0*TdgK+0.343D0*P) &
  &     + (One-F1)*F3 *Two *V1/(V1+V3) *40.0D3 &
  &     -       F2*F3 *Two *V1/(V2+V3) *3.1D3
  !
  vLG(2)= (One-F2)*F1 *Two *V2/(V2+V1) *(25.1D3 -10.8D0*TdgK+0.343D0*P) &
  &     + (One-F2)*F3 *Two *V2/(V2+V3) *3.1D3 &
  &     -       F1*F3 *Two *V2/(V1+V3) *40.0D3
  !
  vLG(3)= (One-F3)*F1 *Two *V3/(V1+V3) *40.0D3 &
  &     + (One-F3)*F2 *Two *V3/(V2+V3) *3.1D3  &
  &     -       F1*F2 *Two *V3/(V2+V1) *(25.1D3 -10.8D0*TdgK+0.343D0*P)
  !
  !LG(1)=X(1)*Dexp(A1/(R*T))
  !LG(2)=X(2)*Dexp(A2/(R*T))
  !LG(3)=X(3)*Dexp(A3/(R*T))
  !
  vActId(:)= X(:)
  vLG(:)= vLG(:)/R_jk/TdgK
  !
end subroutine Felspar_HP03

subroutine Plagio_Avchenko( &
& TdgK, &
& vX,   &
& vActId,vLG)
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: Xb,WI,WC,Iab,Ian
  
  !MIXTURE.MODEL PLAGIO_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE-LOW
  !    ANORTHITE
  !  end
  !end
  
  vActId(:)= vX(:)
  
  WC= 1.07D3
  WI= 9.79D3
  Xb= 0.12D0 + TdgK *0.00038D0
  Ian= (WI-WC) *(One-Xb)**2
  Iab= (WC-WI) *Xb**2
  
  if(vX(2) > Xb) then
    vLG(1)= WI *vX(2)**2 +Iab
    vLG(2)= WI *vX(1)**2
  else
    vLG(1)= WC *vX(2)**2
    vLG(2)= WC *vX(1)**2 +Ian
  end if
  
  vLG(:)= vLG(:)/R_jk/TdgK
  
  return
end subroutine Plagio_Avchenko

subroutine Alkspar_Avchenko( &
& TdgK,Pbar, &
& vX,   &
& vActId,vLG)
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: W_Na,W_K
  !
  !MIXTURE.MODEL ALKSPAR_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE-LOW
  !    MICROCLINE
  !  end
  !end
  !
  vActId(:)= vX(:)
  !
  W_Na= 7.594D3 -TdgK*5.931D0 +Pbar*0.142D0
  W_K=  7.832D3 -TdgK*2.657D0 +Pbar*0.074D0
  !
  vLG(1)= vX(1)**2 *( vX(2)*Two*W_Na + (One-vX(2)*Two)*W_K  )
  vLG(2)= vX(2)**2 *( vX(1)*Two*W_K  + (One-vX(1)*Two)*W_Na )
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Alkspar_Avchenko

subroutine Cordierite_HPWeb( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/thermocalc-cordierites
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL CORDIERITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !  end
  !end
  !
  vActId(1)= ( vX(1)+vX(3) )**2 *( vX(1)+vX(2) ) ! X_Mg**2 *(1-X_H2O)
  vActId(2)=   vX(2)        **2 *( vX(1)+vX(2) ) ! X_Fe**2 *(1-X_H2O)
  vActId(3)= ( vX(1)+vX(3) )**2 *  vX(3)         ! X_Mg**2 *X_H2O
  !
  vLGam(:)= Zero
  !
end subroutine Cordierite_HPWeb

subroutine Cordierite_HP05( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- CORDIERITE, IDEAL RECIPROCAL -- David Dolejs, 13-March-05
!-- code derived from THERIAK/fsol.f90
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL CORDIERITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !    H.FE-CORDIERITE
  !  end
  !end
  !
  vActId(1)= ( vX(1)+vX(3) )**2 *( vX(1)+vX(2) )
  vActId(2)= ( vX(2)+vX(4) )**2 *( vX(1)+vX(2) )
  vActId(3)= ( vX(1)+vX(3) )**2 *( vX(3)+vX(4) )
  vActId(4)= ( vX(2)+vX(4) )**2 *( vX(3)+vX(4) )
  !
  vLGam= Zero
  !
end subroutine Cordierite_HP05

subroutine Cordierite_Avchenko( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!--
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  real(dp):: W_CrdFcrd
  !
  !MIXTURE.MODEL CORDIERITE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !  end
  !end
  !
  !-- crd    Mg2  Al4Si5O18         (crd + hcrd)**2 *(1-hcrd)
  !-- fcrd   Fe2  Al4Si5O18          fcrd**2        *(1-hcrd)
  !-- hcrd   Mg2  Al4Si5O18  H2O    (crd + hcrd)**2 *hcrd
  !
  vActId(1)= (vX(1) + vX(3))**2  *(One-vX(3))
  vActId(2)=  vX(2)**2           *(One-vX(3))
  vActId(3)= (vX(1) + vX(3))**2  *vX(3)
  !
  vLGam(:)= Zero
  !
  W_CrdFcrd= 1.5D3
  vLGam(1)= (One-vX(1)) *vX(2) *W_CrdFcrd
  vLGam(2)= (One-vX(2)) *vX(1) *W_CrdFcrd
  vLGam(3)= -vX(2)      *vX(1) *W_CrdFcrd
  !
  vLGam(:)= vLGam(:) /TdgK/R_jk
  !
end subroutine Cordierite_Avchenko

subroutine BiotiteMgAl_HPWeb( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/biotite
!--
!-- KMASH:
!-- The ideal mixing model used in HP90 involved mixing of
!-- octahedral Al and Mg on the M2 and M3 sites in the trioctahedral
!-- mica structure, M1 being occupied only by Mg.
!-- This may be unrealistic given that Mg-Al ordering is such a strong
!-- control on the behaviour of octahedral site distributions
!-- (e.g. in omphacite and chlorite).
!-- This short range order situation is not easy to model satisfactorily,
!-- but energetically this would be much more like a hypothetical fully
!-- ordered mica in which Al and Mg can mix on only one (M3)
!-- of the three octahedral sites, such that each Al is surrounded
!-- by six Mg nearest neighbours in the octahedral sheet.
!--
!-- Coding in THERMOCALC ----------------------
!  %========= Biotites in KMASH ===============
!  % -- phlogopite 3M-2T model ---
!  bi 2
!  y(bi) 0.01
!  %
!  p(phl)  1   1    1  1 -1  y
!  p(east) 1   1    0  1  1  y
!  %
!  sf
!  W(phl,east)  10  0  0
!  %
!  4     % no of site fractions
!  x(Al,M1)  1 1    0  1  1   y
!  x(Mg,M1)  1 1    1  1 -1   y
!  x(Al,T1)  1 1    1/2 1 1/2  y
!  x(Si,T1)  1 1    1/2 1 -1/2 y
!  %
!  phl      4  3    x(Mg,M1) 1    x(Al,T1) 1   x(Si,T1) 1
!  east     1  2    x(Al,M1) 1    x(Al,T1) 2
!  %
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    EASTONITE
  !  end
  !end
  !
  !--- ideal site mixing activities
  !-- a_Phl- K_A *Mg_M3 *Al_T1 *Si_T1 *4
  !-- a_Eas- K_A *Al_M3 *Al_T1 *Al_T1
  !
  vActId(1)= (One -vX(2))**2 *(One +vX(2))
  vActId(2)=  vX(2)          *(One +vX(2))**2 /4.0D0
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= vX(2) *(One-vX(1)) *10.0D3  ! W-Phl-Eas
  vLG(3)= vX(1) *(One-vX(2)) *10.0D3  ! W-Phl-Eas
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine BiotiteMgAl_HPWeb

subroutine Biotite_HPWeb( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/biotite
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !  end
  !end
  !
  real(dp):: Phl,Ann,Eas,Obi
  real(dp):: W_PhlAnn,W_PhlObi,W_PhlEas,W_AnnEas,W_EasObi,W_AnnObi
  !
  W_PhlAnn=  9.0D3
  W_PhlObi=  3.0D3
  W_PhlEas=  10.0D3
  W_AnnEas= -1.0D3
  W_EasObi=  10.0D3
  W_AnnObi=  6.0D3
  !
  Phl= X(1)
  Ann= X(2)
  Eas= X(3)
  Obi= X(4)
  !
  !--- ideal site mixing activities
  vActId(1)= Phl            *(One -Ann)**2 *(One -Eas) *(One +Eas)
  vActId(2)= (One -Eas-Phl) *Ann**2        *(One -Eas) *(One +Eas)
  vActId(3)= Eas /4.0D0     *(One -Ann)**2 *(One +Eas)**2
  vActId(4)= (One -Eas-Phl) *(One-Ann)**2  *(One -Eas) *(One +Eas)
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= Ann *(One-Phl) *W_PhlAnn &
  &     + Obi *(One-Phl) *W_PhlObi &
  &     + Eas *(One-Phl) *W_PhlEas &
  &     - Eas * Ann      *W_AnnEas &
  &     - Obi * Eas      *W_EasObi &
  &     - Obi * Ann      *W_AnnObi
  !
  vLG(2)= Phl *(One-Ann) *W_PhlAnn  &
  &     + Obi *(One-Ann) *W_AnnObi  &
  &     + Eas *(One-Ann) *W_AnnEas  &
  &     - Obi *Phl       *W_PhlObi  &
  &     - Obi *Eas       *W_EasObi  &
  &     - Eas *Phl       *W_PhlEas
  !
  vLG(3)= Phl *(One-Eas) *W_PhlEas  &
  &     + Ann *(One-Eas) *W_AnnEas  &
  &     + Obi *(One-Eas) *W_EasObi  &
  &     - Phl *Ann       *W_PhlAnn  &
  &     - Obi *Ann       *W_AnnObi  &
  &     - Obi *Phl       *W_PhlObi
  !
  vLG(4)= Ann *(One-Obi) *W_AnnObi  &
  &     + Phl *(One-Obi) *W_PhlObi  &
  &     + Eas *(One-Obi) *W_EasObi  &
  &     - Phl *Ann       *W_PhlAnn  &
  &     - Eas *Phl       *W_PhlEas  &
  &     - Ann *Eas       *W_AnnEas
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Biotite_HPWeb

subroutine Biotite_Avchenko( &
& TdgK, &
& vX,    &
& vActId,vLG)
!--
!--
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !  end
  !end
  !
  real(dp):: Phl,Ann,Eas
  real(dp):: W_PhlAnn,W_PhlEas,W_AnnEas
  !
  W_PhlAnn=  9.0D3
  W_PhlEas=  10.0D3
  W_AnnEas= -1.0D3
  !
  Phl= vX(1)
  Ann= vX(2)
  Eas= vX(3)
  !
  !--                 M1   M2   T1
  !--  phl         K  Mg   Mg2  AlSi  Si2O10(OH)2
  !--  ann         K  Fe   Fe2  AlSi  Si2O10(OH)2
  !--  eas         K  Al   Mg2  Al2   Si2O10(OH)2
  !--
  !--- ideal site mixing activities
  vActId(1)= Phl   *(One -Ann)**2 *(One -Eas) *(One +Eas)
  vActId(2)= Ann   *Ann**2        *(One -Eas) *(One +Eas)
  vActId(3)= Eas   *(One -Ann)**2 *(One +Eas)**2    /4.0D0
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= Ann *(One-Phl) *W_PhlAnn &
  &     + Eas *(One-Phl) *W_PhlEas &
  &     - Eas * Ann      *W_AnnEas
  !
  vLG(2)= Phl *(One-Ann) *W_PhlAnn  &
  &     + Eas *(One-Ann) *W_AnnEas  &
  &     - Eas *Phl       *W_PhlEas
  !
  vLG(3)= Phl *(One-Eas) *W_PhlEas  &
  &     + Ann *(One-Eas) *W_AnnEas  &
  &     - Phl *Ann       *W_PhlAnn
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Biotite_Avchenko

subroutine Muscovite_HPWeb( &
& TdgK,Pbar, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/muscovite
!--
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: X,Y,Z
  real(dp):: Mus,Par,Cel,Fcl
  real(dp):: W_MusPar, W_CelPar, W_ParFcl, I_Par
  !
  !MIXTURE.MODEL MUSCOVITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    MUSCOVITE
  !    PARAGONITE
  !    CELAdoNITE
  !    FE-CELAdoNITE
  !  end
  !end
  !
  Mus= vX(1)
  Par= vX(2)
  Cel= vX(3)
  Fcl= vX(4)
  !
  X= Fcl /(Cel + Fcl) != bulk Fe/(Fe+Mg)
  Y= Mus +Par          != X_Al_M2
  Z= Par                != X_Na_A
  !
  ! Par= Z
  ! Cel= (1-X)(1-Y)
  ! Fcl= X(1-Y)
  ! Mus= Y-Z
  !
  !--- ideal site mixing activities
  vActId(1)= (One-Z) *Y**2             *(Two-Y)            != Mus
  vActId(2)=      Z  *Y**2             *(Two-Y)            != Par
  vActId(3)= (One-Z) *(One-Y) *(One-X) *(Two-Y)**2 /4.0D0  != Cel
  vActId(4)= (One-Z) *(One-Y) *X       *(Two-Y)**2 /4.0D0  != Fcl
  !
  W_MusPar= 12.0D3 +0.4D0*Pbar  ! 0.4D3 *Par/1.D3
  W_CelPar= 14.0D3 +0.2D0*Pbar
  W_ParFcl= 14.0D3 +0.2D0*Pbar
  I_Par=    1.42D3 +0.4D0*Pbar
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= (One -Mus) *Par *W_MusPar &
  &     - Cel        *Par *W_CelPar &
  &     - Fcl        *Par *W_ParFcl
  vLG(2)= (One -Par) *(Mus*W_MusPar  + Cel*W_CelPar + Fcl*W_ParFcl) &
  &     + I_Par
  vLG(3)= (One -Cel) *Par *W_CelPar &
  &     - Par        *Mus *W_MusPar &
  &     - Par        *Fcl *W_ParFcl
  vLG(4)= (One-Fcl)  *Par *W_ParFcl &
  &     - Par        *Mus *W_MusPar &
  &     - Par        *Cel *W_CelPar
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Muscovite_HPWeb

subroutine Talc_HPWeb( & !
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/talc
!--                        M3  M2   T1
!-- 1- ta   talc           Mg  Mg2  Si2  Si2 O10(OH)2
!-- 2- fta  ferrotalc      Fe  Fe2  Si2  Si2 O10(OH)2
!-- 3- tats tschermak-talc Al  Mg2  AlSi Si2 O10(OH)2
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL TALC_HP
  !  MODEL SPECIAL
  !  POLE
  !    TALC
  !    FE-TALC
  !    TSCHERMAK-TALC
  !  end
  !end
  !
  real(dp):: X,Y
  !
  X= vX(2) /(One -vX(3)/3.0D0) !X =Fe/(Fe+Mg) =3.P2/(3.P1+2.P3+3.P2)
  Y= vX(3) != X_Al_M2
  !
  !--- ideal site mixing activities
  vActId(1)= (One-X)**3 *(One-Y)   *(Two-Y)**2 /4.0D0 !XMg_M3*XSi_T1**2
  vActId(2)=      X**3  *(One-Y)   *(Two-Y)**2 /4.0D0 !XFe_M3*XSi_T1**2
  vActId(3)= (One-X)**2 *     Y**2 *(Two-Y)           !XAl_M3*XAl_T1*XSi_T1
  !
  vLGam(:)= Zero
  !
  return
end subroutine Talc_HPWeb

subroutine Talc_Avchenko( & !
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- http://zhurnal.ape.relarn.ru/articles/2007/068.pdf
!--                        M3  M2   T1
!-- 1- ta   talc           Mg  Mg2  Si2  Si2 O10(OH)2
!-- 2- tats tschermak-talc Al  Mg2  AlSi Si2 O10(OH)2
!-- 3- fta  ferrotalc      Fe  Fe2  Si2  Si2 O10(OH)2
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  real(dp):: ta,fta,tats
  !
  !MIXTURE.MODEL TALC_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    TALC
  !    FE-TALC
  !    TSCHERMAK-TALC
  !  end
  !end
  !
  ta=   vX(1)
  fta=  vX(2)
  tats= vX(3)
  !
  !--- ideal site mixing activities
  vActId(1)= ta   *(ta+tats)**2 *(One -tats/Two)**2
  vActId(2)= fta**3             *(One -tats/Two)**2
  vActId(3)= tats *(ta+tats)**2 *(One -tats/Two) *(tats/Two) *4.0D0
  !
  vLGam(:)= Zero
  !
  return
end subroutine Talc_Avchenko

subroutine Epidote_Avchenko( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://zhurnal.ape.relarn.ru/articles/2007/068.pdf
!--               M3  M1
!-- 1- ep    Ca2  Fe  Al  AlSi3O12(OH)
!-- 2- cz    Ca2  Al  Al  AlSi3O12(OH)
!-- 3- fep   Ca2  Fe  Fe  AlSi3O12(OH)
!--
  !
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: Epi,Clz,Fep
  real(dp):: W_EpiFep,W_ClzFep
  !
  !MIXTURE.MODEL EPIdoTE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    EPIdoTE
  !    CLINOZOISITE
  !    FE-EPIdoTE
  !  end
  !end
  !
  Epi=  vX(1)
  Clz=  vX(2)
  Fep=  vX(3)
  !
  !--- ideal site mixing activities
  vActId(1)= (One -Clz) *(One -Fep)
  vActId(2)= Clz        *(One -Fep)
  vActId(3)= (One -Clz) *Fep
  !
  W_EpiFep= 3.0D3
  W_ClzFep= 15.4D3
  !
  vLG(1)= (One-Epi)*Fep *W_EpiFep -Clz*Fep *W_ClzFep
  vLG(2)= (One-Clz)*Fep *W_ClzFep -Epi*Fep *W_EpiFep
  vLG(3)= (One-Fep) *(Epi*W_EpiFep + Clz*W_ClzFep)
  !
  vLG(:)= vLG(:) /TdgK/R_jk
  !
  return
end subroutine Epidote_Avchenko

subroutine Chlorite_HPWeb( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/chlorite
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite 
  !    CLINOCHLORE  
  !    AMESITE      
  !    DAPHNITE     
  !  end
  !end
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_ 
  !endSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite 
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !endPOLE
  !
  !HP parameters X,Y,N=
  !X= X_Fe_M2  = Fe/(Fe+Mg)
  !Y= X_Al_T2_ = 
  !N=(X_Al_M4_ - X_Al_M1_)/2
  !
  !P1=         1 - X_Al_M4_              !X(CHLORITE-ALFREE)
  !P2+P4= 2N =     X_Al_M4_ - X_Al_M1_
  !P3=                        X_Al_M1_   !X(AMESITE)
  !
  !P1= 1   -Y -N
  !P2= 2*N -2*X*(3-Y)/5
  !P3= Y   -N
  !P4= 2*vX*(3-Y_)/5
  !
  real(dp):: X,Y,N
  real(dp):: P1,P2,P3,P4
  real(dp):: W12,W13,W14,W23,W24,W34
  !
  !! integer :: I,J,K
  !! real(dp):: tW(4,4),G
  !! !
  !! tW(:,:)= Zero
  !! tW(1,2)= 18.0D3  ! Afch-Clin
  !! tW(1,3)= 20.0D3  ! Afch-Ames
  !! tW(1,4)= 14.5D3  ! Afch-Daph
  !! tW(2,3)= 18.0D0  ! Clin-Ames
  !! tW(2,4)= 2.50D3  ! Clin-Daph
  !! tW(3,4)= 13.5D3  ! Ames-Daph
  !! !
  !! do I=1,4
  !!   G= Zero
  !!   do J=1,4
  !!     if(J/=I) G= G + vX(J)*tW(I,J)
  !!   end do
  !!   G= G *(One-vX(I))
  !!   do J=1,4
  !!     do K=1,4
  !!       if(J/=I .and. K/=I) G= G - vX(J)*vX(K)*tW(J,K)
  !!     end do
  !!   end do
  !!   vLG(I)= G
  !! end do
  !
  W12= 18.0D3  ! Afch-Clin
  W13= 20.0D3  ! Afch-Ames
  W14= 14.5D3  ! Afch-Daph
  W23= 18.0D0  ! Clin-Ames
  W24= 2.50D3  ! Clin-Daph
  W34= 13.5D3  ! Ames-Daph
  !
  P1= vX(1)  ! CHLORITE-AF  
  P2= vX(2)  ! CLINOCHLORE  
  P3= vX(3)  ! AMESITE      
  P4= vX(4)  ! DAPHNITE     
  !
  !! print *,"P1..P4=",P1,P2,P3,P4
  !! pause
  Y=(One +P3 -P1)/Two
  N=(One -P3 -P1)/Two
  !X= 5*P4 /(5*P4+6*P1+5*P2+4*P3)= 5*P4 /(5 + P1- P3)
  X= P4 / (One + (P1-P3)/5.0D0)
  !
  N= (P2+P4)/Two
  Y= N + P3
  X= P4
  !
  !! print *,"X,Y,Z=",X,Y,N
  !! pause
  !--- ideal activities (mixing on sites)
  vActId(1)= (One-X)**6 *(One-Y+N)*(One-Y-N)*(One-Y)*(One-Y)   ! CHLORITE-ALFREE
  vActId(2)= (One-X)**5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! CLINOCHLORE
  vActId(3)= (One-X)**4 *(    Y-N)*(    Y+N)*Y      *Y         ! AMESITE
  vActId(4)=      X **5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! DAPHNITE
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= (One-P1) *( P2*W12 + P3*W13 + P4*W14) &
  &     - P2 * P3 *W23     &
  &     - P2 * P4 *W24     &
  &     - P3 * P4 *W34      
  vLG(2)= (One-P2) *( P1*W12 + P3*W23 + P4*W24) &
  &     - P1 * P3 *W13     &
  &     - P1 * P4 *W14     &
  &     - P3 * P4 *W34      
  vLG(3)= (One-P3) *( P1*W13 + P2*W23 + P4*W34) &
  &     - P1 * P2 *W12     &
  &     - P1 * P4 *W14     &
  &     - P2 * P4 *W24      
  vLG(4)= (One-P4) *( P1*W14 + P2*W24 + P4*W34) &
  &     - P1 * P2 *W12     &
  &     - P1 * P3 *W13     &
  &     - P2 * P3 *W23      
  !
  vLG(:)= vLG(:) /R_jk/TdgK
  !
  return
end subroutine Chlorite_HPWeb

subroutine Chlorite_HPWeb_( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/chlorite
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !    DAPHNITE
  !  end
  !end
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !end
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !end
  !
  !HP parameters X,Y,N=
  !X= Fe_M2  = Fe/(Fe+Mg)
  !Y= Al_T2  =
  !N=(Al_M4 - Al_M1)/2
  !
  !P1=         1 - Al_M4           !X(CHLORITE-ALFREE)
  !P2+P4= 2N =     Al_M4 - Al_M1
  !P3=                     Al_M1   !X(AMESITE)
  !
  !P1= 1 -Y -N
  !P2= 2*N -2*X*(3-Y)/5
  !P3= Y -N
  !P4= 2*X*(3-Y)/5
  !
  ! integer :: i,j
  ! real(dp):: W(4,4)
  
  real(dp):: X,Y,N
  real(dp):: P1,P2,P3,P4
  real(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,MgM4,FeM4,AlM4,SiT2,AlT2
  real(dp):: W12,W13,W14,W23,W24,W34
  !
  W12= 18.0D3  ! Afch-Clin
  W13= 20.0D3  ! Afch-Ames
  W14= 14.5D3  ! Afch-Daph
  W23= 18.0D0  ! Clin-Ames
  W24= 2.50D3  ! Clin-Daph
  W34= 13.5D3  ! Ames-Daph
  !
  P1= vX(1)  ! CHLORITE-AF
  P2= vX(2)  ! CLINOCHLORE
  P3= vX(3)  ! AMESITE
  P4= vX(4)  ! DAPHNITE
  !
  N= (P2+P4)/Two
  Y= N + P3
  X= P4
  !
  !--- ideal activities (mixing on sites)
  vActId(1)= (One-X)**6 *(One-Y+N)*(One-Y-N)*(One-Y)*(One-Y)   ! CHLORITE-ALFREE
  vActId(2)= (One-X)**5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! CLINOCHLORE
  vActId(3)= (One-X)**4 *(    Y-N)*(    Y+N)*Y      *Y         ! AMESITE
  vActId(4)=      X **5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! DAPHNITE
  
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  
  FeM2= X      ;  MgM2= One -X
  AlT2= Y      ;  SiT2= One -Y
  AlM1= Y-N    ;  FeM1= (One-AlM1) *X  ;  MgM1= (One-AlM1) *(One-X)
  AlM4= Y+N    ;  FeM4= (One-AlM4) *X  ;  MgM4= (One-AlM4) *(One-X)
  
  vActId(1)= MgM2**4 *MgM1 *MgM4* SiT2 *SiT2         ! CHLORITE-AF
  vActId(2)= MgM2**4 *MgM1 *AlM4* AlT2 *SiT2 *4.0D0  ! CLINOCHL
  vActId(4)= MgM2**4 *AlM1 *AlM4* AlT2 *AlT2         ! AMESITE
  vActId(3)= FeM2**4 *FeM1 *AlM4* AlT2 *SiT2 *4.0D0  ! DAPHNITE
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= (One-P1) *( P2*W12 + P3*W13 + P4*W14) &
  &     - P2 * P3 *W23     &
  &     - P2 * P4 *W24     &
  &     - P3 * P4 *W34
  vLG(2)= (One-P2) *( P1*W12 + P3*W23 + P4*W24) &
  &     - P1 * P3 *W13     &
  &     - P1 * P4 *W14     &
  &     - P3 * P4 *W34
  vLG(3)= (One-P3) *( P1*W13 + P2*W23 + P4*W34) &
  &     - P1 * P2 *W12     &
  &     - P1 * P4 *W14     &
  &     - P2 * P4 *W24
  vLG(4)= (One-P4) *( P1*W14 + P2*W24 + P4*W34) &
  &     - P1 * P2 *W12     &
  &     - P1 * P3 *W13     &
  &     - P2 * P3 *W23
  
  ! vLG(:)= Zero
  ! do i=1,4
  !   X= Zero
  !   do j= 1,4
  !     if(j/=i) then
  !       X= X + vX(j) *W(i,j)
  !     end if
  !   end do
  !   vLG(I)= (One -vX(i)) *X
  !   vLG(I)= vLG(I) -X
  ! end do
  !
  vLG(:)= vLG(:) /R_jk/TdgK
  !
  return
end subroutine Chlorite_HPWeb_

subroutine Chlorite_HP_Ideal( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!--
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    DAPHNITE
  !    AMESITE
  !  end
  !end
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !endSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !endPOLE
  !
  !HP parameters X,Y,N=
  !X= X_Fe_M2  = Fe/(Fe+Mg)
  !Y= X_Al_T2_ =
  !N=(X_Al_M4_ - X_Al_M1_)/2
  !
  !P1=         1 - X_Al_M4_              !X(CHLORITE-ALFREE)
  !P2+P4= 2N =     X_Al_M4_ - X_Al_M1_
  !P3=                        X_Al_M1_   !X(AMESITE)
  !
  !P1= 1   -Y -N
  !P2= 2*N -2*X*(3-Y)/5
  !P3= Y   -N
  !P4= 2*vX*(3-Y_)/5
  !
  real(dp):: X,Y,Q
  real(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,MgM4,FeM4,AlM4,SiT2,AlT2
  
  X=  vX(3)
  Y= (vX(2)+vX(3))/2.0D0 +vX(4)
  Q= (vX(2)+vX(3)+vX(4)-vX(4))/2.0D0
  
  FeM2= X
  MgM2= One - X
  AlM1= Y - Q
  FeM1= X*(One-Y+Q)
  MgM1=(One-X)*(One-Y+Q)
  AlM4= Y+Q
  FeM4= X*(One-Y-Q)
  MgM4=(One-X)*(One-Y-Q)
  AlT2= Y
  SiT2= One-Y
  
  !--- ideal activities (mixing on sites)
  vActId(1)= MgM2**4 *MgM1 *MgM4 *SiT2 *SiT2
  vActId(2)= MgM2**4 *MgM1 *AlM4 *AlT2 *SiT2 *4.0D0
  vActId(3)= FeM2**4 *FeM1 *AlM4 *AlT2 *SiT2 *4.0D0
  vActId(4)= MgM2**4 *AlM1 *AlM4 *AlT2 *AlT2
  
  vLG(:)= Zero
  
  return
end subroutine Chlorite_HP_Ideal

subroutine Epidote_HP_Work( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- from Tables of Holland - Powell 2011, JMG
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: f,Q
  real(dp):: FeM1,AlM1,FeM3,AlM3
  
  !MIXTURE.MODEL EPIdoTE_WORK !Fe-free
  !  MODEL SPECIAL
  !  POLE
  !    cz
  !    ep
  !    fep
  !  end
  !end
  
  ! variables
  !   f(ep) = Fe3/(Fe3 + Al)
  !   Q(ep) = order parameter
  ! site fractions
  !   xFeM1 = f - Q
  !   xAlM1 = 1 - f + Q
  !   xFeM3 = f + Q
  !   xAlM3 = 1 - f - Q
  ! proportions
  !   cz = 1 - f - Q
  !   ep = 2Q
  !   fep = f - Q
  ! ideal mixing activities
  !   cz=  xAlM1 xAlM3
  !   ep=  xAlM1 xFeM3
  !   fep= xFeM1 xFeM3
  ! non-ideality by symmetric formalism
  !   W(cz,ep) = 1.0
  !   W(cz,fep) = 3.0
  !   W(ep,fep) = 1.0
  
  Q= vX(2)/Two  ! ep/2
  f= vX(3) +Q   ! fep +ep/2
  
  !-- site occupancy
  FeM1 =       f - Q
  AlM1 = One - f + Q
  FeM3 =       f + Q
  AlM3 = One - f - Q

  !--- ideal activities (mixing on sites)
  vActId(1)= AlM1 *AlM3  ! cz
  vActId(2)= AlM1 *FeM3  ! ep
  vActId(3)= FeM1 *FeM3  ! fep
  
  vLG(:)= Zero
  
  return
end subroutine Epidote_HP_Work

subroutine Chlorite_HP_Work( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- from Tables of Holland - Powell 2011, JMG
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: Q,Y
  real(dp):: MgM1,AlM1,MgM4,AlM4,SiT2,AlT2
  !
  !MIXTURE.MODEL CHLORITE_WORK !Fe-free
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !  end
  !end
  
  ! chlorite
  !
  ! variables
  !   y(chl) = Al/(Al + Mg)
  !   Q(chl) = order parameter
  !
  ! site fractions
  !   x(Al,M1) = y -Q/2
  !   x(Mg,M1) = 1 -y +Q/2
  !   x(Al,M4) = y +Q/2
  !   x(Mg,M4) = 1 -y -Q/2
  !   x(Al,T2) = y
  !   x(Si,T2) = 1 -y
  !
  ! proportions
  !   afchl = 1 - y - 1/2 Q
  !   clin = Q
  !   ames = y - 1/2 Q
  !
  ! ideal mixing activities
  !   afchl = MgM1 *MgM4 *SiT2 *SiT2
  !   clin =  MgM1 *AlM4 *AlT2 *SiT2 *4.0D0
  !   ames =  AlM1 *AlM4 *AlT2 *AlT2
  !
  ! non-ideality by symmetric formalism
  !   W(afchl,clin) = 17.0
  !   W(afchl,ames) = 20.0
  !   W(clin,ames) = 17.0
  
  Q= vX(2)            ! CLIN
  Y= vX(3) +vX(2)/Two ! AMES + CLIN/Two
  
  !-- site occupancy
  AlM1=      Y -Q/Two
  MgM1= One -Y +Q/Two
  AlM4=      Y +Q/Two
  MgM4= One -Y -Q/Two
  AlT2=      Y
  SiT2= One -Y

  !--- ideal activities (mixing on sites)
  vActId(1)= MgM1 *MgM4 *SiT2 *SiT2
  vActId(2)= MgM1 *AlM4 *AlT2 *SiT2 *4.0D0
  vActId(3)= AlM1 *AlM4 *AlT2 *AlT2
  
  vLG(:)= Zero
  
  return
end subroutine Chlorite_HP_Work

subroutine Biotite_HP05( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- BIO-HP- BIOTITE, Holland & Powell, David Dolejs 13-Mar-05
!-- from THERIAK-doMINO
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !    SIDEROPHYLLITE
  !  end
  !end
  !
  real(dp):: X,Y,N
  real(dp):: P1,P2,P3,P4
  !
  !--- compositional variables
  !-- X- bulk Fe/(Fe+Mg)
  !-- Y- X_Al_M1
  !-- N- (X -X_Fe_M2) *3
  X= ( 3*vX(2) +vX(4) +2*vX(5) ) &
  &  / ( 3*vX(2) +vX(4) +2*vX(5) +3*vX(1) +2*vX(3) +2*vX(4) )
  !
  Y= vX(3)+vX(5)
  !
  N= (X-vX(2)-vX(5)) *3
  !---/compositional variables

  !--- ideal mixing activities
  vActId(1)= ((One-X)*(One-Y) -Two/3.0D0*N) &
  &      * ( One-X+N/3.0D0)**2        &
  &      * ( One+Y )*( One-Y )
  !
  vActId(2)= (X*(One-Y)+Two/3.0D0*N) &
  &      * (X-N/3.0D0)**2 &
  &      * ( One+Y )*( One-Y )
  !
  vActId(3)=  Y &
  &      * (One-X+N/3.0D0)**2 &
  &      * (One+Y)**2 /4.0D0
  !
  vActId(4)= (X*(One-Y)+Two/3.0D0*N) &
  &      * (One-X+N/3.0D0)**2 &
  &      * (One+Y)*(One-Y)
  !
  vActId(5)=  Y &
  &      * (X-N/3.0D0)**2 &
  &      * (One+Y)**2 /4.0D0
  !---/ideal mixing activities
  
  !--- independent mole fractions
  P1= (One-X)*(One-Y) -N*Two/3.0D0    ! Phl
  P2= X -N/3.0D0                      !
  P3= Y
  P4=-X*Y +N
  !
  !--- RTlog(activity coefficients) --
  vLG(1)= (One-P1) &
  &       * ( P2  *9.0D3     & ! W-Phl-Ann
  &         + P4  *3.0D3     & ! W-Phl-Obi
  &         + P3  *10.0D3 )  & ! W-Phl-Eas
  &     - P3 * P2 *(-1.0D3)  & ! W-Ann-Eas
  &     - P4 * P3 *10.0D3    & ! W-Eas-Obi
  &     - P4 * P2 *6.0D3       ! W-Ann-Obi
  !
  vLG(2)= (One-P2)             &
  &       * (P1 *9.0D3         &
  &        + P4 *6.0D3         &
  &        + P3 *(-1.0D3) )    &
  &     - P4 *P1       *3.0D3  &
  &     - P4 *P3       *10.0D3 &
  &     - P3 *P1       *10.0D3
  !
  !! vLG(3)= P1 *(One-P3) *10.0D3 &
  !! &     + P2 *(One-P3) *(-1.0D3) &
  !! &     + P4 *(One-P3) *10.0D3 &
  !! &     - P1 *P2       *9.0D3  &
  !! &     - P4 *P2       *6.0D3  &
  !! &     - P4 *P1       *3.0D3
  !
  vLG(3)= (One-P3)             &
  &       *( P1 *10.0D3        &
  &        + P2 *(-1.0D3)      &
  &        + P4 *10.0D3 )      &
  &     - P1 *P2       *9.0D3  &
  &     - P4 *P2       *6.0D3  &
  &     - P4 *P1       *3.0D3
  !
  vLG(4)= P2 *(One-P4) *6.0D3  &
  &     + P1 *(One-P4) *3.0D3  &
  &     + P3 *(One-P4) *10.0D3 &
  &     - P1 *P2       *9.0D3  &
  &     - P3 *P1       *10.0D3 &
  &     - P2 *P3       *(-1.0D3)
  !
  vLG(5)= P1     *(One-P2) *9.0D3  &
  &     + P1     *(One-P3) *10.0D3 &
  &     - P1     *(One+P4) *3.0D3  &
  &     -(One-P2)*(One-P3) *(-1.0D3) &
  &     +(One-P2)*(One+P4) *6.0D3 &
  &     +(One-P3)*(One+P4) *10.0D3
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Biotite_HP05

subroutine Opx_HP05( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- OPX-HP- ORTHOPYROXENE, Holland & Powell -- David Dolejs, 13-March-05
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL OPX_HP05
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    FM.PYX
  !    MG-TSCHER.PYX
  !    FE-TSCHER.PYX
  !  end
  !end
  !
  real(dp):: XHP,YHP,QHP
  real(dp):: P1,P2,P3,P4
  !
  !--- compositional variables
  YHP= X(4)+X(5)
  XHP= ( X(2)*Two +X(3) +X(5) ) /( Two -X(5) -X(4) )
  QHP= ( X(2) +X(3) +X(5) - (X(2)*Two +X(3) +X(5)) / (Two-X(5)-X(4)) )*Two
  !
  !--- ideal mixing activities
  vActId(1)= (One-XHP -QHP/Two) *(One-YHP +QHP/Two -XHP*(One-YHP))
  vActId(2)= (XHP     +QHP/Two) *(        -QHP/Two +XHP*(One-YHP))
  vActId(3)= (XHP     +QHP/Two) *(One-YHP +QHP/Two -XHP*(One-YHP))
  vActId(4)= (One-XHP -QHP/Two) *YHP
  vActId(5)= (XHP     +QHP/Two) *YHP
  !
  !--- independent mole fractions
  P1= One -XHP -YHP -QHP/Two
  P2= XHP*(One -YHP) -QHP/Two
  P3= QHP +XHP*YHP
  P4= YHP
  !
  !--- RT*ln(activity coefficients)
  vLG(1)=   (One-P1) *P2 *6800.0 &
  &       + (One-P1) *P3 *4500.0 &
  &       -      P2  *P4 *(-1000.0) &
  &       -      P2  *P3 *4500.0 &
  &       -      P4  *P3 *1200.0
  !
  vLG(2)=    P1 *(One-P2) *6800.0 &
  &       -  P1 *     P3  *4500.0 &
  &       + (One-P2) *P4  *(-1000.0) &
  &       + (One-P2) *P3  *4500.0 &
  &       -      P4  *P3  *1200.0
  !
  vLG(3)= - P1      *P2  *6800.0 &
  &       + P1 *(One-P3) *4500.0 &
  &       - P2      *P4  *(-1000.0) &
  &       + P2 *(One-P3) *4500.0 &
  &       + P4 *(One-P3) *1200.0
  !
  vLG(4)= - P1      *P2  *6800.0 &
  &       - P1      *P3  *4500.0 &
  &       + P2 *(One-P4) *(-1000.0) &
  &       - P2      *P3  *4500.0 &
  &       +(One-P4) *P3  *1200.0
  !
  vLG(5)= - (One+P1)      *P2  *6800.0 &
  &       + (One+P1) *(One-P3) *4500.0 &
  &       +      P2  *(One-P4) *(-1000.0) &
  &       +      P2  *(One-P3) *4500 &
  &       - (One-P4) *(One-P3) *1200
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Opx_HP05

subroutine Opx_HP11( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- OPX-HP- ORTHOPYROXENE, Holland & Powell, 2011
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL OPX_HP11
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    MG-TSCHER.PYX
  !    FM.PYX
  !  end
  !end
  !
  real(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,SiT1,AlT1
  ! real(dp):: P1,P2,P3,P4
  !
  !--- site occupancy
  MgM1= X(1) +X(4)
  FeM1= X(2)
  AlM1= X(3)
  !
  MgM2= X(1) +X(3)
  FeM2= X(2) +X(4)
  !
  SiT1= X(1) +X(2) +X(3)/2.0D0 +X(4)
  AlT1= X(3)/2.0D0
  !---/site occupancy
  
  vActId(1)= MgM1 *MgM2 *SiT1**0.5D0
  vActId(2)= FeM1 *FeM2 *SiT1**0.5D0
  vActId(3)= AlM1 *MgM2 *SiT1**0.25D0 *AlT1**0.25D0 *sqrt(2.0D0)
  vActId(4)= MgM1 *FeM2 *SiT1**0.5D0
      
  vLG(:)= Zero
  !
  return
end subroutine Opx_HP11

subroutine Opx_Avchenko( &
& TdgK, &
& vX,    &
& vActId,vLG)
!--
!--
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  real(dp):: W_EnFs,W_FsMgts
  !
  !MIXTURE.MODEL OPX_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    MG-TSCHER.PYX
  !  end
  !end
  !
  !-- en   Mg Mg Si2  O6
  !-- fs   Fe Fe Si2  O6
  !-- mgts Al Mg AlSi O6
  !
  !--- ideal mixing activities
  vActId(1)= vX(1)*(vX(1)+vX(3)) ! en*(en+mgts)
  vActId(2)= vX(2)**2            ! fs**2
  vActId(3)= vX(3)*(vX(1)+vX(3)) ! mgts*(en+mgts)
  !
  W_EnFs=    6.8D3
  W_FsMgts= -1.0D3
  !
  vLG(1)= (One-vX(1))*vX(2)*W_EnFs  -vX(2)*vX(3)*W_FsMgts
  vLG(2)= (One-vX(2))*(vX(1)*W_EnFs +vX(3)*W_FsMgts)
  vLG(3)= (One-vX(3))*vX(2)*W_FsMgts -vX(2)*vX(1)*W_EnFs
  !
  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Opx_Avchenko

subroutine Mica_Vidal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- MICA4- Vidal
!--
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL MICA_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CELAdoNITE      !KMgAlSi4O10(OH)2
  !    MUSCOVITE       !KAl3Si3O10(OH)2
  !    PYROPHYLLITE    !Al2Si4O10(OH)2
  !  end
  !end
  !
  real(dp):: Si_T,Al_T,Al_O,Mg_O,K_A,VaA
  real(dp):: gamAl,gamMg,gamAlc,gamv,gamK,gamMgc
  real(dp):: W_AlMg,W_KKv,W_Kvv
  !
  !-- site     A  O     T
  !--             M2    T2
  !-- celado   K  MgAl  Si2  Si2  O10(OH)2
  !-- muscov   K  Al2   AlSi Si2  O10(OH)2
  !-- pyroph      Al2   Si2  Si2  O10(OH)2
  !
  !--- site occupancy --
  Si_T=  vX(1) +vX(2)/Two +vX(3)
  Al_T=         vX(2)/Two
  !
  Al_O=  vX(1)/Two +vX(2) +vX(3)
  Mg_O=  vX(1)/Two
  !
  K_A=    vX(1) +vX(2)
  VaA=    vX(3)
  !---/

  !--- ideal activities --
  vActId(1)= Si_T**2         *Al_O     *Mg_O  *K_A *4.0D0 ! CELAdoNITE
  vActId(2)= Si_T      *Al_T *Al_O**2         *K_A *4.0D0 ! MUSCOVITE
  vActId(3)= Si_T**2         *Al_O**2         *VaA        ! PYROPHYLLITE
  !---/

  !--- T,P-dependent interaction param's --
  W_AlMg= -30.500D3 +15.0D0*TdgK +78.0D-2*Pbar !!$$$!!! Bar of kBar ???
  W_KKv=   35.000D3 -25.0D0*TdgK -85.0D-2*Pbar
  W_Kvv=   45.000D3 -10.0D0*TdgK -85.0D-2*Pbar
  !---/

  !--- RTlngamma's --
  gamMgc= W_AlMg *(Al_O -0.5D0) *(0.5D0 -Mg_O)
  gamAlc= 0.0D0
  gamAl=  W_AlMg*Mg_O*(1.0D0-Al_O)
  gamK=   W_KKv*(2.0D0*K_A*VaA-2.0D0*K_A*K_A*VaA) + W_Kvv*(VaA*VaA-2.0D0*K_A*VaA*VaA)
  gamv=   W_Kvv*(2.0D0*K_A*VaA-2.0D0*K_A*VaA*VaA) + W_KKv*(K_A*K_A-2.0D0*K_A*K_A*VaA)
  gamK=   0.0D0
  gamv=   0.0D0
  !---/
  !
  vLG(1)= (gamAlc+gamMgc+gamK) ! gamCel
  vLG(2)= (gamAl+gamK)         ! gamMus
  vLG(3)= (gamAl+gamv)         ! gamPrl
  !
  vLG(:)= vLG(:) /R_jk/TdgK
  !
  return
end subroutine Mica_Vidal

subroutine Chlorite_Vidal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- CHLVIDAL- CHLORITE Vidal, AJS301, 2001
!--
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE   !Mg5Al2Si3O10(OH)8
  !    DAPHNITE      !Fe5Al2Si3O10(OH)8
  !    AMESITE       !Mg4Al4Si2O10(OH)8  ! = AM-VID
  !    SUdoITE       !Mg2Al4Si3O10(OH)8
  !  end
  !end
  !
  real(dp):: SiT2,AlT2,AlM1,MgM1,FeM1,VaM1,AlM2,MgM2,FeM2
  real(dp):: W_AlMg, W_AlFe, W_VaMg, W_VaAl, W_VaFe
  !
  !            T1     T2     M1   M2+M3     M4
  ! 1=CLIN     Si,Si  Si,Al  Mg   Mg4       Al   !Mg5Al2Si3O10(OH)8
  ! 2=DAPH     Si,Si  Si,Al  Fe   Fe4       Al   !Fe5Al2Si3O10(OH)8
  ! 3=AMESITE  Si,Si  Al,Al  Al   Mg4       Al   !Mg4Al4Si2O10(OH)8
  ! 4=SUdoITE  Si,Si  Si,Al  Va   Al2,Mg2   Al   !Mg2Al4Si3O10(OH)8
  !
  !--- site occupancy --
  SiT2= ( vX(1)+vX(2)+vX(4) )/Two
  AlT2= SiT2 +vX(3)
  !
  MgM1=  vX(1)
  FeM1=  vX(2)
  AlM1=  vX(3)
  VaM1=  vX(4) ! vacancy
  !
  MgM2=  vX(1) +vX(3) +vX(4)/Two
  FeM2=  vX(2)
  AlM2=  vX(4)/Two
  !---/site occupancy --

  !--- ideal activities --
  !-- cf Vidal et al, AJS 301, p.563
  vActId(1)= SiT2 *AlT2    *MgM1 *MgM2**4 *4.0D0           ! CLINOCHLORE
  vActId(2)= SiT2 *AlT2    *FeM1 *FeM2**4 *4.0D0           ! DAPHNITE
  vActId(3)=       AlT2**2 *AlM1 *MgM2**4                  ! AM-VID
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUdoITE
  !---/

  !--- interaction param's, table 2, Vidal,2001 --
  W_AlMg= -9.40D3 - 30.0D0*TdgK - 0.20D-3*Pbar
  W_AlFe=  12.0D3 + 35.0D0*TdgK - 0.50D-3*Pbar
  W_VaMg=  5.00D3 - 25.0D0*TdgK + 0.90D-3*Pbar
  W_VaAl= -10.0D3 - 30.0D0*TdgK + 0.50D-3*Pbar
  W_VaFe=  2.00D3 - 15.0D0*TdgK + 0.40D-3*Pbar
  !---/

  !--- RT*ln(Gamma)'s --
  !-- cf Vidal et al, AJS 301, p.565
  vLG(1)=   W_AlMg *AlM1 *(One-MgM1) &
  &       - W_AlFe *AlM1 *FeM1       &
  &       + W_VaMg *VaM1 *(One-MgM1) &
  &       - W_VaAl *VaM1 *AlM1       &
  &       - W_VaFe *VaM1 *FeM1
  vLG(2)= - W_AlMg *AlM1 *MgM1       &
  &       + W_AlFe *AlM1 *(One-FeM1) &
  &       - W_VaMg *VaM1 *MgM1       &
  &       - W_VaAl *VaM1 *AlM1       &
  &       + W_VaFe *VaM1 *(One-FeM1)
  vLG(3)= - W_AlMg *AlM1 *MgM1       &
  &       - W_AlFe *AlM1 *FeM1       &
  &       + W_VaMg *MgM1 *(One-VaM1) &
  &       + W_VaAl *AlM1 *(One-VaM1) &
  &       + W_VaFe *FeM1 *(One-VaM1)
  vLG(4)= + W_AlMg *MgM1 *(One-AlM1) &
  &       + W_AlFe *FeM1 *(One-AlM1) &
  &       - W_VaMg *VaM1 *MgM1       &
  &       + W_VaAl *VaM1 *(One-AlM1) &
  &       - W_VaFe *VaM1 *FeM1
  !---/

  vLG(:)= vLG(:)/R_jk/TdgK
  !
  return
end subroutine Chlorite_Vidal

subroutine Chlorite_Vidal_Ideal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- CHLVIDAL- CHLORITE Vidal, AJS301, 2001
!--
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE   !Mg5Al2Si3O10(OH)8
  !    DAPHNITE      !Fe5Al2Si3O10(OH)8
  !    AMESITE       !Mg4Al4Si2O10(OH)8  ! = AM-VID
  !    SUdoITE       !Mg2Al4Si3O10(OH)8
  !  end
  !end
  !
  real(dp):: SiT2,AlT2,AlM1,MgM1,FeM1,VaM1,AlM2,MgM2,FeM2
  ! real(dp):: W_AlMg, W_AlFe, W_VaMg, W_VaAl, W_VaFe
  !
  !            T1     T2     M1   M2+M3     M4
  ! 1=CLIN     Si,Si  Si,Al  Mg   Mg4       Al   !Mg5Al2Si3O10(OH)8
  ! 2=DAPH     Si,Si  Si,Al  Fe   Fe4       Al   !Fe5Al2Si3O10(OH)8
  ! 3=AMESITE  Si,Si  Al,Al  Al   Mg4       Al   !Mg4Al4Si2O10(OH)8
  ! 4=SUdoITE  Si,Si  Si,Al  Va   Al2,Mg2   Al   !Mg2Al4Si3O10(OH)8
  !
  !--- site occupancy --
  SiT2= ( vX(1)+vX(2)+vX(4) )/Two
  AlT2= SiT2 +vX(3)
  !
  MgM1=  vX(1)
  FeM1=  vX(2)
  AlM1=  vX(3)
  VaM1=  vX(4) ! vacancy
  !
  MgM2=  vX(1) +vX(3) +vX(4)/Two
  FeM2=  vX(2)
  AlM2=  vX(4)/Two
  !---/site occupancy --

  !--- ideal activities --
  !-- cf Vidal et al, AJS 301, p.563
  vActId(1)= SiT2 *AlT2    *MgM1 *MgM2**4 *4.0D0           ! CLINOCHLORE
  vActId(2)= SiT2 *AlT2    *FeM1 *FeM2**4 *4.0D0           ! DAPHNITE
  vActId(3)=       AlT2**2 *AlM1 *MgM2**4                  ! AM-VID
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUdoITE
  !---/

  vLG(:)= Zero
  
  return
end subroutine Chlorite_Vidal_Ideal

subroutine Chlorite_Fe( &
& TdgK,  &
& vX,    &
& vActId,vLGam)
!--
!-- CHLFe- Vidal
!--
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL CHLORITE_FE
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE
  !    DAPHNITE
  !    FE-AMESITE
  !    SUdoITE
  !  end
  !end
  !
  !           T1     T2     M1   M2+M3     M4
  ! 1=CLIN    Si,Si  Si,Al  Mg   Mg4       Al
  ! 2=DAPH    Si,Si  Si,Al  Fe   Fe4       Al
  ! 3=FE-AM   Si,Si  Al,Al  Al   Fe4       Al
  ! 4=MG-SUD  Si,Si  Si,Al  Va   Al2,Mg2   Al
  !
  real(dp):: SiT2,AlT2
  real(dp):: AlM1,MgM1,FeM1,VaM1
  real(dp):: AlM2,MgM2,FeM2
  ! real(dp):: SiC,AltC,AlM1C,MgM1C,FeM1C,vM1C,AlM2C,MgM2C,FeM2C
  !
  !--- site occupancy --
  SiT2= ( vX(1)+vX(2)+vX(4) )/Two
  AlT2= SiT2 +vX(3)
  !
  MgM1=  vX(1)
  FeM1=  vX(2)
  AlM1=  vX(3)
  VaM1=  vX(4)
  !
  MgM2= vX(1) +vX(4)/Two
  FeM2= vX(2) +vX(3)
  AlM2= vX(4)/Two
  !---/site occupancy --

  !--- ideal activities --
  vActId(1)= SiT2 *AlT2    *MgM1 *MgM2**4 *4.0D0           ! CLINOCHLORE
  vActId(2)= SiT2 *AlT2    *FeM1 *FeM2**4 *4.0D0           ! DAPHNITE
  vActId(3)=       AlT2**2 *AlM1 *FeM2**4                  ! FE-AMESITE
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUdoITE
  !
  vLGam(:)= Zero
  !
  return
end subroutine Chlorite_Fe

subroutine Chlorite_HP98( & !
& TdgK, &
& vX,   &
& vActId,vLGam)
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !    DAPHNITE
  !  end
  !end
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !endSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !endPOLE
  !
  !HP parameters X,Y,N=
  !X= X_Fe_M2  = Fe/(Fe+Mg)
  !Y= X_Al_T2_ =
  !N=(X_Al_M4_ - X_Al_M1_)/2
  !
  !vX(1)=        1 - X_Al_M4_                  !X(CHLORITE-ALFREE)
  !vX(2)+vX(4)= 2N = X_Al_M4_ - X_Al_M1_
  !vX(3)=            X_Al_M1_                  !X(AMESITE)
  !
  !vX(1)= 1   -Y -N
  !vX(2)= 2*N -2*X*(3-Y)/5
  !vX(3)= Y   -N
  !vX(4)= 2*vX*(3-Y_)/5
  !
  real(dp):: X,Y,N
  !real(dp):: G(4)
  !real(dp):: W(4,4)
  !integer :: I,J,K
  !W daph ames
  !clin 2.5 18
  !daph 20.5
  !
  Y=(One +vX(3) -vX(1))/Two
  N=(One -vX(3) -vX(1))/Two
  X= vX(4) / (One - (vX(3)-vX(1))/5.0D0)
  !
  !mu(i)= mu°(i) + RT.ln(a_id(i)) + RTLnGam(i)
  !
  !--- ideal activities (mixing on sites)
  vActId(1)= (One-X)**6 *(One-Y+N)*(One-Y-N)*(One-Y)*(One-Y)   ! CHLORITE-ALFREE
  vActId(2)= (One-X)**5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! CLINOCHLORE
  vActId(3)= (One-X)**4 *(    Y-N)*(    Y+N)*Y      *Y         ! AMESITE
  vActId(4)=      X **5 *(One-Y+N)*(    Y+N)*(One-Y)*Y  *4.0D0 ! DAPHNITE
  !
  vLGam(:)= Zero
  !
  !RTLnGam(i)= sum( xj.(1-xi).tDtbMltWij ) - sum( xj.xk .Wjk )
  !            j/=i                   j/=i & k/=i
  !!
  !! if(S%NMarg>0) then
  !!   W=Zero
  !!   ! retrieve values of W(I,J) from DeCapitani formattted coeffs
  !!   ! this code restricted to binary 12 margules coeffs,
  !!   ! which is general case in HP98 models
  !!   !
  !!   do I=1,S%NMarg
  !!     M= S%vMarg(I)
  !!     !! W(M%vIPole(1),M%vIPole(2))= M%WG
  !!   end do
  !!   !
  !!   do I=1,4
  !!     G(I)=Zero
  !!     do J=1,4
  !!       if(J/=I) G(I)=G(I) + X(J)*(One-X(I))*W(I,J)
  !!     end do
  !!     do J=1,4
  !!       do K=1,4
  !!         if(J/=I .and. K/=I) G(I)=G(I) - X(J)*X(K)*W(J,K)
  !!       end do
  !!     end do
  !!     F%vLMarg(I)= G(I)/R_jk/TdgK
  !!     F%vLGam(I)=  F%vLGam(I)+F%vLMarg(I)
  !!   end do
  !! end if
  !
  return
end subroutine Chlorite_HP98

subroutine Biotite_HP98( &
& TdgK, &
& X,    &
& vActId,vLGam)
  real(dp),intent(in) :: TdgK
  real(dp),intent(in) :: X(:)
  real(dp),intent(out):: vActId(:)
  real(dp),intent(out):: vLGam(:)
  !
  !MIXTURE.MODEL BIOTITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !  end
  !end
  !
  !biotite, KFMASH
  !Fe-Mg ordering -> M1 the Fe-loving site, M2 the other
  !OBIOTITE=PHLOGOPITE[2/3]ANNITE[1/3]
  !SITE
  !  M1 1 MG_FE_AL_
  !  M2 2 MG_FE_
  !  T1 2 SI_AL_
  !endSITE
  !POLE
  !   PHLOGOPITE     MG_   MG_MG_  AL_SI_
  !   ANNITE         FE_   FE_FE_  AL_SI_
  !   EASTONITE      AL_   MG_MG_  AL_AL_
  !   OBIOTITE       FE_   MG_MG_  AL_SI_
  !endPOLE
  !
  real(dp):: X1,X2,Y
  !real(dp):: W(4,4)
  !integer :: I,J,K
  !
  !IN=X1,X2,Y -> OUT= site distribution
  !X_Mg_M1_=(One-Y)*X1;   X_Fe_M1_=(One-Y)*(One-X1);  X_Al_M1_=Y
  !X_Mg_M2_= One-X2;       X_Fe_M2_= X2
  !X_Si_T1_=(One-Y)/2.D0; X_Al_T1_=(One+Y)/2.D0;
  !
  Y=  X(3)  !X(3)= X_Al_M1_= Y
  X2= X(2)  !X(2)= X_Fe_M2_= X2
  X1= (One -X(3)-X(1)) / (One -X(3))
  !
  vActId(1)= (One-Y)**2 *(One-X1) *(One-X2)**2 *(One+Y)           !PHLOGOPITE
  vActId(2)= (One-Y)**2 *     X1  *     X2     *(One+Y)           !ANNITE
  vActId(3)=                       (One-X2)**2 *(One+Y)**2 /4.0D0 !EASTONITE
  vActId(4)= (One-Y)**2 *     X1  *(One-X2)**2 *(One+Y)           !OBIOTITE
  !
  vLGam(:)= Zero
  !
  !! if(S%NMarg>0) then
  !!   W=Zero
  !!   !retrieve values of W(I,J) from DeCapitani formattted coeffs
  !!   !this code restricted to binary 12 margules coeffs,
  !!   !which is general case in HP98 models
  !!   do I=1,S%NMarg
  !!     M=S%vMarg(I)
  !!     ! W(M%vIPole(1),M%vIPole(2))=M%WG
  !!   end do
  !!   !RTLnGam(i)= sum( xj.(1-xi).tDtbMltWij ) - sum( xj.xk .Wjk )
  !!   !            j/=i                   j/=i & k/=i
  !!   do I=1,4
  !!     G=Zero
  !!     do J=1,4
  !!       if(J/=I) G= G + X(J)*(One-X(I))*W(I,J)
  !!     end do
  !!     do J=1,4
  !!       do K=1,4
  !!         if(J/=I .and. K/=I) G= G - X(J)*X(K)*W(J,K)
  !!       end do
  !!     end do
  !!     F%vLMarg(I)= G /R_jk/TdgK
  !!     F%vLGam(I)=  F%vLGam(I)+F%vLMarg(I)
  !!   end do
  !! end if
  !
  return
end subroutine Biotite_HP98

end module M_MixModel_Special

