MODULE M_MixModel_Special
  USE M_Kinds
  USE M_Dtb_Const, ONLY: R_jk
  USE M_Trace,ONLY: Stop_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: MixModel_Special_Activities

CONTAINS

SUBROUTINE MixModel_Special_Activities( &
& sKey,   &
& TdgK, Pbar, &
& vX,     &
& vLPole, &
& vLActIdl,  &
& vLGam)
  CHARACTER(LEN=*),INTENT(IN) :: sKey
  REAL(dp),INTENT(IN) :: TdgK, Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  LOGICAL, INTENT(IN) :: vLPole(:)
  REAL(dp),INTENT(OUT):: vLActIdl(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  INTEGER:: N
  REAL(dp),ALLOCATABLE:: vActId(:)
  !
  N= SIZE(vX)
  ALLOCATE(vActId(N))
  !
  SELECT CASE(TRIM(sKey))

    CASE("SPECIAL_2")        ;  CALL Special_2(TdgK,vX,vActId,vLGam)
    CASE("SPECIAL_3")        ;  CALL Special_3(TdgK,vX,vActId,vLGam)

    CASE("FLUID_KERJAC")     ;  CALL Fluid_KerJac(TdgK,Pbar,vX,vActId,vLGam)

    CASE("FELSPAR_GH84")     ;  CALL Felspar_Gh84(TdgK,vX,vActId,vLGam)
    CASE("FELSPAR_HP03A")    ;  CALL Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)
    CASE("FELSPAR_HP03B")    ;  CALL Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)
    CASE("FELSPAR_HP03C")    ;  CALL Felspar_HP03(TdgK,Pbar,vX,vActId,vLGam)

    CASE("CORDIERITE_HP05")  ;  CALL Cordierite_HP05(TdgK,vX,vActId,vLGam)
    CASE("CORDIERITE_HPWEB") ;  CALL Cordierite_HPWeb(TdgK,vX,vActId,vLGam)

    CASE("BIOTITE_HP05")     ;  CALL Biotite_HP05(TdgK,vX,vActId,vLGam)
    CASE("BIOTITE_HP98")     ;  CALL Biotite_HP98(TdgK,vX,vActId,vLGam)
    CASE("BIOTITE_HPWEB")    ;  CALL Biotite_HPWeb(TdgK,vX,vActId,vLGam)
    CASE("MUSCOVITE_HP05")   ;  CALL Muscovite_HPWeb(TdgK,Pbar,vX,vActId,vLGam)
    CASE("MICA_VIDAL")       ;  CALL Mica_Vidal(TdgK,Pbar,vX,vActId,vLGam)

    CASE("OPX_HP05")         ;  CALL Opx_HP05(TdgK,vX,vActId,vLGam)
    CASE("CHLORITE_VIDAL")   ;  CALL Chlorite_Vidal(TdgK,Pbar,vX,vActId,vLGam)
    CASE("CHLORITE_VIDAL2")  ;  CALL Chlorite_Vidal_Ideal(TdgK,Pbar,vX,vActId,vLGam)
    CASE("CHLORITE_FE")      ;  CALL Chlorite_Fe(TdgK,vX,vActId,vLGam)
    CASE("CHLORITE_HP05")    ;  CALL Chlorite_HPWeb(TdgK,vX,vActId,vLGam)

    !--- AVCHENKO MODELS
    CASE("BIOTITE_AVCH")     ;  CALL Biotite_Avchenko(TdgK,vX,vActId,vLGam)
    CASE("CORDIERITE_AVCH")  ;  CALL Cordierite_Avchenko(TdgK,vX,vActId,vLGam)
    CASE("OPX_AVCH")         ;  CALL Opx_Avchenko(TdgK,vX,vActId,vLGam)
    CASE("TALC_AVCH")        ;  CALL Talc_Avchenko(TdgK,vX,vActId,vLGam)
    CASE("EPIDOTE_AVCH")     ;  CALL Epidote_Avchenko(TdgK,vX,vActId,vLGam)
    CASE("ALKSPAR_AVCH")     ;  CALL Alkspar_Avchenko(TdgK,Pbar,vX,vActId,vLGam)
    CASE("PLAGIO_AVCH")      ;  CALL Plagio_Avchenko(TdgK,vX,vActId,vLGam)
    !---/AVCHENKO MODELS

    CASE DEFAULT
      CALL Stop_(TRIM(sKey)//" is NOT A SPECIAL MODEL")

  END SELECT
  !
  WHERE(vLPole(:))  ;  vLActIdl(:)= LOG(vActId(:))
  ELSEWHERE         ;  vLActIdl(:)= Zero
  END WHERE
  !
  DEALLOCATE(vActId)
  !
  RETURN
END SUBROUTINE MixModel_Special_Activities

SUBROUTINE Special_2( &
& TdgK,               &
& X,                  &
& vActId,vLGam)
  !MIXTURE.MODEL SPECIAL_01
  !  MODEL SPECIAL
  !  POLE
  !    AA_M1
  !    AB_M1
  !  END
  !END
  !---------------------------------------------------------------------
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !---------------------------------------------------------------------
  REAL(dp):: W12
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
END SUBROUTINE Special_2

SUBROUTINE Special_3( &
& TdgK, &
& X,    &
& vActId,vLGam)
  !MIXTURE.MODEL SPECIAL_03
  !  MODEL SPECIAL
  !  POLE
  !    MIN1
  !    MIN2
  !    MIN3
  !  END
  !END
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  REAL(dp):: W12,W23,W31
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
END SUBROUTINE Special_3

SUBROUTINE Fluid_KerJac(TdgK,Pbar,vX,vActId,vLGam)
  USE M_Eos_KerJac
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)

  vActId(:)= vX(:)
  vLGam(:)=  Zero !vLGam(:) /R_jk/TdgK

  RETURN
END SUBROUTINE

SUBROUTINE Plagio_Wrk( &
& TdgK, &
& X,    &
& vActId,vLGam)
!--
!-- G-FSP- Ghiorso, 1984
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL FELSPAR_WORK
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE
  !    ANORTHITE
  !  END
  !END

  vActId(1)= X(1) *( One -X(2)**2)       ! albite
  vActId(2)= X(2) *((One +X(2))/2.0D0)**2 ! anorthite

  vLGam(1)= X(2)**2 *( 28211.0974D0 -39485.4998D0*X(1) )
  vLGam(2)= X(1)**2 *( 8468.3475D0  +39485.4998D0*X(2) )
  vLGam(3)= Zero
  !
  vLGam(:)= vLGam(:) /R_jk/TdgK

  RETURN
END SUBROUTINE Plagio_Wrk

SUBROUTINE Felspar_Gh84( &
& TdgK, &
& X,    &
& vActId,vLGam)
!--
!-- G-FSP- Ghiorso, 1984
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL FELSPAR_GH84
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE
  !    ANORTHITE
  !    SANIDINE
  !  END
  !END
  !
  vActId(1)= X(1) *(Two -X(1))
  vActId(2)= X(2) *(One +X(2)) /4.0D0
  vActId(3)= X(3)
  !
  vLGam(1)= X(2)**2 *( 28211.0974D0 -39485.4998D0*X(1) )
  vLGam(2)= X(1)**2 *( 8468.3475D0  +39485.4998D0*X(2) )
  vLGam(3)= Zero
  !
  vLGam(:)= vLGam(:) /R_jk/TdgK
  !
  RETURN
END SUBROUTINE Felspar_Gh84

SUBROUTINE Felspar_HP03( &
& TdgK, Pbar, &
& X, &
& vActId,vLG)
!--
!-- (Holland and Powell (2003)) -- David Dolejs, 13-March-05
!--
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL FELSPAR_HP03A
  !  MODEL SPECIAL
  !  POLE
  !    MICROCLINE
  !    ALBITE
  !    ANORTHITE-C1
  !  END
  !END
  !! ..OR..
  !MIXTURE.MODEL FELSPAR_HP03B
  !  MODEL SPECIAL
  !  POLE
  !    SANIDINE
  !    ALBITE
  !    ANORTHITE-C1
  !  END
  !END
  !! ..OR..
  !MIXTURE.MODEL FELSPAR_HP03C
  !  MODEL SPECIAL
  !  POLE
  !    SANIDINE
  !    HIGH-ALBITE
  !    ANORTHITE-C1
  !  END
  !END
  !
  REAL(dp):: V1,V2,V3,F1,F2,F3,P
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
  !LG(1)=X(1)*DEXP(A1/(R*T))
  !LG(2)=X(2)*DEXP(A2/(R*T))
  !LG(3)=X(3)*DEXP(A3/(R*T))
  !
  vActId(:)= X(:)
  vLG(:)= vLG(:)/R_jk/TdgK
  !
END SUBROUTINE Felspar_HP03

SUBROUTINE Plagio_Avchenko( &
& TdgK, &
& vX,   &
& vActId,vLG)
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: Xb,WI,WC,Iab,Ian
  
  !MIXTURE.MODEL PLAGIO_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE-LOW
  !    ANORTHITE
  !  END
  !END
  
  vActId(:)= vX(:)
  
  WC= 1.07D3
  WI= 9.79D3
  Xb= 0.12D0 + TdgK *0.00038D0
  Ian= (WI-WC) *(One-Xb)**2
  Iab= (WC-WI) *Xb**2
  
  IF(vX(2) > Xb) THEN
    vLG(1)= WI *vX(2)**2 +Iab
    vLG(2)= WI *vX(1)**2
  ELSE
    vLG(1)= WC *vX(2)**2
    vLG(2)= WC *vX(1)**2 +Ian
  ENDIF
  
  vLG(:)= vLG(:)/R_jk/TdgK
  
  RETURN
END SUBROUTINE Plagio_Avchenko

SUBROUTINE Alkspar_Avchenko( &
& TdgK,Pbar, &
& vX,   &
& vActId,vLG)
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: W_Na,W_K
  !
  !MIXTURE.MODEL ALKSPAR_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ALBITE-LOW
  !    MICROCLINE
  !  END
  !END
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
  RETURN
END SUBROUTINE Alkspar_Avchenko

SUBROUTINE Cordierite_HPWeb( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/thermocalc-cordierites
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL CORDIERITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !  END
  !END
  !
  vActId(1)= ( vX(1)+vX(3) )**2 *( vX(1)+vX(2) ) ! X_Mg**2 *(1-X_H2O)
  vActId(2)=   vX(2)        **2 *( vX(1)+vX(2) ) ! X_Fe**2 *(1-X_H2O)
  vActId(3)= ( vX(1)+vX(3) )**2 *  vX(3)         ! X_Mg**2 *X_H2O
  !
  vLGam(:)= Zero
  !
END SUBROUTINE Cordierite_HPWeb

SUBROUTINE Cordierite_HP05( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!-- CORDIERITE, IDEAL RECIPROCAL -- David Dolejs, 13-March-05
!-- code derived from THERIAK/fsol.f90
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL CORDIERITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !    H.FE-CORDIERITE
  !  END
  !END
  !
  vActId(1)= ( vX(1)+vX(3) )**2 *( vX(1)+vX(2) )
  vActId(2)= ( vX(2)+vX(4) )**2 *( vX(1)+vX(2) )
  vActId(3)= ( vX(1)+vX(3) )**2 *( vX(3)+vX(4) )
  vActId(4)= ( vX(2)+vX(4) )**2 *( vX(3)+vX(4) )
  !
  vLGam= Zero
  !
END SUBROUTINE Cordierite_HP05

SUBROUTINE Cordierite_Avchenko( &
& TdgK, &
& vX,   &
& vActId,vLGam)
!--
!--
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  REAL(dp):: W_CrdFcrd
  !
  !MIXTURE.MODEL CORDIERITE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    CORDIERITE
  !    FE-CORDIERITE
  !    HYDR.CORDIERITE
  !  END
  !END
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
END SUBROUTINE Cordierite_Avchenko

SUBROUTINE BiotiteMgAl_HPWeb( &
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
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    EASTONITE
  !  END
  !END
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
  RETURN
ENDSUBROUTINE BiotiteMgAl_HPWeb

SUBROUTINE Biotite_HPWeb( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/biotite
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !  END
  !END
  !
  REAL(dp):: Phl,Ann,Eas,Obi
  REAL(dp):: W_PhlAnn,W_PhlObi,W_PhlEas,W_AnnEas,W_EasObi,W_AnnObi
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
  RETURN
ENDSUBROUTINE Biotite_HPWeb

SUBROUTINE Biotite_Avchenko( &
& TdgK, &
& vX,    &
& vActId,vLG)
!--
!--
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !  END
  !END
  !
  REAL(dp):: Phl,Ann,Eas
  REAL(dp):: W_PhlAnn,W_PhlEas,W_AnnEas
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
  RETURN
ENDSUBROUTINE Biotite_Avchenko

SUBROUTINE Muscovite_HPWeb( &
& TdgK,Pbar, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/muscovite
!--
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: X,Y,Z
  REAL(dp):: Mus,Par,Cel,Fcl
  REAL(dp):: W_MusPar, W_CelPar, W_ParFcl, I_Par
  !
  !MIXTURE.MODEL MUSCOVITE_HPWEB
  !  MODEL SPECIAL
  !  POLE
  !    MUSCOVITE
  !    PARAGONITE
  !    CELADONITE
  !    FE-CELADONITE
  !  END
  !END
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
  RETURN
END SUBROUTINE Muscovite_HPWeb

SUBROUTINE Talc_HPWeb( & !
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
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL TALC_HP
  !  MODEL SPECIAL
  !  POLE
  !    TALC
  !    FE-TALC
  !    TSCHERMAK-TALC
  !  END
  !END
  !
  REAL(dp):: X,Y
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
  RETURN
END SUBROUTINE Talc_HPWeb

SUBROUTINE Talc_Avchenko( & !
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
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  REAL(dp):: ta,fta,tats
  !
  !MIXTURE.MODEL TALC_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    TALC
  !    FE-TALC
  !    TSCHERMAK-TALC
  !  END
  !END
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
  RETURN
END SUBROUTINE Talc_Avchenko

SUBROUTINE Epidote_Avchenko( & !
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
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: Epi,Clz,Fep
  REAL(dp):: W_EpiFep,W_ClzFep
  !
  !MIXTURE.MODEL EPIDOTE_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    EPIDOTE
  !    CLINOZOISITE
  !    FE-EPIDOTE
  !  END
  !END
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
  RETURN
END SUBROUTINE Epidote_Avchenko

SUBROUTINE Chlorite_HPWeb( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/chlorite
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite 
  !    CLINOCHLORE  
  !    AMESITE      
  !    DAPHNITE     
  !  END
  !END
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_ 
  !ENDSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite 
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !ENDPOLE
  !
  !HP PARAMETERs X,Y,N=
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
  REAL(dp):: X,Y,N
  REAL(dp):: P1,P2,P3,P4
  REAL(dp):: W12,W13,W14,W23,W24,W34
  !
  !! INTEGER :: I,J,K
  !! REAL(dp):: tW(4,4),G
  !! !
  !! tW(:,:)= Zero
  !! tW(1,2)= 18.0D3  ! Afch-Clin
  !! tW(1,3)= 20.0D3  ! Afch-Ames
  !! tW(1,4)= 14.5D3  ! Afch-Daph
  !! tW(2,3)= 18.0D0  ! Clin-Ames
  !! tW(2,4)= 2.50D3  ! Clin-Daph
  !! tW(3,4)= 13.5D3  ! Ames-Daph
  !! !
  !! DO I=1,4
  !!   G= Zero
  !!   DO J=1,4
  !!     IF(J/=I) G= G + vX(J)*tW(I,J)
  !!   ENDDO
  !!   G= G *(One-vX(I))
  !!   DO J=1,4
  !!     DO K=1,4
  !!       IF(J/=I .AND. K/=I) G= G - vX(J)*vX(K)*tW(J,K)
  !!     ENDDO
  !!   ENDDO
  !!   vLG(I)= G
  !! ENDDO
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
  !~ print *,"P1..P4=",P1,P2,P3,P4
  !~ pause
  Y=(One +P3 -P1)/Two
  N=(One -P3 -P1)/Two
  !X= 5*P4 /(5*P4+6*P1+5*P2+4*P3)= 5*P4 /(5 + P1- P3)
  X= P4 / (One + (P1-P3)/5.0D0)
  !
  N= (P2+P4)/Two
  Y= N + P3
  X= P4
  !
  !~ print *,"X,Y,Z=",X,Y,N
  !~ pause
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
  RETURN
ENDSUBROUTINE Chlorite_HPWeb

SUBROUTINE Chlorite_HPWeb_( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- http://wserv2.esc.cam.ac.uk/research/research-groups/holland/ &
!-- & thermocalc/chlorite
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !    DAPHNITE
  !  END
  !END
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !END
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !END
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
  ! INTEGER :: i,j
  ! REAL(dp):: W(4,4)
  
  REAL(dp):: X,Y,N
  REAL(dp):: P1,P2,P3,P4
  REAL(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,MgM4,FeM4,AlM4,SiT2,AlT2
  REAL(dp):: W12,W13,W14,W23,W24,W34
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
  ! DO i=1,4
  !   X= Zero
  !   DO j= 1,4
  !     IF(j/=i) THEN
  !       X= X + vX(j) *W(i,j)
  !     END IF
  !   END DO
  !   vLG(I)= (One -vX(i)) *X
  !   vLG(I)= vLG(I) -X
  ! END DO
  !
  vLG(:)= vLG(:) /R_jk/TdgK
  !
  RETURN
ENDSUBROUTINE Chlorite_HPWeb_

SUBROUTINE Chlorite_HP_Ideal( & !
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!--
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    DAPHNITE
  !    AMESITE
  !  END
  !END
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !ENDSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !ENDPOLE
  !
  !HP PARAMETERs X,Y,N=
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
  REAL(dp):: X,Y,Q
  REAL(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,MgM4,FeM4,AlM4,SiT2,AlT2
  
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
  
  RETURN
ENDSUBROUTINE Chlorite_HP_Ideal

SUBROUTINE Epidote_HP_Work( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- from Tables of Holland - Powell 2011, JMG
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: f,Q
  REAL(dp):: FeM1,AlM1,FeM3,AlM3
  
  !MIXTURE.MODEL EPIDOTE_WORK !Fe-free
  !  MODEL SPECIAL
  !  POLE
  !    cz
  !    ep
  !    fep
  !  END
  !END
  
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
  
  RETURN
END SUBROUTINE Epidote_HP_Work

SUBROUTINE Chlorite_HP_Work( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- from Tables of Holland - Powell 2011, JMG
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: Q,Y
  REAL(dp):: MgM1,AlM1,MgM4,AlM4,SiT2,AlT2
  !
  !MIXTURE.MODEL CHLORITE_WORK !Fe-free
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !  END
  !END
  
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
  
  RETURN
END SUBROUTINE Chlorite_HP_Work

SUBROUTINE Biotite_HP05( &
& TdgK, &
& vX,   &
& vActId,vLG)
!--
!-- BIO-HP- BIOTITE, Holland & Powell, David Dolejs 13-Mar-05
!-- from THERIAK-DOMINO
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL BIOTITE_HP05
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !    SIDEROPHYLLITE
  !  END
  !END
  !
  REAL(dp):: X,Y,N
  REAL(dp):: P1,P2,P3,P4
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
  !~ vLG(3)= P1 *(One-P3) *10.0D3 &
  !~ &     + P2 *(One-P3) *(-1.0D3) &
  !~ &     + P4 *(One-P3) *10.0D3 &
  !~ &     - P1 *P2       *9.0D3  &
  !~ &     - P4 *P2       *6.0D3  &
  !~ &     - P4 *P1       *3.0D3
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
  RETURN
ENDSUBROUTINE Biotite_HP05

SUBROUTINE Opx_HP05( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- OPX-HP- ORTHOPYROXENE, Holland & Powell -- David Dolejs, 13-March-05
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL OPX_HP05
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    FM.PYX
  !    MG-TSCHER.PYX
  !    FE-TSCHER.PYX
  !  END
  !END
  !
  REAL(dp):: XHP,YHP,QHP
  REAL(dp):: P1,P2,P3,P4
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
  RETURN
END SUBROUTINE Opx_HP05

SUBROUTINE Opx_HP11( &
& TdgK, &
& X,    &
& vActId,vLG)
!--
!-- OPX-HP- ORTHOPYROXENE, Holland & Powell, 2011
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL OPX_HP11
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    MG-TSCHER.PYX
  !    FM.PYX
  !  END
  !END
  !
  REAL(dp):: MgM1,FeM1,AlM1,MgM2,FeM2,SiT1,AlT1
  ! REAL(dp):: P1,P2,P3,P4
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
  vActId(3)= AlM1 *MgM2 *SiT1**0.25D0 *AlT1**0.25D0 *SQRT(2.0D0)
  vActId(4)= MgM1 *FeM2 *SiT1**0.5D0
      
  vLG(:)= Zero
  !
  RETURN
END SUBROUTINE Opx_HP11

SUBROUTINE Opx_Avchenko( &
& TdgK, &
& vX,    &
& vActId,vLG)
!--
!--
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  REAL(dp):: W_EnFs,W_FsMgts
  !
  !MIXTURE.MODEL OPX_AVCH
  !  MODEL SPECIAL
  !  POLE
  !    ENSTATITE
  !    FERROSILITE
  !    MG-TSCHER.PYX
  !  END
  !END
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
  RETURN
END SUBROUTINE Opx_Avchenko

SUBROUTINE Mica_Vidal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- MICA4- Vidal
!--
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL MICA_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CELADONITE      !KMgAlSi4O10(OH)2
  !    MUSCOVITE       !KAl3Si3O10(OH)2
  !    PYROPHYLLITE    !Al2Si4O10(OH)2
  !  END
  !END
  !
  REAL(dp):: Si_T,Al_T,Al_O,Mg_O,K_A,VaA
  REAL(dp):: gamAl,gamMg,gamAlc,gamv,gamK,gamMgc
  REAL(dp):: W_AlMg,W_KKv,W_Kvv
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
  vActId(1)= Si_T**2         *Al_O     *Mg_O  *K_A *4.0D0 ! CELADONITE
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
  RETURN
END SUBROUTINE Mica_Vidal

SUBROUTINE Chlorite_Vidal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- CHLVIDAL- CHLORITE Vidal, AJS301, 2001
!--
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE   !Mg5Al2Si3O10(OH)8
  !    DAPHNITE      !Fe5Al2Si3O10(OH)8
  !    AMESITE       !Mg4Al4Si2O10(OH)8  ! = AM-VID
  !    SUDOITE       !Mg2Al4Si3O10(OH)8
  !  END
  !END
  !
  REAL(dp):: SiT2,AlT2,AlM1,MgM1,FeM1,VaM1,AlM2,MgM2,FeM2
  REAL(dp):: W_AlMg, W_AlFe, W_VaMg, W_VaAl, W_VaFe
  !
  !            T1     T2     M1   M2+M3     M4
  ! 1=CLIN     Si,Si  Si,Al  Mg   Mg4       Al   !Mg5Al2Si3O10(OH)8
  ! 2=DAPH     Si,Si  Si,Al  Fe   Fe4       Al   !Fe5Al2Si3O10(OH)8
  ! 3=AMESITE  Si,Si  Al,Al  Al   Mg4       Al   !Mg4Al4Si2O10(OH)8
  ! 4=SUDOITE  Si,Si  Si,Al  Va   Al2,Mg2   Al   !Mg2Al4Si3O10(OH)8
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
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUDOITE
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
  RETURN
END SUBROUTINE Chlorite_Vidal

SUBROUTINE Chlorite_Vidal_Ideal( &
& TdgK,Pbar, &
& vX,    &
& vActId,vLG)
!--
!-- CHLVIDAL- CHLORITE Vidal, AJS301, 2001
!--
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLG(:)
  !
  !MIXTURE.MODEL CHLORITE_VIDAL
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE   !Mg5Al2Si3O10(OH)8
  !    DAPHNITE      !Fe5Al2Si3O10(OH)8
  !    AMESITE       !Mg4Al4Si2O10(OH)8  ! = AM-VID
  !    SUDOITE       !Mg2Al4Si3O10(OH)8
  !  END
  !END
  !
  REAL(dp):: SiT2,AlT2,AlM1,MgM1,FeM1,VaM1,AlM2,MgM2,FeM2
  ! REAL(dp):: W_AlMg, W_AlFe, W_VaMg, W_VaAl, W_VaFe
  !
  !            T1     T2     M1   M2+M3     M4
  ! 1=CLIN     Si,Si  Si,Al  Mg   Mg4       Al   !Mg5Al2Si3O10(OH)8
  ! 2=DAPH     Si,Si  Si,Al  Fe   Fe4       Al   !Fe5Al2Si3O10(OH)8
  ! 3=AMESITE  Si,Si  Al,Al  Al   Mg4       Al   !Mg4Al4Si2O10(OH)8
  ! 4=SUDOITE  Si,Si  Si,Al  Va   Al2,Mg2   Al   !Mg2Al4Si3O10(OH)8
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
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUDOITE
  !---/

  vLG(:)= Zero
  
  RETURN
END SUBROUTINE Chlorite_Vidal_Ideal

SUBROUTINE Chlorite_Fe( &
& TdgK,  &
& vX,    &
& vActId,vLGam)
!--
!-- CHLFe- Vidal
!--
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL CHLORITE_FE
  !  MODEL SPECIAL
  !  POLE
  !    CLINOCHLORE
  !    DAPHNITE
  !    FE-AMESITE
  !    SUDOITE
  !  END
  !END
  !
  !           T1     T2     M1   M2+M3     M4
  ! 1=CLIN    Si,Si  Si,Al  Mg   Mg4       Al
  ! 2=DAPH    Si,Si  Si,Al  Fe   Fe4       Al
  ! 3=FE-AM   Si,Si  Al,Al  Al   Fe4       Al
  ! 4=MG-SUD  Si,Si  Si,Al  Va   Al2,Mg2   Al
  !
  REAL(dp):: SiT2,AlT2
  REAL(dp):: AlM1,MgM1,FeM1,VaM1
  REAL(dp):: AlM2,MgM2,FeM2
  ! REAL(dp):: SiC,AltC,AlM1C,MgM1C,FeM1C,vM1C,AlM2C,MgM2C,FeM2C
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
  vActId(4)= SiT2 *AlT2    *VaM1 *MgM2**2 *AlM2**2 *64.0D0 ! SUDOITE
  !
  vLGam(:)= Zero
  !
  RETURN
END SUBROUTINE Chlorite_Fe

SUBROUTINE Chlorite_HP98( & !
& TdgK, &
& vX,   &
& vActId,vLGam)
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL CHLORITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    CHLORITE-AF !Al-free-chlorite
  !    CLINOCHLORE
  !    AMESITE
  !    DAPHNITE
  !  END
  !END
  !
  !MODEL CHLORITE_HP98 !SITE
  !SITE
  !  M2 4 MG_FE_
  !  M1 1 MG_FE_AL_
  !  M4 1 MG_AL_
  !  T2 2 AL_SI_
  !ENDSITE
  !POLE
  !  !            M2             M1    M4    T2
  !  CHLORITE-AF  MG_MG_MG_MG_   MG_   MG_   SI_SI_ !Al-free-chlorite
  !  CLINOCHLORE  MG_MG_MG_MG_   MG_   AL_   SI_AL_
  !  AMESITE      MG_MG_MG_MG_   AL_   AL_   AL_AL_
  !  DAPHNITE     FE_FE_FE_FE_   FE_   AL_   SI_AL_
  !ENDPOLE
  !
  !HP PARAMETERs X,Y,N=
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
  REAL(dp):: X,Y,N
  !REAL(dp):: G(4)
  !REAL(dp):: W(4,4)
  !INTEGER :: I,J,K
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
  !
  !~ IF(S%NMarg>0) THEN
    !~ W=Zero
    !~ ! retrieve values of W(I,J) from DeCapitani formattted coeffs
    !~ ! this code restricted to binary 12 margules coeffs,
    !~ ! which is general case in HP98 models
    !~ !
    !~ DO I=1,S%NMarg
      !~ M= S%vMarg(I)
      !~ !! W(M%vIPole(1),M%vIPole(2))= M%WG
    !~ ENDDO
    !~ !
    !~ DO I=1,4
      !~ G(I)=Zero
      !~ DO J=1,4
        !~ IF(J/=I) G(I)=G(I) + X(J)*(One-X(I))*W(I,J)
      !~ ENDDO
      !~ DO J=1,4
        !~ DO K=1,4
          !~ IF(J/=I .AND. K/=I) G(I)=G(I) - X(J)*X(K)*W(J,K)
        !~ ENDDO
      !~ ENDDO
      !~ F%vLMarg(I)= G(I)/R_jk/TdgK
      !~ F%vLGam(I)=  F%vLGam(I)+F%vLMarg(I)
    !~ ENDDO
  !~ ENDIF
  !
  RETURN
ENDSUBROUTINE Chlorite_HP98

SUBROUTINE Biotite_HP98( &
& TdgK, &
& X,    &
& vActId,vLGam)
  REAL(dp),INTENT(IN) :: TdgK
  REAL(dp),INTENT(IN) :: X(:)
  REAL(dp),INTENT(OUT):: vActId(:)
  REAL(dp),INTENT(OUT):: vLGam(:)
  !
  !MIXTURE.MODEL BIOTITE_HP98
  !  MODEL SPECIAL
  !  POLE
  !    PHLOGOPITE
  !    ANNITE
  !    EASTONITE
  !    OBIOTITE
  !  END
  !END
  !
  !biotite, KFMASH
  !Fe-Mg ordering -> M1 the Fe-loving site, M2 the other
  !OBIOTITE=PHLOGOPITE[2/3]ANNITE[1/3]
  !SITE
  !  M1 1 MG_FE_AL_
  !  M2 2 MG_FE_
  !  T1 2 SI_AL_
  !ENDSITE
  !POLE
  !   PHLOGOPITE     MG_   MG_MG_  AL_SI_
  !   ANNITE         FE_   FE_FE_  AL_SI_
  !   EASTONITE      AL_   MG_MG_  AL_AL_
  !   OBIOTITE       FE_   MG_MG_  AL_SI_
  !ENDPOLE
  !
  REAL(dp):: X1,X2,Y
  !REAL(dp):: W(4,4)
  !INTEGER :: I,J,K
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
  !! IF(S%NMarg>0) THEN
  !!   W=Zero
  !!   !retrieve values of W(I,J) from DeCapitani formattted coeffs
  !!   !this code restricted to binary 12 margules coeffs,
  !!   !which is general CASE in HP98 models
  !!   DO I=1,S%NMarg
  !!     M=S%vMarg(I)
  !!     ! W(M%vIPole(1),M%vIPole(2))=M%WG
  !!   ENDDO
  !!   !RTLnGam(i)= sum( xj.(1-xi).tDtbMltWij ) - sum( xj.xk .Wjk )
  !!   !            j/=i                   j/=i & k/=i
  !!   DO I=1,4
  !!     G=Zero
  !!     DO J=1,4
  !!       IF(J/=I) G= G + X(J)*(One-X(I))*W(I,J)
  !!     ENDDO
  !!     DO J=1,4
  !!       DO K=1,4
  !!         IF(J/=I .AND. K/=I) G= G - X(J)*X(K)*W(J,K)
  !!       ENDDO
  !!     ENDDO
  !!     F%vLMarg(I)= G /R_jk/TdgK
  !!     F%vLGam(I)=  F%vLGam(I)+F%vLMarg(I)
  !!   ENDDO
  !! ENDIF
  !
  RETURN
ENDSUBROUTINE Biotite_HP98

END MODULE M_MixModel_Special

