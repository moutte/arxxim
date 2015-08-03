MODULE M_Dtb_Calc
!--
!-- routines to compute species properties from databases
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,Stop_,T_,Pause_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Species_TP_Update_fromDtb
  PUBLIC:: SpeciesMin_TP_Update_fromDtb
  PUBLIC:: Dtb_TP_Check
  PUBLIC:: DtbSpc_GrtTable_Build
  PUBLIC:: SpeciesDtb_ToSpecies
  !
CONTAINS

SUBROUTINE SpeciesDtb_ToSpecies( &
& vSpcDtb, & ! IN
& vSpc)      ! OUT

  USE M_T_Species,ONLY: T_Species,T_SpeciesDtb,Species_Index
  !
  USE M_Dtb_Vars, ONLY: &
  & vDtbMinHkf,vDtbMinThr,vDtbAquHkf,vDtbLogKtbl,vDtblogKanl
  !
  TYPE(T_SpeciesDtb),INTENT(IN)   :: vSpcDtb(:)
  TYPE(T_Species),   INTENT(INOUT):: vSpc(:)
  !
  TYPE(T_Species):: S
  INTEGER:: I,J,N
  !
  N= SIZE(vSpcDtb)
  !
  DO J= 1,N

    I= vSpcDtb(J)%Indx

    SELECT CASE(vSpcDtb(J)%DtbModel)

    CASE("H2O_HGK")  ! HGK= Haar-Gallagher-Kell
      S%NamSp=    "H2O"
      S%WeitKg=   0.0180152D0
      S%Formula=  "H(2)O(1)"
      S%Typ=      "AQU"

    CASE("AQU_HKF")
      S%NamSp=    TRIM(vDtbAquHkf(I)%Name)
      S%Formula=  TRIM(vDtbAquHkf(I)%Formula)
      S%WeitKg=   vDtbAquHkf(I)%WeitKg
      !S%AquSize=  vDtbAquHkf(I)%AquSize
      S%Typ=      "AQU"

    !! CASE("AQU_THR")

    CASE("MIN_HKF","GAS_HKF")
      S%NamSp=    TRIM(vDtbMinHkf(I)%Name)
      S%Formula=  TRIM(vDtbMinHkf(I)%Formula)
      S%WeitKg=   vDtbMinHkf(I)%WeitKg
      S%Typ=      vDtbMinHkf(I)%Typ

    CASE("MIN_THR","GAS_THR")
      S%NamSp=    TRIM(vDtbMinThr(I)%Name)
      S%Formula=  TRIM(vDtbMinThr(I)%Formula)
      S%WeitKg=   vDtbMinThr(I)%WeitKg
      S%Typ=      vDtbMinThr(I)%Typ

    CASE("LOGKTBL")
      S%NamSp=    TRIM(vDtbLogKtbl(I)%Name)
      S%Formula=  TRIM(vDtbLogKtbl(I)%Formula)
      S%WeitKg=   vDtbLogKtbl(I)%WeitKg
      S%AquSize=  vDtbLogKtbl(I)%AquSize
      S%Typ=      vDtbLogKtbl(I)%Typ

    CASE("LOGKANL")
      S%NamSp=    TRIM(vDtbLogKanl(I)%Name)
      S%Formula=  TRIM(vDtblogKanl(I)%Formula)
      S%WeitKg=   vDtbLogKanl(I)%WeitKg
      S%AquSize=  vDtbLogKanl(I)%AquSize
      S%Typ=      vDtbLogKanl(I)%Typ

    ENDSELECT

    S%iDiscret= 0
    S%iDtb=     J
    vSpc(J)=    S

    ! print *,S%NamSp
    
  ENDDO

  ! pause_
  
  !to prevent some bugs,
  !put solvent species H2O as species 1 in vSpc
  N= Species_Index("H2O",vSpc)
  IF(N/=0 .AND. N/=1) THEN
    S=       vSpc(1)
    vSpc(1)= vSpc(N)
    vSpc(N)= S
  ENDIF

  IF(iDebug<4) RETURN

  PRINT *,"SpeciesDtb_ToSpecies"
  DO i=1,N
    PRINT *,vSpc(i)%NamSp,vSpc(i)%Typ,"=",vSpcDtb(vSpc(i)%iDtb)%DtbModel
  ENDDO

  RETURN
ENDSUBROUTINE SpeciesDtb_ToSpecies

SUBROUTINE Species_TP_Update_fromDtb( &
& TdgK,Pbar,PropsH2O,vSpcDtb, &
& S)
!--
!-- compute thermodyn.prop's of species S at TdgK,Pbar,PropsH2O
!-- according to S%Model and S%iDtb
!-- should be called only when S%iDtb/-0
!--
  USE M_T_Species,  ONLY: T_Species,T_SpeciesDtb
  !
  USE M_T_DtbH2OHkf, ONLY: T_DtbH2OHkf
  USE M_T_DtbAquHkf, ONLY: DtbAquHkf_Calc,DtbAquHkf_CalcThr
  USE M_T_DtbMinHkf, ONLY: DtbMinHkf_Calc
  USE M_T_DtbMinThr, ONLY: DtbMinThr_Calc
  USE M_T_DtbLogKtbl,ONLY: DtbLogKtbl_Calc
  USE M_T_DtbLogKanl,ONLY: DtbLogKanl_Calc
  !
  USE M_Dtb_Vars,   ONLY: vDtbAquHkf,vDtbMinHkf,vDtbMinThr  ! data-
  USE M_Dtb_Vars,   ONLY: vDtbLogKtbl,vDtbLogKanl           ! bases
  !
  USE M_Fluid_Calc, ONLY: Eos_H2O_Haar_Ghiorso
  !!USE M_Fluid_Calc, ONLY: CalcGH2O_Supcrt !,Eos_H2O_Haar
  
  !---------------------------------------------------------------------
  REAL(dp),          INTENT(IN) :: TdgK,Pbar
  TYPE(T_DtbH2OHkf), INTENT(IN) :: PropsH2O
  TYPE(T_SpeciesDtb),INTENT(IN) :: vSpcDtb(:)
  !
  TYPE(T_Species),   INTENT(INOUT):: S
  !---------------------------------------------------------------------
  INTEGER:: i,j
  !---------------------------------------------------------------------
  
  i= S%iDtb
  IF(i==0) RETURN
  !
  j= vSpcDtb(i)%Indx
  !
  SELECT CASE(TRIM(vSpcDtb(i)%DtbModel))

  CASE("H2O_HGK"); CALL Eos_H2O_Haar_Ghiorso(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
  ! CASE("H2O_HGK") ; CALL CalcGH2O_Supcrt(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)

  CASE("AQU_HKF") ; CALL DtbAquHkf_Calc(vDtbAquHkf(j),PropsH2O,S)
  CASE("AQU_THR") ; CALL DtbAquHkf_CalcThr(vDtbAquHkf(j),PropsH2O,S)

  CASE("MIN_HKF") ; CALL DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)
  CASE("GAS_HKF") ; CALL DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)

  CASE("MIN_THR") ; CALL DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)
  CASE("GAS_THR") ; CALL DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)

  CASE("LOGKTBL") ; CALL DtbLogKtbl_Calc(vDtbLogKtbl(j),TdgK,Pbar,S)
  CASE("LOGKANL") ; CALL DtbLogKanl_Calc(vDtbLogKanl(j),TdgK,Pbar,S)

  END SELECT
  
  RETURN
ENDSUBROUTINE Species_TP_Update_fromDtb

SUBROUTINE SpeciesMin_TP_Update_fromDtb(TdgK,Pbar,vSpcDtb, S)
!--
!-- calculate thermodyn.properties of species S at TdgK,Pbar
!-- same as DtbSpc_TP_Update, for non'aqu'species ONLY -> no need for PropsH2O
!-- should be called ONLY when S%iDtb/-0
!--
  USE M_T_Species,  ONLY: T_Species,T_SpeciesDtb
  !
  USE M_T_DtbH2OHkf, ONLY: T_DtbH2OHkf
  USE M_T_DtbMinHkf, ONLY: DtbMinHkf_Calc
  USE M_T_DtbMinThr, ONLY: DtbMinThr_Calc
  USE M_T_DtbLogKtbl,ONLY: DtbLogKtbl_Calc
  USE M_T_DtbLogKanl,ONLY: DtbLogKanl_Calc
  !
  USE M_Dtb_Vars,   ONLY: &
  & vDtbAquHkf,vDtbMinHkf,vDtbMinThr,vDtbLogKtbl, vDtbLogKanl !-> the databases
  !
  USE M_Fluid_Calc, ONLY: Eos_H2O_Haar_Ghiorso
  !!USE M_Fluid_Calc, ONLY: CalcGH2O_Supcrt !,Eos_H2O_Haar
  !---------------------------------------------------------------------
  REAL(dp),          INTENT(IN) :: TdgK,Pbar
  TYPE(T_SpeciesDtb),INTENT(IN) :: vSpcDtb(:)
  !
  TYPE(T_Species),   INTENT(INOUT):: S
  !---------------------------------------------------------------------
  INTEGER:: i,j
  !---------------------------------------------------------------------
  i= S%iDtb
  !
  IF(i==0) RETURN
  !
  j= vSpcDtb(i)%Indx
  !
  SELECT CASE(vSpcDtb(i)%DtbModel)
    CASE("H2O_HGK"); CALL Eos_H2O_Haar_Ghiorso(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
    ! CASE("H2O_HGK"); CALL CalcGH2O_Supcrt(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
    CASE("MIN_HKF"); CALL DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)
    CASE("MIN_THR"); CALL DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)
    CASE("GAS_HKF"); CALL DtbMinHkf_Calc(vDtbMinHkf(j),TdgK,Pbar,S)
    CASE("GAS_THR"); CALL DtbMinThr_Calc(vDtbMinThr(j),TdgK,Pbar,S)
    CASE("LOGKTBL"); CALL DtbLogKtbl_Calc(vDtbLogKtbl(j),TdgK,Pbar,S)
    CASE("LOGKANL"); CALL DtbLogKanl_Calc(vDtbLogKanl(j),TdgK,Pbar,S)
  END SELECT
  !
  RETURN
ENDSUBROUTINE SpeciesMin_TP_Update_fromDtb

SUBROUTINE Dtb_TP_Check( &
& DtbFormat,DtbLogK_vTPCond,Psat_Auto, &
& TdgK,Pbar,Ok,Msg)
!--
!-- check that TdgK is in the validity range of the database,
!-- and (in case of logK database) compute corresponding pressure
!--
  USE M_T_Tpcond,  ONLY: T_TPCond
  USE M_Dtb_Const, ONLY: T_CK, PminHSV, PmaxHSV, TCminHSV, TCmaxHSV, Pref
  USE M_Fluid_Calc,ONLY: Eos_H2O_Psat
  !
  CHARACTER(LEN=*),INTENT(IN)   :: DtbFormat
  TYPE(T_TPCond),  INTENT(IN)   :: DtbLogK_vTPCond(:)
  LOGICAL,         INTENT(IN)   :: Psat_Auto
  REAL(dp),        INTENT(INOUT):: TdgK,Pbar
  LOGICAL,         INTENT(OUT)  :: Ok
  CHARACTER(LEN=*),INTENT(OUT)  :: Msg
  !
  REAL(dp):: Tmin,Tmax,Psat !,Pmin,Pmax
  !
  Ok=  .true.
  Msg= "OK"
  !
  !------------------------------------------------ check Temperature --
  SELECT CASE(DtbFormat)

  CASE("LOGKTBL")
    Tmin= DtbLogK_vTPCond(1)%TdgC +T_CK
    Tmax= DtbLogK_vTPCond(SIZE(DtbLogK_vTPCond))%TdgC +T_CK

  CASE DEFAULT
    Tmin= TCminHSV + T_CK
    Tmax= TCmaxHSV + T_CK

  END SELECT
  !
  IF(TdgK <Tmin .or. TdgK>Tmax) THEN
    Ok= .false.
    Msg= "Temperature outside validity range of the database"
    RETURN !==================================================< error ==
  ENDIF
  !-----------------------------------------------/ check Temperature --
  !
  !--------------------------------------------------- check Pressure --
  SELECT CASE(DtbFormat)
  !--------------------------------------------------- logK databases --
  CASE("LOGKTBL","LOGKANL")
    IF(TdgK <= 100.0D0+T_CK) THEN
      IF(Pbar<Pref) THEN
        IF(Psat_Auto) THEN
          Pbar= Pref
        ELSE
          Ok= .false.
          Msg= "Pressure < Psat(T) = not valid for HKF model"
          RETURN !============================================= error ==
        ENDIF
      ENDIF
    ELSE
      CALL Eos_H2O_Psat(TdgK,Psat)
      IF(Pbar<Psat) THEN
        IF(Psat_Auto) THEN
          Pbar= Psat
        ELSE
          Ok= .false.
          Msg= "Pressure < Psat(T) = not valid for HKF model"
          RETURN
        ENDIF
      ENDIF
    ENDIF
  !-------------------------------------------------- other databases --
  CASE DEFAULT
    IF(Pbar <PminHSV .or. Pbar>PmaxHSV) THEN
      Ok= .false.
      Msg= "Pressure outside validity range of the database"
      RETURN !================================================= error ==
    ENDIF
    !
    !------------------------------------- if P-PsatH2O(T) then error --
    CALL Eos_H2O_Psat(TdgK,Psat)
    !! Pbar= max(Pbar,Psat)
    IF(Pbar < Psat) THEN
      IF(Psat_Auto) THEN
        Pbar= Psat
      ELSE
        Ok= .false.
        Msg= "Pressure < Psat(T) = not valid for HKF model"
        RETURN !=============================================== error ==
      ENDIF
    ENDIF
    !
  END SELECT
  !--------------------------------------------------/ check Pressure --
  !
ENDSUBROUTINE Dtb_TP_Check

SUBROUTINE DtbSpc_GrtTable_Build( &
& vTPCond,vSpcDtb,vSpc, &
& tGrt)

  USE M_T_Species,ONLY: T_Species,T_SpeciesDtb
  USE M_T_Tpcond, ONLY: T_TPCond
  USE M_Dtb_Const,ONLY: T_CK
  USE M_T_DtbH2OHkf,ONLY: DtbH2OHkf_Calc,T_DtbH2OHkf
  
  TYPE(T_TPCond),    INTENT(IN) :: vTPCond(:)
  TYPE(T_SpeciesDtb),INTENT(IN) :: vSpcDtb(:)
  TYPE(T_Species),   INTENT(IN) :: vSpc(:)
  REAL(dp),          INTENT(OUT):: tGrt(:,:)
  
  INTEGER :: iTP,jSp
  REAL(dp):: TdgK,Pbar
  TYPE(T_DtbH2OHkf)::PropsH2O
  TYPE(T_Species):: S
  
  DO iTP=1,SIZE(vTPCond)
    !
    TdgK= vTPcond(iTP)%TdgC + T_CK
    Pbar= vTPcond(iTP)%Pbar
    !--- solvent properties, for aqu'species
    CALL DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
    !---
    !
    DO jSp=1,SIZE(vSpc)
      !
      S= vSpc(jSp)
      !
      IF(S%iDtb>0) THEN
        CALL Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,S)
        !-> update vSpc(jSp)%G0rt !-> =G/RT
        tGrt(jSp,iTP)= S%G0rt
      !ELSEIF(vSpc(jSp)%iDtb>0) THEN
      ENDIF
      !
    ENDDO
    !
  ENDDO
  
ENDSUBROUTINE DtbSpc_GrtTable_Build

ENDMODULE M_Dtb_Calc

