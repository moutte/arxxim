MODULE M_Dtb_Write
  
  USE M_Kinds
  USE M_IoTools,ONLY: GetUnit
  USE M_Trace,  ONLY: fTrc,T_,iDebug,fHtm,Warning_
  
  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC:: &
  & DtbAquHkf_Tabulate, &
  & DtbAquThr_Tabulate, &
  & DtbMinHkf_Tabulate, &
  & DtbMin_Tabulate, &
  & DtbSys_Tabulate, &
  & DtbH2OHkf_Tabulate
  !
  PUBLIC:: &
  & DtbAquHkf_Write, &
  & DtbMinHkf_Write, &
  & DtbMinThr_Write
  !
  INTEGER:: fLogK

CONTAINS

SUBROUTINE DtbAquHkf_Tabulate(vTPCond)
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Dtb_Const,  ONLY: R_jk,T_CK,CalToJoule
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_Species,  ONLY: T_Species
  USE M_T_DtbH2OHkf,ONLY: T_DtbH2OHkf,DtbH2OHkf_Calc
  USE M_T_DtbAquHkf,ONLY: T_DtbAquHkf,DtbAquHkf_Calc
  !
  USE M_Dtb_Vars,   ONLY: vDtbAquHkf
  !
  TYPE(T_TPCond),DIMENSION(:),INTENT(IN):: vTPCond
  !
  INTEGER :: N,NN,jTP,iSp,f1,f2,f3 !,f4
  REAL(dp):: T,XX
  REAL(dp),ALLOCATABLE:: tG0rt1(:,:),tG0rt2(:,:),vH0(:)
  REAL(dp),ALLOCATABLE:: Rho(:),Alfa(:),dAlfdT(:),Beta(:)
  !
  TYPE(T_DtbAquHkf):: M
  TYPE(T_Species)  :: S
  TYPE(T_Species),  ALLOCATABLE:: tS(:,:)
  TYPE(T_DtbH2OHkf),ALLOCATABLE:: pW(:)
  !
  REAL(dp):: Rho_ref,Alfa_ref,dAlfdT_ref
  REAL(dp):: T_ref,P_ref,H_ref
  TYPE(T_DtbH2OHkf):: pW_ref
  TYPE(T_Species)  :: S_ref
  TYPE(T_Species),ALLOCATABLE:: vS_ref(:)
  !
  N=  SIZE(vDtbAquHkf)
  NN= SIZE(vTPCond)
  !
  ALLOCATE(Rho(NN))
  ALLOCATE(Alfa(NN))
  ALLOCATE(dAlfdT(NN))
  ALLOCATE(Beta(NN))
  !
  ALLOCATE(tS(1:N,1:NN))
  ALLOCATE(pW(1:N))
  ALLOCATE(vS_ref(1:N))
  ALLOCATE(tG0rt1(1:N,1:NN))
  ALLOCATE(tG0rt2(1:N,1:NN))
  ALLOCATE(vH0(1:N))
  !
  P_ref= 1.0D0
  T_ref= 25.0D0 +T_CK
  !
  !--- compute H0R from G0R and S0_
  DO iSp=1,SIZE(vDtbAquHkf)
    M= vDtbAquHkf(iSp)
    XX= M%G0R    &
    & - T_ref *M%S0Ele /CalToJoule &
    & + T_ref *M%S0_
    WRITE(21,'(A,3G15.6)') M%Name,XX,M%H0R,XX-M%H0R
    ! vDtbAquHkf(iSp)%H0R= XX
    vH0(iSp)= XX *CalToJoule
  END DO
  !---/
  !
  !------------------------ compute properties at standard T,P cond'n --
  CALL DtbH2OHkf_Calc( &
  & T_ref,P_ref,&
  & pW_ref) !out
  Rho_ref=    pW_ref%Rho
  Alfa_ref=   pW_ref%Alfa
  dAlfdT_ref= pW_ref%dAlfdT
  !
  DO iSp=1,N
    CALL DtbAquHkf_Calc( &
    & vDtbAquHkf(iSp), & !in
    & pW_ref,          & !in
    & vS_ref(iSp))
  END DO
  !
  ! print *,'(A,G15.6)',"RHO_REF=    ",Rho_ref
  ! print *,'(A,G15.6)',"Alfa_ref=   ",Alfa_ref
  ! print *,'(A,G15.6)',"dAlfdT_ref= ",dAlfdT_ref
  ! pause
  !------------------------/compute properties at standard T,P cond'n --
  !
  !-------------------- compute Thermodyn. properties for all species --
  DO jTP=1,NN
    !
    CALL DtbH2OHkf_Calc( &
    & vTPCond(jTP)%TdgC+T_CK, & !in
    & vTPCond(jTP)%Pbar, &      !in
    & pW(jTP))
    !
    Rho(jTP)=    pW(jTP)%Rho
    Alfa(jTP)=   pW(jTP)%Alfa
    dAlfdT(jTP)= pW(jTP)%dAlfdT
    Beta(jTP)=   pW(jTP)%Beta
    !
  ENDDO
  !
  DO iSp=1,SIZE(vDtbAquHkf)
    DO jTP=1,NN
      CALL DtbAquHkf_Calc( &
      & vDtbAquHkf(iSp), pW(jTP), & !in
      & tS(iSp,jTP))
    ENDDO
  ENDDO
  !
  !----------------- compute logK according to Anderson Density Model --
  DO iSp=1,SIZE(vDtbAquHkf)
    DO jTP=1,NN
      !
      S= tS(iSp,jTP)
      S_ref= vS_ref(iSp)
      T= vTPCond(jTP)%TdgC+T_CK
      ! Rho= pW(jTP)%Rho
      !
      !--- compute using standard enthalpy tabulated in database
      H_ref= S_ref%H0
      XX= H_ref *(T-T_ref) /T_ref &
      & - S_ref%Cp0 /dAlfdT_ref /T_ref &
      &   *(Alfa_ref *(T-T_ref) + LOG(Rho(jTP)/Rho_ref))
      tG0rt1(iSp,jTP)= vS_ref(iSp)%G0rt - XX /R_jk /T
      !
      !--- compute using enthalpy computed from G0 and S0
      H_ref= vH0(iSp)
      XX= H_ref *(T-T_ref) /T_ref &
      & - S_ref%Cp0 /dAlfdT_ref /T_ref &
      &   *(Alfa_ref *(T-T_ref) + LOG(Rho(jTP)/Rho_ref))
      tG0rt2(iSp,jTP)= vS_ref(iSp)%G0rt - XX /R_jk /T
      !
    ENDDO
  ENDDO
  !-----------------/compute logK according to Anderson Density Model --
  !
  !----------------------------------------------------------/compute --
  !
  !----------------------------------------------------------- OUTPUT --
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirDtbOut)//"hkf_aqu_volumes.tab")
  CALL WriteTPsequence(f1,vTPCond) 
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_aqu_volumes.tab",&
  & "DTB: Hkf/ AquSpecies/ Volumes")
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirDtbOut)//"hkf_aqu_enthalpy.tab")
  CALL WriteTPsequence(f2,vTPCond)
  !-> WRITE logK on fLogK
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_aqu_enthalpy.tab",&
  & "DTB: Hkf/ AquSpecies/ Enthalpy/ KiloJoule/Mole")
  !
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirDtbOut)//"hkf_aqu_gibbs.tab")
  CALL WriteTPsequence(f3,vTPCond) 
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_aqu_gibbs.tab",&
  & "DTB: Hkf/ AquSpecies/ GibbsEnergy")
  !
  !--------------------------------------------- tabulate aqu'species --
  
  WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,"NUM",T_,"NAME",T_,"Rho",T_
  DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') Rho(jTP), T_;   ENDDO
  WRITE(f1,*)
  WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,"NUM",T_,"NAME",T_,"Alfa",T_
  DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') Alfa(jTP), T_;   ENDDO
  WRITE(f1,*)
  WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,"NUM",T_,"NAME",T_,"Beta",T_
  DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') Beta(jTP), T_;   ENDDO
  WRITE(f1,*)
  WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,"NUM",T_,"NAME",T_,"dAlfdT",T_
  DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') dAlfdT(jTP), T_;   ENDDO
  WRITE(f1,*)
  
  DO iSp=1,SIZE(vDtbAquHkf)
  
    M=vDtbAquHkf(iSp)
    
    !!-- write volumes
    ! WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"Vs",T_
    ! DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%Vs, T_;  ENDDO
    ! WRITE(f1,*)
    WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"Vr",T_
    DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%V0, T_;   ENDDO
    WRITE(f1,*)
    WRITE(f1,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"Cp",T_
    DO jTP=1,NN; WRITE(f1,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%Cp0, T_;   ENDDO
    WRITE(f1,*)
    
    !-- enthalpies
    WRITE(f2,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"Hr",T_
    DO jTP=1,NN; WRITE(f2,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%H0/1.0D3, T_;   ENDDO
    WRITE(f2,*)
    !
    !-- Gibbs/RT/ln10
    WRITE(f3,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"logK",T_
    DO jTP=1,NN
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%G0rt/LOG(10.), T_
    ENDDO
    WRITE(f3,*)
    !
    WRITE(f3,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"logK_1",T_
    DO jTP=1,NN
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') tG0rt1(iSp,jTP)/LOG(10.), T_
    ENDDO
    WRITE(f3,*)
    !
    WRITE(f3,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"logK_2",T_
    DO jTP=1,NN
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') tG0rt2(iSp,jTP)/LOG(10.), T_
    ENDDO
    WRITE(f3,*)
    !
    WRITE(f3,'(4(A,A1))', ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,"DELTA",T_
    DO jTP=1,NN
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') (tG0rt2(iSp,jTP)-tG0rt1(iSp,jTP))/LOG(10.), T_
    ENDDO
    WRITE(f3,*)
    !
    !!!  !WRITE SIZE PARAMETER !-> not much variable ... -> not used ...
    !!!  WRITE(f4,'(A,A1,  A15,A1,   A23,A1)',ADVANCE='NO') &
    !!!  &      "AQU",T_,M%Num,T_,M%Name,T_
    !!!  DO jTP=1,NN; WRITE(f4,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%AquSize, T_; ENDDO
    !!!  WRITE(f4,*)
    !
    !!!  !WRITE logK's
    !!!  WRITE(fLogK,'(A,A1,  A15,A1,   A23,A1,      A39,A1)', ADVANCE='NO') &
    !!!  &        "AQU",T_,M%Num,T_,M%Name,T_,M%Formula,T_
    !!!  WRITE(fLogK,'(F7.2,A1)',ADVANCE='NO') tS(iSp,1)%AquSize, T_
    !!!  DO jTP=1,NN; WRITE(fLogK,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%logK, T_; ENDDO
    !!!  WRITE(fLogK,*)
    !
  ENDDO
  !---------------------------------------------/tabulate aqu'species --
  !-----------------------------------------------------------/OUTPUT --
  
  CLOSE(f1)
  CLOSE(f2)
  CLOSE(f3)
  
  DEALLOCATE(Rho)
  DEALLOCATE(Alfa)
  DEALLOCATE(dAlfdT)
  !
  DEALLOCATE(tS)
  DEALLOCATE(vS_ref)
  DEALLOCATE(pW)
  DEALLOCATE(tG0rt1)
  DEALLOCATE(tG0rt2)
  DEALLOCATE(vH0)
  
  PRINT '(A)',"Aqu'Species Volumes in .... "//TRIM(DirDtbOut)//"hkf_aqu_volumes.tab"
  PRINT '(A)',"Aqu'Species Enthalpies in . "//TRIM(DirDtbOut)//"hkf_aqu_enthalpy.tab"
  PRINT '(A)',"Aqu'Species Gibbs in ...... "//TRIM(DirDtbOut)//"hkf_aqu_gibbs.tab"

  RETURN
ENDSUBROUTINE DtbAquHkf_Tabulate

SUBROUTINE DtbMinHkf_Tabulate(vTPCond)
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_Numeric_Const, ONLY: Ln10
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_DtbMinHkf,ONLY: T_DtbMinHkf,DtbMinHkf_Calc
  USE M_T_Species,  ONLY: T_Species
  !
  USE M_Dtb_Vars,   ONLY: vDtbMinHkf
  !
  TYPE(T_TPCond), DIMENSION(:),INTENT(IN):: vTPCond
  !
  TYPE(T_DtbMinHkf):: M
  TYPE(T_Species)  :: S
  INTEGER :: jTP,iSp,f2,f3,N,NN
  !
  TYPE(T_Species),ALLOCATABLE,DIMENSION(:,:)::tS
  !
  N=  SIZE(vDtbMinHkf)
  NN= SIZE(vTPCond)
  ALLOCATE(tS(1:N,1:NN))
  !
  !-------------------- compute Thermodyn. properties for all species --
  DO iSp=1,SIZE(vDtbMinHkf)
    DO jTP=1,NN
      !tS(iSp,jTP)=
      CALL DtbMinHkf_Calc( &
      & vDtbMinHkf(iSp),       &
      & vTPCond(jTP)%TdgC+T_CK,&
      & vTPCond(jTP)%Pbar,     &
      & tS(iSp,jTP))
    ENDDO
  ENDDO
  !---------------------------------------------------------/ compute --
  !
  !----------------------------------------------------------- OUTPUT --
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirDtbOut)//"hkf_min_cp.tab")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_min_cp.tab",&
  & "Cp of Minerals, Hkf model")
  ! 
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirDtbOut)//"hkf_min_gibbs.tab")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_min_gibbs.tab",&
  & "DTB: Gibbs energy of Minerals, Hkf model")
  !
  !
  CALL WriteTPsequence(f2,vTPCond)
  !
  CALL WriteTPsequence(f3,vTPCond)
  !
  DO iSp=1,SIZE(vDtbMinHkf)
    !
    M= vDtbMinHkf(iSp)
    S= tS(iSp,1)
    !
    WRITE(f2,'(3(A,A1),I3,A1)',ADVANCE='NO') "MIN",T_,M%Num,T_,M%Name,T_,M%Div,T_
    WRITE(f3,'(3(A,A1))',ADVANCE='NO') "MIN",T_,M%Num,T_,M%Name,T_
    DO jTP=1,NN
      WRITE(f2,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%Cp0, T_
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%G0rt, T_
    ENDDO
    WRITE(f2,*)
    WRITE(f3,*)
    !
  ENDDO
  !-----------------------------------------------------------/OUTPUT --
  !
  !CLOSE(f1)
  CLOSE(f2)
  CLOSE(f3) 
  !
  DEALLOCATE(tS)
  
  RETURN
ENDSUBROUTINE DtbMinHkf_Tabulate

SUBROUTINE DtbAquThr_Tabulate(vTPCond)
  USE M_Files,        ONLY: DirDtbOut,Files_Index_Write
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Dtb_Const,    ONLY: R_jk,T_CK
  USE M_T_Tpcond,     ONLY: T_TPCond
  USE M_T_Species,    ONLY: T_Species
  USE M_T_DtbAquHkf
  USE M_T_DtbH2OHkf,ONLY: T_DtbH2OHkf,DtbH2OHkf_Calc
  USE M_Fluid_Calc, ONLY: Eos_H2O_Haar
  !
  USE M_Dtb_Vars,   ONLY: vDtbAquHkf
  !
  TYPE(T_TPCond), DIMENSION(:),INTENT(IN):: vTPCond
  !
  INTEGER :: jTP,iSp,f1,f2,f3,N,NN
  REAL(dp):: Pbar,TdgK,G_H2O,V_H2O !for H2O tabulation
  !
  REAL(dp):: Rho_ref,Alfa_ref,dAlfdT_ref
  REAL(dp):: TdgK_ref,Pbar_ref
  TYPE(T_DtbH2OHkf):: pW_ref
  TYPE(T_Species),ALLOCATABLE:: vS_ref(:)
  !
  TYPE(T_DtbAquHkf):: M
  TYPE(T_Species),  ALLOCATABLE:: tS(:,:)
  TYPE(T_DtbH2OHkf),ALLOCATABLE:: pW(:)
  !
  N=  SIZE(vDtbAquHkf)
  NN= SIZE(vTPCond)
  ALLOCATE(tS(1:N,1:NN))
  ALLOCATE(pW(1:NN))
  ALLOCATE(vS_ref(1:N))
  !
  !------------------------ compute properties at standard T,P cond'n --
  Pbar_ref= 1.0D0
  TdgK_ref= 25.0D0 +T_CK
  !
  CALL DtbH2OHkf_Calc( &
  & TdgK_ref,Pbar_ref,&
  & pW_ref) !out
  Rho_ref=    pW_ref%Rho
  Alfa_ref=   pW_ref%Alfa
  dAlfdT_ref= pW_ref%dAlfdT
  !
  DO iSp=1,N
    CALL DtbAquHkf_CalcThr( &
    & vDtbAquHkf(iSp), & !in
    & pW_ref,          & !in
    & vS_ref(iSp))
  END DO
  !------------------------/compute properties at standard T,P cond'n --
  !
  !-------------------- compute Thermodyn. properties for all species --
  DO jTP=1,NN
    CALL DtbH2OHkf_Calc( &
    & vTPCond(jTP)%TdgC+T_CK,vTPCond(jTP)%Pbar,&
    & pW(jTP)) !out
  ENDDO
  DO iSp=1,N
    DO jTP=1,NN
      CALL DtbAquHkf_CalcThr( &
      & vDtbAquHkf(iSp), & !in
      & pW(jTP),         & !in
      & tS(iSp,jTP))
    ENDDO
  ENDDO
  !---------------------------------------------------------- compute --
  !
  !----------------------------------------------------------- OUTPUT --
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirDtbOut)//"thr_aqu_volume.tab")
  CALL WriteTPsequence(f1,vTPCond)
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"thr_aqu_volume.tab",&
  & "DTB: Aqu.Species,HkfData,ThrConvention/ Volumes")
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirDtbOut)//"thr_aqu_enthalpy.tab")
  CALL WriteTPsequence(f2,vTPCond) 
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"thr_aqu_enthalpy.tab",&
  & "DTB: AquSpecies,HkfData,ThrConvention/  Enthalpy/ KiloJoule/Mole")
  !
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirDtbOut)//"thr_aqu_gibbs.tab")
  CALL WriteTPsequence(f3,vTPCond)
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"thr_aqu_gibbs.tab",&
  & "DTB: logK,AquSpecies,HkfData,ThrConvention")
  !
  !------------------------------------------------ tabulate logK_H2O --
  WRITE(fLogK,'(A3,A1,A15,A1,A23,A1,A39,A1,A7,A1)',ADVANCE='NO') &
  & "AQU",T_,"THR_H2O",T_,"H2O",T_,"O(1)H(2)",T_,"   0.00",T_
  DO jTP=1,NN
    Pbar=vTPCond(jTP)%Pbar
    TdgK=vTPCond(jTP)%TdgC+T_CK
    CALL Eos_H2O_Haar(Pbar,TdgK,G_H2O,V_H2O)
    !! G_H2O= G_H2O + 69544.987825D0
    WRITE(fLogK,'(F12.6,A1)',ADVANCE='NO') -G_H2O/R_jk/TdgK/Ln10,T_
  ENDDO
  WRITE(fLogK,*)
  !-----------------------------------------------/ tabulate logK_H2O --
  !
  !WRITE(fLogK,'(A,A1,A,A1,A7,A1)',ADVANCE='NO') &
  !& "H2O",T_,"O(1)H(2)",T_,"   0.00",T_
  !DO jTP=1,NN
  !  Pbar=vTPCond(jTP)%Pbar; TdgK=vTPCond(jTP)%TdgC+273.150D0
  !  CALL Eos_H2O_HolPow91(Pbar,TdgK,G_H2O,V); K_H2O=-G_H2O/R_jk/TdgK/LOG(10.0)
  !  WRITE(fLogK,'(F12.6,A1)',ADVANCE='NO') K_H2O,T_
  !ENDDO
  !WRITE(fLogK,*)
  !
  !--------------------------------------------- tabulate aqu'species --
  DO iSp=1,SIZE(vDtbAquHkf)
    !
    M=vDtbAquHkf(iSp)
    !
    WRITE(f1,'(3(A,A1))',ADVANCE='NO')    "AQU",T_,M%Num,T_,M%Name,T_
    WRITE(f2,'(3(A,A1))',ADVANCE='NO')    "AQU",T_,M%Num,T_,M%Name,T_
    WRITE(f3,'(3(A,A1))',ADVANCE='NO')    "AQU",T_,M%Num,T_,M%Name,T_
    WRITE(fLogK,'(4(A,A1))',ADVANCE='NO') "AQU",T_,M%Num,T_,M%Name,T_,M%Formula,T_
    WRITE(fLogK,'(F7.2,A1)',ADVANCE='NO') tS(iSp,1)%AquSize, T_
    !
    DO jTP=1,NN
      !write volumes
      WRITE(f1,'(F12.6,A1)',ADVANCE='NO') tS(iSp,jTP)%V0, T_
      !write enthalpies
      WRITE(f2,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%H0/1.0D3, T_
      !write gibbs/RT
      WRITE(f3,'(G15.8,A1)',ADVANCE='NO') tS(iSp,jTP)%G0rt, T_
      !write LogK
      WRITE(fLogK,'(G15.8,A1)',ADVANCE='NO') -tS(iSp,jTP)%G0rt/Ln10, T_
    ENDDO
    !
    WRITE(f1,*)
    WRITE(f2,*)
    WRITE(f3,*)
    WRITE(fLogK,*)
    !
  ENDDO !iSp
  !--------------------------------------------/ tabulate aqu'species --
  !
  CLOSE(f1)
  CLOSE(f2)
  CLOSE(f3)
  CLOSE(fLogK)
  !-----------------------------------------------------------/OUTPUT --
  !
  DEALLOCATE(tS)
  DEALLOCATE(pW)
  
  RETURN
ENDSUBROUTINE DtbAquThr_Tabulate

SUBROUTINE DtbMin_Tabulate(sCode,vTPCond) !
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_Numeric_Const, ONLY: Ln10
  USE M_Dtb_Const,  ONLY: T_CK, Tref, R_jK
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_Species,  ONLY: T_Species
  USE M_T_DtbMinHkf,ONLY: T_DtbMinHkf,DtbMinHkf_Calc
  !
  USE M_T_DtbMinThr
  !
  USE M_Dtb_Vars,   ONLY: vDtbMinThr, vDtbMinHkf
  !
  CHARACTER(LEN=*),INTENT(IN):: sCode
  TYPE(T_TPCond),  INTENT(IN):: vTPCond(:)
  !
  TYPE(T_DtbMinThr):: Mt
  TYPE(T_DtbMinHkf):: Mh
  TYPE(T_Species),ALLOCATABLE:: tS(:,:)
  TYPE(T_Species):: S
  INTEGER :: iMin,jTP,nTP,nMin
  INTEGER :: f1,f1b,f2,f3,f4,f5,f
  !! INTEGER :: ff(1:5)
  REAL(dp):: TdgK, Pbar
  CHARACTER(LEN=63):: Str
  !
  ! must have called InitTPsequence beforehand ...
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_density.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_density.tab",&
  & "DTB: minerals,Thr,Density,kg/m3")
  CALL WriteTPsequence(f1,vTPCond)
  !
  CALL GetUnit(f1b)
  OPEN(f1b,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_volum.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_volum.tab",&
  & "DTB: minerals,Thr,Volume,cm3")
  CALL WriteTPsequence(f1b,vTPCond)
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_gibbs.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_gibbs.tab",&
  & "DTB: minerals, Thr, Gibbs, kiloJoules, Benson-Helgeson Conv.")
  CALL WriteTPsequence(f2,vTPCond)  
  !
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_cp.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_cp.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,Cp")
  CALL WriteTPsequence(f3,vTPCond)
  !
  CALL GetUnit(f4)
  OPEN(f4,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_enthalpy.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_enthalpy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,enthalpy")
  CALL WriteTPsequence(f4,vTPCond)
  !
  CALL GetUnit(f5)
  OPEN(f5,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_entropy.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_entropy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,entropy")
  CALL WriteTPsequence(f5,vTPCond)
  !
  nTP=  SIZE(vTPCond)
  IF(sCode=="THR")  nMin= SIZE(vDtbMinThr)
  IF(sCode=="HKF")  nMin= SIZE(vDtbMinHkf)
  ALLOCATE(tS(1:nMin,1:nTP))
  !
  !-------------------- compute Thermodyn. properties for all species --
  DO iMin=1,nMin!____________
    DO jTP=1,nTP
      TdgK= vTPCond(jTP)%TdgC+T_CK
      Pbar= vTPCond(jTP)%Pbar
      !
      IF(sCode=="THR") CALL DtbMinThr_Calc( &
      & vDtbMinThr(iMin), TdgK, Pbar, &
      & S)
      !
      IF(sCode=="HKF") CALL DtbMinHkf_Calc( &
      & vDtbMinHkf(iMin), TdgK, Pbar, &
      & S)
      !
      S%G0rt= S%G0rt*TdgK*R_jK !! -Tref*S%S0Ele
      !
      tS(iMin,jTP)= S
    ENDDO
  ENDDO
  !-------------------/ compute Thermodyn. properties for all species --
  !
  DO iMin=1,nMin
    IF(sCode=="THR") THEN 
      Mt= vDtbMinThr(iMin)
      WRITE(f1, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f1b,'(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f2, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f3, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f4, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f5, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
    ENDIF
    IF(sCode=="HKF") THEN
      Mh= vDtbMinHkf(iMin)
      WRITE(f1, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f1b,'(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f2, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f3, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f4, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f5, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
    ENDIF
    DO jTP=1,nTP
      S= tS(iMin,jTP)
      WRITE(f1, '(G15.6,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_ 
      WRITE(f1b,'(G15.6,A1)',ADVANCE='NO') S%V0*1.E6, T_ 
      WRITE(f2, '(G15.6,A1)',ADVANCE='NO') S%G0rt/1.E3, T_
      WRITE(f3, '(G15.6,A1)',ADVANCE='NO') S%CP0, T_ !/nOx -> "Sdatcaled" to NrOfOxygenSdat
      WRITE(f4, '(G15.6,A1)',ADVANCE='NO') S%H0/1.E3, T_
      WRITE(f5, '(G15.6,A1)',ADVANCE='NO') S%S0/1.E3, T_
    ENDDO
    WRITE(f1, *)
    WRITE(f1b,*)
    WRITE(f2, *)
    WRITE(f3, *)
    WRITE(f4, *)
    WRITE(f5, *)
    !! WRITE(fLogK, *)
  ENDDO
  CLOSE(f1); CLOSE(f1b); CLOSE(f2)
  CLOSE(f3); CLOSE(f4);  CLOSE(f5)
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_all.tab")
  !
  DO iMin=1,nMin
    !
    IF(sCode=="THR") Str= TRIM(vDtbMinThr(iMin)%Name)
    IF(sCode=="HKF") Str= TRIM(vDtbMinHkf(iMin)%Name)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "VOL",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f,'(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%V0*1.E6, T_ 
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "Gibbs",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%G0rt/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "Cp",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%CP0, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "H",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%H0/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "S",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%S0/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
  ENDDO
  !
  CLOSE(f)
  !
  DEALLOCATE(tS)
  
  RETURN
ENDSUBROUTINE DtbMin_Tabulate

SUBROUTINE DtbSys_Tabulate(sCode,vTPCond) !
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_Numeric_Const, ONLY: Ln10
  USE M_Dtb_Const,  ONLY: T_CK, Tref, R_jK
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_Species,  ONLY: T_Species
  USE M_T_DtbMinHkf,ONLY: T_DtbMinHkf,DtbMinHkf_Calc
  !
  USE M_T_DtbMinThr
  !
  USE M_Dtb_Vars,   ONLY: vDtbMinThr, vDtbMinHkf
  !
  CHARACTER(LEN=*),INTENT(IN):: sCode
  TYPE(T_TPCond),  INTENT(IN):: vTPCond(:)
  !
  TYPE(T_DtbMinThr):: Mt
  TYPE(T_DtbMinHkf):: Mh
  TYPE(T_Species),ALLOCATABLE:: tS(:,:)
  TYPE(T_Species):: S
  INTEGER :: iMin,jTP,nTP,nMin
  INTEGER :: f1,f1b,f2,f3,f4,f5,f
  !! INTEGER :: ff(1:5)
  REAL(dp):: TdgK, Pbar
  CHARACTER(LEN=63):: Str
  !
  ! must have called InitTPsequence beforehand ...
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_density.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_density.tab",&
  & "DTB: minerals,Thr,Density,kg/m3")
  CALL WriteTPsequence(f1,vTPCond)
  !
  CALL GetUnit(f1b)
  OPEN(f1b,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_volum.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_volum.tab",&
  & "DTB: minerals,Thr,Volume,cm3")
  CALL WriteTPsequence(f1b,vTPCond)
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_gibbs.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_gibbs.tab",&
  & "DTB: minerals, Thr, Gibbs, kiloJoules, Benson-Helgeson Conv.")
  CALL WriteTPsequence(f2,vTPCond)  
  !
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_cp.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_cp.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,Cp")
  CALL WriteTPsequence(f3,vTPCond)
  !
  CALL GetUnit(f4)
  OPEN(f4,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_enthalpy.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_enthalpy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,enthalpy")
  CALL WriteTPsequence(f4,vTPCond)
  !
  CALL GetUnit(f5)
  OPEN(f5,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_entropy.tab")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//TRIM(sCode)//"_min_entropy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,entropy")
  CALL WriteTPsequence(f5,vTPCond)
  !
  nTP=  SIZE(vTPCond)
  IF(sCode=="THR")  nMin= SIZE(vDtbMinThr)
  IF(sCode=="HKF")  nMin= SIZE(vDtbMinHkf)
  ALLOCATE(tS(1:nMin,1:nTP))
  !
  !-------------------- compute Thermodyn. properties for all species --
  DO iMin=1,nMin!____________
    DO jTP=1,nTP
      TdgK= vTPCond(jTP)%TdgC+T_CK
      Pbar= vTPCond(jTP)%Pbar
      !
      IF(sCode=="THR") CALL DtbMinThr_Calc( &
      & vDtbMinThr(iMin), TdgK, Pbar, &
      & S)
      !
      IF(sCode=="HKF") CALL DtbMinHkf_Calc( &
      & vDtbMinHkf(iMin), TdgK, Pbar, &
      & S)
      !
      S%G0rt= S%G0rt*TdgK*R_jK !! -Tref*S%S0Ele
      !
      tS(iMin,jTP)= S
    ENDDO
  ENDDO
  !-------------------/ compute Thermodyn. properties for all species --
  !
  DO iMin=1,nMin
    IF(sCode=="THR") THEN 
      Mt= vDtbMinThr(iMin)
      WRITE(f1, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f1b,'(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f2, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f3, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f4, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      WRITE(f5, '(3(A,A1),I3,A1)',ADVANCE='NO') Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
    ENDIF
    IF(sCode=="HKF") THEN
      Mh= vDtbMinHkf(iMin)
      WRITE(f1, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f1b,'(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f2, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f3, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f4, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      WRITE(f5, '(3(A,A1),I3,A1)',ADVANCE='NO') Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
    ENDIF
    DO jTP=1,nTP
      S= tS(iMin,jTP)
      WRITE(f1, '(G15.6,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_ 
      WRITE(f1b,'(G15.6,A1)',ADVANCE='NO') S%V0*1.E6, T_ 
      WRITE(f2, '(G15.6,A1)',ADVANCE='NO') S%G0rt/1.E3, T_
      WRITE(f3, '(G15.6,A1)',ADVANCE='NO') S%CP0, T_ !/nOx -> "Sdatcaled" to NrOfOxygenSdat
      WRITE(f4, '(G15.6,A1)',ADVANCE='NO') S%H0/1.E3, T_
      WRITE(f5, '(G15.6,A1)',ADVANCE='NO') S%S0/1.E3, T_
    ENDDO
    WRITE(f1, *)
    WRITE(f1b,*)
    WRITE(f2, *)
    WRITE(f3, *)
    WRITE(f4, *)
    WRITE(f5, *)
    !! WRITE(fLogK, *)
  ENDDO
  CLOSE(f1); CLOSE(f1b); CLOSE(f2)
  CLOSE(f3); CLOSE(f4);  CLOSE(f5)
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirDtbOut)//TRIM(sCode)//"_min_all.tab")
  !
  DO iMin=1,nMin
    !
    IF(sCode=="THR") Str= TRIM(vDtbMinThr(iMin)%Name)
    IF(sCode=="HKF") Str= TRIM(vDtbMinHkf(iMin)%Name)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "VOL",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f,'(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%V0*1.E6, T_ 
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "Gibbs",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%G0rt/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "Cp",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%CP0, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "H",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%H0/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') "S",T_,Trim(Str),T_
    DO jTP=1,nTP
      WRITE(f, '(G15.6,A1)',ADVANCE='NO') tS(iMin,jTP)%S0/1.E3, T_
    ENDDO
    WRITE(f, *)
    !
  ENDDO
  !
  CLOSE(f)
  !
  DEALLOCATE(tS)
  
  RETURN
ENDSUBROUTINE DtbSys_Tabulate

SUBROUTINE DtbH2OHkf_Tabulate(vTPCond)
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_DtbH2OHkf,ONLY: T_DtbH2OHkf,DtbH2OHkf_Calc
  USE M_T_DtbAquHkf
  USE M_IoTools
  !
  TYPE(T_TPCond), DIMENSION(:),INTENT(IN):: vTPCond
  !
  INTEGER::i,f1,NN
  TYPE(T_DtbH2OHkf):: p
  TYPE(T_DtbH2OHkf),ALLOCATABLE:: vDtbH2OHkf(:)
  
  !----------------------------------------- write Solvent properties --
  !
  !------------------------------------- results of HKF-Haar routines --
  NN= SIZE(vTPCond)
  ALLOCATE(vDtbH2OHkf(1:NN))
  
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirDtbOut)//"propsh2o.restab")
  
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"propsh2o.restab",&
  & "DTB: Solvent properties, results of HKF-Haar routines")
  
  WRITE(f1,'(3(A,A1))',ADVANCE='NO')  "TdgC",T_,"Pbar", T_,"RH2O", T_
  WRITE(f1,'(4(A,A1))',ADVANCE='NO')  "G",   T_,"dG_dP",T_,"dG_dT",T_,"d2G_dT2",T_
  !WRITE(f1,'(4(A,A1))',ADVANCE='NO') "Z",   T_,"Q",    T_,"Y",    T_,"X",T_
  WRITE(f1,'(2(A,A1))')               "dhA", T_,"dhB"
  
  DO I=1,NN
    CALL DtbH2OHkf_Calc( &
    & vTPCond(i)%TdgC+T_CK,vTPCond(i)%Pbar, &
    & vDtbH2OHkf(i))
  ENDDO
  
  DO I=1,NN
    p= vDtbH2OHkf(i)
    WRITE(f1,'(F7.3,A1,F7.3,A1,G15.8,A1)',ADVANCE='NO') &
    & p%TdgK-T_CK,T_, &
    & p%Pbar,     T_, &
    & p%Rho,      T_
    WRITE(f1,'(4(G15.8,A1))',ADVANCE='NO') &
    & p%Gshok,T_,p%dGShdP,T_,p%dGShdT,T_,p%d2GShdT2,T_
    WRITE(f1,*)
  ENDDO !I
  
  CLOSE(f1)
  !----------------------------------------/ write Solvent properties --
  !
  WRITE(fLogK,'(/,A)') "TP.TABLE"
  !
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%TdgK-T_CK,Opt_S="TdgC")
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Pbar,     Opt_S="Pbar")
  !
  WRITE(fLogK,'(A,/)') "END TP.TABLE"
  !
  WRITE(fLogK,'(/,A)') "SOLVENT"
  !
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Rho, Opt_S="Rho")
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Eps, Opt_S="Eps")
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%dhA, Opt_S="DHA")
  CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%dhB, Opt_S="DHB")
  !! CALL OutStrVec(fLogK,vDtbH2OHkf(1:NN)%bDOt,Opt_S="BDOT")
  !
  WRITE(fLogK,'(A,/)') "END"
  !
  !CALL OutStrVec(f1,S="Gshok__",vDtbH2OHkf(1:NN)%Gshok)
  !CALL OutStrVec(f1,S="dGShdP_",vDtbH2OHkf(1:NN)%dGShdP)
  !CALL OutStrVec(f1,S="dGShdT_",vDtbH2OHkf(1:NN)%dGShdT)
  !CALL OutStrVec(f1,S="d2GSdT2",vDtbH2OHkf(1:NN)%d2GShdT2)
  !CALL OutStrVec(f1,S="dEpsdP_",vDtbH2OHkf(1:NN)%dEpsdP)
  !CALL OutStrVec(f1,S="dEpsdT_",vDtbH2OHkf(1:NN)%dEpsdT)
  !CALL OutStrVec(f1,S="d2EpsdT",vDtbH2OHkf(1:NN)%d2EpsdT2)
  !CALL OutStrVec(f1,S="Qeps___",vDtbH2OHkf(1:NN)%Qeps)
  !CALL OutStrVec(f1,S="Xeps___",vDtbH2OHkf(1:NN)%Xeps)
  !CALL OutStrVec(f1,S="Yeps___",vDtbH2OHkf(1:NN)%Yeps)
  !CALL OutStrVec(f1,S="Zeps___",vDtbH2OHkf(1:NN)%Zeps)
  !
ENDSUBROUTINE DtbH2OHkf_Tabulate

SUBROUTINE WriteTPsequence(f,vTPCond)
  
  USE M_T_Tpcond,ONLY: T_TPCond
  !
  TYPE(T_TPCond), DIMENSION(:),INTENT(IN):: vTPCond
  !
  INTEGER:: f
  INTEGER:: i,NN
  
  NN= SIZE(vTPCond)
  
  WRITE(f,'(4(A6,A1))',ADVANCE='NO') &
  & ".TYP",T_,".SOURCE",T_,".SPECIES",T_,".TDGC",T_
  DO I=1,NN; WRITE(f,'(G15.4, A1)',ADVANCE='NO') vTPCond(I)%TdgC,T_; ENDDO; WRITE(f,*)
  
  WRITE(f,'(4(A6,A1))',ADVANCE='NO') &
  & ".TYP",T_,".SOURCE",T_,".SPECIES",T_,".PBAR",T_
  DO I=1,NN; WRITE(f,'(G15.4, A1)',ADVANCE='NO') vTPCond(I)%Pbar,T_; ENDDO; WRITE(f,*)
  
  RETURN
ENDSUBROUTINE WriteTPsequence

SUBROUTINE DtbMinHkf_Write
!--
!-- check result DtbMinHKF_Read
!--
  USE M_Files,    ONLY: DirDtbOut,Files_Index_Write
  USE M_Dtb_Const,ONLY: Tref
  USE M_T_DtbMinHkf
  USE M_Dtb_Vars, ONLY: vDtbMinHkf
  !
  TYPE(T_DtbMinHkf):: M
  INTEGER       :: I,f
  REAL(dp)      :: X
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbMinHkf_Write"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirDtbOut)//"hkf_min.restab")
  !DO I=1,nEleDtb; WRITE(f,'(G15.8,A1)',ADVANCE='NO') vS0Ele(I),T_; ENDDO; WRITE(f,*)
  
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"hkf_min.restab",&
  & "DTB: table of data for Min.species, HKF")
  
  WRITE(f,'(3(A,A1))',ADVANCE='NO') &
  & "Num",     T_,&
  & "Name",    T_,&
  & "Formula", T_
  WRITE(f,'(6(A,A1))',ADVANCE='NO') &
  & "(G0R-H0R)/Tr+S0",T_, &
  & "S0Ele",          T_, &
  & "G0R",            T_, &
  & "H0R",            T_, &
  & "S0",             T_, &
  & "V0R",            T_
  WRITE(f,'(4(A,A1))') &
  & "MK1_1__",T_,&
  & "MK1_2__",T_,&
  & "MK1_3__",T_,&
  & "MK1_4__",T_
  
  DO I=1,SIZE(vDtbMinHkf)
  
    M= vDtbMinHkf(I)
    X= (M%G0R -M%H0R)/Tref +M%S0_ 
    
    WRITE(f,'(3(A,A1))',ADVANCE='NO') &
    & M%Num,    T_, &
    & M%Name,   T_, &
    & M%Formula,T_
    WRITE(f,'(6(G15.8,A1))',ADVANCE='NO') &
    & X,       T_,  &
    & M%S0Ele, T_,  &
    & M%G0R,   T_,  &
    & M%H0R,   T_,  &
    & M%S0_,   T_,  &
    & M%V0R,   T_
    WRITE(f,'(3(G15.8,A1))') &
    & M%MK1(1),T_, &
    & M%MK1(2),T_, &
    & M%MK1(3),T_
  
  ENDDO
  
  CLOSE(f)
  
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbMinHkf_Write"
ENDSUBROUTINE DtbMinHkf_Write

SUBROUTINE DtbAquHkf_Write
!--
!-- check result DtbAquHKF_Read
!-- make a new version of Hkf base, one species per line
!--

  USE M_Dtb_Const,ONLY: Tref
  USE M_Files,    ONLY: DirDtbOut,Files_Index_Write
  USE M_T_DtbAquHkf
  USE M_Dtb_Vars, ONLY: vDtbAquHkf
  !
  TYPE(T_DtbAquHkf):: M
  INTEGER:: I,f
  
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbAquHkf_Write"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirDtbOut)//"_aqu_hkf.restab")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"aqu_hkf.restab",&
  & "DTB: table of DATA for Aqu.species, HKF")
  !
  WRITE(f,'(14(A,A1))') &
  & "Num",T_,"Name",T_,"Formula",T_, &
  !!& "(G0R-H0R)/Tr+S0",T_, "S0Ele",T_,&
  & "GfAq",T_, "HfAq",T_, "SPrTrAq",T_, &
  & "A1",T_,"A2",T_,"A3",T_,     "A4",T_, &
  & "C1",T_,"C2",T_,"Omega__",T_,"Chg",T_
  !
  DO I=1,SIZE(vDtbAquHkf)
    M=vDtbAquHkf(I)
    WRITE(f,'(3(A,A1))',ADVANCE='NO') &
    & M%Num,T_,&
    & M%Name,T_,&
    & M%Formula,T_
    !!Y=DOT_PRODUCT(M%Stoik(1:nEl),vEle(1:nEl)%S0)
    !!WRITE(f,'(5(G15.8,A1))',ADVANCE='NO') X,   T_,Y,   T_,M%G0R, T_,M%H0R,T_,M%S0_,T_
    WRITE(f,'(3(G15.8,A1))',ADVANCE='NO') M%G0R, T_,M%H0R,T_,M%S0_,T_
    WRITE(f,'(5(G15.8,A1))',ADVANCE='NO') M%A1,T_,M%A2,T_,M%A3,  T_,M%A4, T_
    WRITE(f,'(3(G15.8,A1),I3,A1)')        M%C1,T_,M%C2,T_,M%wref,T_,M%Chg,T_
  ENDDO
  !
  CLOSE(f)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbAquHkf_Write"
  
  RETURN
ENDSUBROUTINE DtbAquHkf_Write

SUBROUTINE DtbMinThr_Write
!--
!--check result DtbMinHKF_Read
!--
  
  USE M_Files,      ONLY: DirDtbOut,Files_Index_Write
  USE M_Dtb_Vars,   ONLY: vDtbMinThr
  USE M_T_DtbMinThr,ONLY: T_DtbMinThr
  
  TYPE(T_DtbMinThr):: M
  INTEGER          :: I,f
  
  !!=PARAMETERs=========================================================
  !  INTEGER :: nLanda !number of landau transitions
  !  INTEGER :: codVol !code for volume FUNCTION
  !  REAL(dp)::& !for landau transitions
  !  & TQ1B(1:3),TRE(1:3), ASPK(1:3),BSPK(1:3),DHTR(1:3),&
  !  & TEQ(1:3), DVTR(1:3),DVDT(1:3),DVDP(1:3)
  !  REAL(dp):: G0R,H0R,S0_,V0R
  !  REAL(dp):: S0Ele,WeitKg
  !  REAL(dp):: AA0,AAT,BB0,BBT
  !  !          !for VDW/RDK/... cubic EoS of Gases (VanDerWaals,RedlichKwong,...)
  !  REAL(dp):: Tcrit,Pcrit,Acentric !for SRK/PRS/... cubic EoS of gases
  !  REAL(dp):: K1,K2,K3,K4,K5,K6,K7,K8,K9  !coeff's for Cp formulas
  !  REAL(dp):: D1,D2,D3,D4,D5,D6,D7,D8,D9  !disorder /Berman (HP:D1,D2,D3)
  !  REAL(dp):: TD0,TDMax,VAdj              !disorder /Berman 
  !  REAL(dp):: VTA,VTB,VPA,VPB             !volume FUNCTION
  !  !__________________________a la Berman: USE VTA,VTB,VPA,VPB
  !  !__________________________a la HP:     USE VTA,VTB,VPA,TKRI,SMA
  !  !VL0,VLA,VLN,VL2,VAA,VAB,VB, & !older formats ??
  !  REAL(dp):: TKri,SMA
  !  LOGICAL::Dis !,Min,VdW,RdK !,TL1
  !!--------------------------------------------------------------------
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbMinThr_Write"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(DirDtbOut)//"_min_thr.restab")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirDtbOut)//"_min_thr.restab",&
  & "DTB: DATA for Min.species, Theriak, in spREADsheet form")
  !
  WRITE(f,'(2(A,A1))',ADVANCE='NO') &
  & "Name                ",T_,& !LEN=20
  & "Formula                                ",T_
  WRITE(f,'(4(A,A1))',ADVANCE='NO') &
  & "G0R(J)         ",T_,&
  & "H0R(J)         ",T_,&
  & "S0(J/K)        ",T_,&
  & "V0R(J/bar)     ",T_
  WRITE(f,'(9(A,A1))',ADVANCE='NO') &
  & "K1             ",T_,&
  & "K4             ",T_,&
  & "K3             ",T_,&
  & "K8             ",T_,&
  & "K6             ",T_,&
  & "K2             ",T_,&
  & "K5             ",T_,&
  & "K7             ",T_,&
  & "K9             ",T_
  WRITE(f,'(6(A,A1))',ADVANCE='NO') &
  & "VTA*1E5         ",T_,&
  & "VTB*1E5         ",T_,&
  & "VPA*1E5         ",T_,&
  & "VPB*1E8         ",T_,&
  & "Tkri           ",T_,&
  & "SMA            ",T_
  WRITE(f,'(9(A,A1))',ADVANCE='NO') &
  & "D1             ",T_,&
  & "D4             ",T_,&
  & "D3             ",T_,&
  & "D8             ",T_,&
  & "D6             ",T_,&
  & "D2             ",T_,&
  & "D5             ",T_,&
  & "D7             ",T_,&
  & "D9             ",T_
  WRITE(f,*)
  !
  DO I=1,SIZE(vDtbMinThr)
    M=vDtbMinThr(I)
    !
    WRITE(f,'(2(A,A1))',ADVANCE='NO') M%Name,T_,M%Formula,T_
    !
    WRITE(f,'(4(G15.8,A1))',ADVANCE='NO')  &
    & M%G0R,T_,M%H0R,T_,M%S0_,T_,M%V0R,T_
    !
    WRITE(f,'(9(G15.8,A1))',ADVANCE='NO')  &
    & M%K1, T_,M%K4, T_,M%K3, T_,M%K8, T_, &
    & M%K6, T_,M%K2, T_,M%K5, T_,M%K7, T_,M%K9, T_
    !
    WRITE(f,'(6(G15.8,A1))',ADVANCE='NO')  &
    & M%VTA*1.0D5,T_,M%VTB*1.0D5,T_, &
    & M%VPA*1.0D5,T_,M%VPB*1.0D8,T_, &
    & M%Tkri,T_,     M%SMA,T_
    !
    WRITE(f,'(9(G15.8,A1))',ADVANCE='NO')  &
    & M%D1,T_,M%D2,T_,M%D3,T_,M%D4,T_, &
    & M%D5,T_,M%D6,T_,M%D7,T_,M%D8,T_,M%D9,T_ 
    !
    WRITE(f,*)
  ENDDO
  !
  CLOSE(f)
  !
  CALL Warning_("Data written to "//TRIM(DirDtbOut)//"_min_thr.restab")
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbMinThr_Write"
ENDSUBROUTINE DtbMinThr_Write

ENDMODULE M_Dtb_Write
