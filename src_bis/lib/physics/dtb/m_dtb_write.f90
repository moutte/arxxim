module M_Dtb_Write
  
  use M_Kinds
  use M_IoTools,only: GetUnit
  use M_Trace,  only: fTrc,T_,iDebug,fHtm,Warning_
  
  implicit none

  private
  !
  public:: &
  & DtbAquHkf_Tabulate, &
  & DtbAquThr_Tabulate, &
  & DtbMinHkf_Tabulate, &
  & DtbMin_Tabulate, &
  & DtbSys_Tabulate, &
  & DtbH2OHkf_Tabulate
  !
  public:: &
  & DtbAquHkf_Write, &
  & DtbMinHkf_Write, &
  & DtbMinThr_Write
  !
  integer:: fLogK

contains

subroutine DtbAquHkf_Tabulate(vTPCond)
  use M_Numeric_Const,only: Ln10
  use M_Dtb_Const,  only: R_jk,T_CK,CalToJoule
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_T_Tpcond,   only: T_TPCond
  use M_T_Species,  only: T_Species
  use M_T_DtbH2OHkf,only: T_DtbH2OHkf,DtbH2OHkf_Calc
  use M_T_DtbAquHkf,only: T_DtbAquHkf,DtbAquHkf_Calc
  !
  use M_Dtb_Vars,   only: vDtbAquHkf
  !
  type(T_TPCond),dimension(:),intent(in):: vTPCond
  !
  integer :: N,NN,jTP,iSp,f1,f2,f3 !,f4
  real(dp):: T,XX
  real(dp),allocatable:: tG0rt1(:,:),tG0rt2(:,:),vH0(:)
  real(dp),allocatable:: Rho(:),Alfa(:),dAlfdT(:),Beta(:)
  !
  type(T_DtbAquHkf):: M
  type(T_Species)  :: S
  type(T_Species),  allocatable:: tS(:,:)
  type(T_DtbH2OHkf),allocatable:: pW(:)
  !
  real(dp):: Rho_ref,Alfa_ref,dAlfdT_ref
  real(dp):: T_ref,P_ref,H_ref
  type(T_DtbH2OHkf):: pW_ref
  type(T_Species)  :: S_ref
  type(T_Species),allocatable:: vS_ref(:)
  !
  N=  size(vDtbAquHkf)
  NN= size(vTPCond)
  !
  allocate(Rho(NN))
  allocate(Alfa(NN))
  allocate(dAlfdT(NN))
  allocate(Beta(NN))
  !
  allocate(tS(1:N,1:NN))
  allocate(pW(1:N))
  allocate(vS_ref(1:N))
  allocate(tG0rt1(1:N,1:NN))
  allocate(tG0rt2(1:N,1:NN))
  allocate(vH0(1:N))
  !
  P_ref= 1.0D0
  T_ref= 25.0D0 +T_CK
  !
  !--- compute H0R from G0R and S0_
  do iSp=1,size(vDtbAquHkf)
    M= vDtbAquHkf(iSp)
    XX= M%G0R    &
    & - T_ref *M%S0Ele /CalToJoule &
    & + T_ref *M%S0_
    write(21,'(A,3G15.6)') M%Name,XX,M%H0R,XX-M%H0R
    ! vDtbAquHkf(iSp)%H0R= XX
    vH0(iSp)= XX *CalToJoule
  end do
  !---/
  !
  !------------------------ compute properties at standard T,P cond'n --
  call DtbH2OHkf_Calc( &
  & T_ref,P_ref,&
  & pW_ref) !out
  Rho_ref=    pW_ref%Rho
  Alfa_ref=   pW_ref%Alfa
  dAlfdT_ref= pW_ref%dAlfdT
  !
  do iSp=1,N
    call DtbAquHkf_Calc( &
    & vDtbAquHkf(iSp), & !in
    & pW_ref,          & !in
    & vS_ref(iSp))
  end do
  !
  ! print *,'(A,G15.6)',"RHO_REF=    ",Rho_ref
  ! print *,'(A,G15.6)',"Alfa_ref=   ",Alfa_ref
  ! print *,'(A,G15.6)',"dAlfdT_ref= ",dAlfdT_ref
  ! pause
  !------------------------/compute properties at standard T,P cond'n --
  !
  !-------------------- compute Thermodyn. properties for all species --
  do jTP=1,NN
    !
    call DtbH2OHkf_Calc( &
    & vTPCond(jTP)%TdgC+T_CK, & !in
    & vTPCond(jTP)%Pbar, &      !in
    & pW(jTP))
    !
    Rho(jTP)=    pW(jTP)%Rho
    Alfa(jTP)=   pW(jTP)%Alfa
    dAlfdT(jTP)= pW(jTP)%dAlfdT
    Beta(jTP)=   pW(jTP)%Beta
    !
  enddo
  !
  do iSp=1,size(vDtbAquHkf)
    do jTP=1,NN
      call DtbAquHkf_Calc( &
      & vDtbAquHkf(iSp), pW(jTP), & !in
      & tS(iSp,jTP))
    enddo
  enddo
  !
  !----------------- compute logK according to Anderson Density Model --
  do iSp=1,size(vDtbAquHkf)
    do jTP=1,NN
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
      &   *(Alfa_ref *(T-T_ref) + log(Rho(jTP)/Rho_ref))
      tG0rt1(iSp,jTP)= vS_ref(iSp)%G0rt - XX /R_jk /T
      !
      !--- compute using enthalpy computed from G0 and S0
      H_ref= vH0(iSp)
      XX= H_ref *(T-T_ref) /T_ref &
      & - S_ref%Cp0 /dAlfdT_ref /T_ref &
      &   *(Alfa_ref *(T-T_ref) + log(Rho(jTP)/Rho_ref))
      tG0rt2(iSp,jTP)= vS_ref(iSp)%G0rt - XX /R_jk /T
      !
    enddo
  enddo
  !-----------------/compute logK according to Anderson Density Model --
  !
  !----------------------------------------------------------/compute --
  !
  !----------------------------------------------------------- OUTPUT --
  call GetUnit(f1)
  open(f1,file=trim(DirDtbOut)//"hkf_aqu_volumes.tab")
  call WriteTPsequence(f1,vTPCond) 
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_aqu_volumes.tab",&
  & "DTB: Hkf/ AquSpecies/ Volumes")
  !
  call GetUnit(f2)
  open(f2,file=trim(DirDtbOut)//"hkf_aqu_enthalpy.tab")
  call WriteTPsequence(f2,vTPCond)
  !-> write logK on fLogK
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_aqu_enthalpy.tab",&
  & "DTB: Hkf/ AquSpecies/ Enthalpy/ KiloJoule/Mole")
  !
  call GetUnit(f3)
  open(f3,file=trim(DirDtbOut)//"hkf_aqu_gibbs.tab")
  call WriteTPsequence(f3,vTPCond) 
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_aqu_gibbs.tab",&
  & "DTB: Hkf/ AquSpecies/ GibbsEnergy")
  !
  !--------------------------------------------- tabulate aqu'species --
  
  write(f1,'(4(A,A1))', advance="no") "AQU",T_,"NUM",T_,"NAME",T_,"Rho",T_
  do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") Rho(jTP), T_;   enddo
  write(f1,*)
  write(f1,'(4(A,A1))', advance="no") "AQU",T_,"NUM",T_,"NAME",T_,"Alfa",T_
  do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") Alfa(jTP), T_;   enddo
  write(f1,*)
  write(f1,'(4(A,A1))', advance="no") "AQU",T_,"NUM",T_,"NAME",T_,"Beta",T_
  do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") Beta(jTP), T_;   enddo
  write(f1,*)
  write(f1,'(4(A,A1))', advance="no") "AQU",T_,"NUM",T_,"NAME",T_,"dAlfdT",T_
  do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") dAlfdT(jTP), T_;   enddo
  write(f1,*)
  
  do iSp=1,size(vDtbAquHkf)
  
    M=vDtbAquHkf(iSp)
    
    !!-- write volumes
    ! write(f1,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"Vs",T_
    ! do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") tS(iSp,jTP)%Vs, T_;  enddo
    ! write(f1,*)
    write(f1,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"Vr",T_
    do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") tS(iSp,jTP)%V0, T_;   enddo
    write(f1,*)
    write(f1,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"Cp",T_
    do jTP=1,NN; write(f1,'(G15.8,A1)',advance="no") tS(iSp,jTP)%Cp0, T_;   enddo
    write(f1,*)
    
    !-- enthalpies
    write(f2,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"Hr",T_
    do jTP=1,NN; write(f2,'(G15.8,A1)',advance="no") tS(iSp,jTP)%H0/1.0D3, T_;   enddo
    write(f2,*)
    !
    !-- Gibbs/RT/ln10
    write(f3,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"logK",T_
    do jTP=1,NN
      write(f3,'(G15.8,A1)',advance="no") tS(iSp,jTP)%G0rt/log(10.), T_
    enddo
    write(f3,*)
    !
    write(f3,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"logK_1",T_
    do jTP=1,NN
      write(f3,'(G15.8,A1)',advance="no") tG0rt1(iSp,jTP)/log(10.), T_
    enddo
    write(f3,*)
    !
    write(f3,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"logK_2",T_
    do jTP=1,NN
      write(f3,'(G15.8,A1)',advance="no") tG0rt2(iSp,jTP)/log(10.), T_
    enddo
    write(f3,*)
    !
    write(f3,'(4(A,A1))', advance="no") "AQU",T_,M%Num,T_,M%Name,T_,"DELTA",T_
    do jTP=1,NN
      write(f3,'(G15.8,A1)',advance="no") (tG0rt2(iSp,jTP)-tG0rt1(iSp,jTP))/log(10.), T_
    enddo
    write(f3,*)
    !
    !!!  !write size parameter !-> not much variable ... -> not used ...
    !!!  write(f4,'(A,A1,  A15,A1,   A23,A1)',advance="no") &
    !!!  &      "AQU",T_,M%Num,T_,M%Name,T_
    !!!  do jTP=1,NN; write(f4,'(G15.8,A1)',advance="no") tS(iSp,jTP)%AquSize, T_; enddo
    !!!  write(f4,*)
    !
    !!!  !write logK's
    !!!  write(fLogK,'(A,A1,  A15,A1,   A23,A1,      A39,A1)', advance="no") &
    !!!  &        "AQU",T_,M%Num,T_,M%Name,T_,M%Formula,T_
    !!!  write(fLogK,'(F7.2,A1)',advance="no") tS(iSp,1)%AquSize, T_
    !!!  do jTP=1,NN; write(fLogK,'(G15.8,A1)',advance="no") tS(iSp,jTP)%logK, T_; enddo
    !!!  write(fLogK,*)
    !
  enddo
  !---------------------------------------------/tabulate aqu'species --
  !-----------------------------------------------------------/OUTPUT --
  
  close(f1)
  close(f2)
  close(f3)
  
  deallocate(Rho)
  deallocate(Alfa)
  deallocate(dAlfdT)
  !
  deallocate(tS)
  deallocate(vS_ref)
  deallocate(pW)
  deallocate(tG0rt1)
  deallocate(tG0rt2)
  deallocate(vH0)
  
  print '(A)',"Aqu'Species Volumes in .... "//trim(DirDtbOut)//"hkf_aqu_volumes.tab"
  print '(A)',"Aqu'Species Enthalpies in . "//trim(DirDtbOut)//"hkf_aqu_enthalpy.tab"
  print '(A)',"Aqu'Species Gibbs in ...... "//trim(DirDtbOut)//"hkf_aqu_gibbs.tab"

  return
end subroutine DtbAquHkf_Tabulate

subroutine DtbMinHkf_Tabulate(vTPCond)
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_Numeric_Const, only: Ln10
  use M_Dtb_Const,  only: T_CK
  use M_T_Tpcond,   only: T_TPCond
  use M_T_DtbMinHkf,only: T_DtbMinHkf,DtbMinHkf_Calc
  use M_T_Species,  only: T_Species
  !
  use M_Dtb_Vars,   only: vDtbMinHkf
  !
  type(T_TPCond), dimension(:),intent(in):: vTPCond
  !
  type(T_DtbMinHkf):: M
  type(T_Species)  :: S
  integer :: jTP,iSp,f2,f3,N,NN
  !
  type(T_Species),allocatable,dimension(:,:)::tS
  !
  N=  size(vDtbMinHkf)
  NN= size(vTPCond)
  allocate(tS(1:N,1:NN))
  !
  !-------------------- compute Thermodyn. properties for all species --
  do iSp=1,size(vDtbMinHkf)
    do jTP=1,NN
      !tS(iSp,jTP)=
      call DtbMinHkf_Calc( &
      & vDtbMinHkf(iSp),       &
      & vTPCond(jTP)%TdgC+T_CK,&
      & vTPCond(jTP)%Pbar,     &
      & tS(iSp,jTP))
    enddo
  enddo
  !---------------------------------------------------------/ compute --
  !
  !----------------------------------------------------------- OUTPUT --
  call GetUnit(f2)
  open(f2,file=trim(DirDtbOut)//"hkf_min_cp.tab")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_min_cp.tab",&
  & "Cp of Minerals, Hkf model")
  ! 
  call GetUnit(f3)
  open(f3,file=trim(DirDtbOut)//"hkf_min_gibbs.tab")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_min_gibbs.tab",&
  & "DTB: Gibbs energy of Minerals, Hkf model")
  !
  !
  call WriteTPsequence(f2,vTPCond)
  !
  call WriteTPsequence(f3,vTPCond)
  !
  do iSp=1,size(vDtbMinHkf)
    !
    M= vDtbMinHkf(iSp)
    S= tS(iSp,1)
    !
    write(f2,'(3(A,A1),I3,A1)',advance="no") "MIN",T_,M%Num,T_,M%Name,T_,M%Div,T_
    write(f3,'(3(A,A1))',advance="no") "MIN",T_,M%Num,T_,M%Name,T_
    do jTP=1,NN
      write(f2,'(G15.8,A1)',advance="no") tS(iSp,jTP)%Cp0, T_
      write(f3,'(G15.8,A1)',advance="no") tS(iSp,jTP)%G0rt, T_
    enddo
    write(f2,*)
    write(f3,*)
    !
  enddo
  !-----------------------------------------------------------/OUTPUT --
  !
  !close(f1)
  close(f2)
  close(f3) 
  !
  deallocate(tS)
  
  return
end subroutine DtbMinHkf_Tabulate

subroutine DtbAquThr_Tabulate(vTPCond)
  use M_Files,        only: DirDtbOut,Files_Index_Write
  use M_Numeric_Const,only: Ln10
  use M_Dtb_Const,    only: R_jk,T_CK
  use M_T_Tpcond,     only: T_TPCond
  use M_T_Species,    only: T_Species
  use M_T_DtbAquHkf
  use M_T_DtbH2OHkf,only: T_DtbH2OHkf,DtbH2OHkf_Calc
  use M_Fluid_Calc, only: Eos_H2O_Haar
  !
  use M_Dtb_Vars,   only: vDtbAquHkf
  !
  type(T_TPCond), dimension(:),intent(in):: vTPCond
  !
  integer :: jTP,iSp,f1,f2,f3,N,NN
  real(dp):: Pbar,TdgK,G_H2O,V_H2O !for H2O tabulation
  !
  real(dp):: Rho_ref,Alfa_ref,dAlfdT_ref
  real(dp):: TdgK_ref,Pbar_ref
  type(T_DtbH2OHkf):: pW_ref
  type(T_Species),allocatable:: vS_ref(:)
  !
  type(T_DtbAquHkf):: M
  type(T_Species),  allocatable:: tS(:,:)
  type(T_DtbH2OHkf),allocatable:: pW(:)
  !
  N=  size(vDtbAquHkf)
  NN= size(vTPCond)
  allocate(tS(1:N,1:NN))
  allocate(pW(1:NN))
  allocate(vS_ref(1:N))
  !
  !------------------------ compute properties at standard T,P cond'n --
  Pbar_ref= 1.0D0
  TdgK_ref= 25.0D0 +T_CK
  !
  call DtbH2OHkf_Calc( &
  & TdgK_ref,Pbar_ref,&
  & pW_ref) !out
  Rho_ref=    pW_ref%Rho
  Alfa_ref=   pW_ref%Alfa
  dAlfdT_ref= pW_ref%dAlfdT
  !
  do iSp=1,N
    call DtbAquHkf_CalcThr( &
    & vDtbAquHkf(iSp), & !in
    & pW_ref,          & !in
    & vS_ref(iSp))
  end do
  !------------------------/compute properties at standard T,P cond'n --
  !
  !-------------------- compute Thermodyn. properties for all species --
  do jTP=1,NN
    call DtbH2OHkf_Calc( &
    & vTPCond(jTP)%TdgC+T_CK,vTPCond(jTP)%Pbar,&
    & pW(jTP)) !out
  enddo
  do iSp=1,N
    do jTP=1,NN
      call DtbAquHkf_CalcThr( &
      & vDtbAquHkf(iSp), & !in
      & pW(jTP),         & !in
      & tS(iSp,jTP))
    enddo
  enddo
  !---------------------------------------------------------- compute --
  !
  !----------------------------------------------------------- OUTPUT --
  call GetUnit(f1)
  open(f1,file=trim(DirDtbOut)//"thr_aqu_volume.tab")
  call WriteTPsequence(f1,vTPCond)
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"thr_aqu_volume.tab",&
  & "DTB: Aqu.Species,HkfData,ThrConvention/ Volumes")
  !
  call GetUnit(f2)
  open(f2,file=trim(DirDtbOut)//"thr_aqu_enthalpy.tab")
  call WriteTPsequence(f2,vTPCond) 
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"thr_aqu_enthalpy.tab",&
  & "DTB: AquSpecies,HkfData,ThrConvention/  Enthalpy/ KiloJoule/Mole")
  !
  call GetUnit(f3)
  open(f3,file=trim(DirDtbOut)//"thr_aqu_gibbs.tab")
  call WriteTPsequence(f3,vTPCond)
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"thr_aqu_gibbs.tab",&
  & "DTB: logK,AquSpecies,HkfData,ThrConvention")
  !
  !------------------------------------------------ tabulate logK_H2O --
  write(fLogK,'(A3,A1,A15,A1,A23,A1,A39,A1,A7,A1)',advance="no") &
  & "AQU",T_,"THR_H2O",T_,"H2O",T_,"O(1)H(2)",T_,"   0.00",T_
  do jTP=1,NN
    Pbar=vTPCond(jTP)%Pbar
    TdgK=vTPCond(jTP)%TdgC+T_CK
    call Eos_H2O_Haar(Pbar,TdgK,G_H2O,V_H2O)
    !! G_H2O= G_H2O + 69544.987825D0
    write(fLogK,'(F12.6,A1)',advance="no") -G_H2O/R_jk/TdgK/Ln10,T_
  enddo
  write(fLogK,*)
  !-----------------------------------------------/ tabulate logK_H2O --
  !
  !write(fLogK,'(A,A1,A,A1,A7,A1)',advance="no") &
  !& "H2O",T_,"O(1)H(2)",T_,"   0.00",T_
  !do jTP=1,NN
  !  Pbar=vTPCond(jTP)%Pbar; TdgK=vTPCond(jTP)%TdgC+273.150D0
  !  call Eos_H2O_HolPow91(Pbar,TdgK,G_H2O,V); K_H2O=-G_H2O/R_jk/TdgK/log(10.0)
  !  write(fLogK,'(F12.6,A1)',advance="no") K_H2O,T_
  !enddo
  !write(fLogK,*)
  !
  !--------------------------------------------- tabulate aqu'species --
  do iSp=1,size(vDtbAquHkf)
    !
    M=vDtbAquHkf(iSp)
    !
    write(f1,'(3(A,A1))',advance="no")    "AQU",T_,M%Num,T_,M%Name,T_
    write(f2,'(3(A,A1))',advance="no")    "AQU",T_,M%Num,T_,M%Name,T_
    write(f3,'(3(A,A1))',advance="no")    "AQU",T_,M%Num,T_,M%Name,T_
    write(fLogK,'(4(A,A1))',advance="no") "AQU",T_,M%Num,T_,M%Name,T_,M%Formula,T_
    write(fLogK,'(F7.2,A1)',advance="no") tS(iSp,1)%AquSize, T_
    !
    do jTP=1,NN
      !write volumes
      write(f1,'(F12.6,A1)',advance="no") tS(iSp,jTP)%V0, T_
      !write enthalpies
      write(f2,'(G15.8,A1)',advance="no") tS(iSp,jTP)%H0/1.0D3, T_
      !write gibbs/RT
      write(f3,'(G15.8,A1)',advance="no") tS(iSp,jTP)%G0rt, T_
      !write LogK
      write(fLogK,'(G15.8,A1)',advance="no") -tS(iSp,jTP)%G0rt/Ln10, T_
    enddo
    !
    write(f1,*)
    write(f2,*)
    write(f3,*)
    write(fLogK,*)
    !
  enddo !iSp
  !--------------------------------------------/ tabulate aqu'species --
  !
  close(f1)
  close(f2)
  close(f3)
  close(fLogK)
  !-----------------------------------------------------------/OUTPUT --
  !
  deallocate(tS)
  deallocate(pW)
  
  return
end subroutine DtbAquThr_Tabulate

subroutine DtbMin_Tabulate(sCode,vTPCond) !
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_Numeric_Const, only: Ln10
  use M_Dtb_Const,  only: T_CK, Tref, R_jK
  use M_T_Tpcond,   only: T_TPCond
  use M_T_Species,  only: T_Species
  use M_T_DtbMinHkf,only: T_DtbMinHkf,DtbMinHkf_Calc
  !
  use M_T_DtbMinThr
  !
  use M_Dtb_Vars,   only: vDtbMinThr, vDtbMinHkf
  !
  character(len=*),intent(in):: sCode
  type(T_TPCond),  intent(in):: vTPCond(:)
  !
  type(T_DtbMinThr):: Mt
  type(T_DtbMinHkf):: Mh
  type(T_Species),allocatable:: tS(:,:)
  type(T_Species):: S
  integer :: iMin,jTP,nTP,nMin
  integer :: f1,f1b,f2,f3,f4,f5,f
  !! integer :: ff(1:5)
  real(dp):: TdgK, Pbar
  character(len=63):: Str
  !
  ! must have called InitTPsequence beforehand ...
  !
  call GetUnit(f1)
  open(f1,file=trim(DirDtbOut)//trim(sCode)//"_min_density.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_density.tab",&
  & "DTB: minerals,Thr,Density,kg/m3")
  call WriteTPsequence(f1,vTPCond)
  !
  call GetUnit(f1b)
  open(f1b,file=trim(DirDtbOut)//trim(sCode)//"_min_volum.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_volum.tab",&
  & "DTB: minerals,Thr,Volume,cm3")
  call WriteTPsequence(f1b,vTPCond)
  !
  call GetUnit(f2)
  open(f2,file=trim(DirDtbOut)//trim(sCode)//"_min_gibbs.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_gibbs.tab",&
  & "DTB: minerals, Thr, Gibbs, kiloJoules, Benson-Helgeson Conv.")
  call WriteTPsequence(f2,vTPCond)  
  !
  call GetUnit(f3)
  open(f3,file=trim(DirDtbOut)//trim(sCode)//"_min_cp.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_cp.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,Cp")
  call WriteTPsequence(f3,vTPCond)
  !
  call GetUnit(f4)
  open(f4,file=trim(DirDtbOut)//trim(sCode)//"_min_enthalpy.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_enthalpy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,enthalpy")
  call WriteTPsequence(f4,vTPCond)
  !
  call GetUnit(f5)
  open(f5,file=trim(DirDtbOut)//trim(sCode)//"_min_entropy.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_entropy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,entropy")
  call WriteTPsequence(f5,vTPCond)
  !
  nTP=  size(vTPCond)
  if(sCode=="THR")  nMin= size(vDtbMinThr)
  if(sCode=="HKF")  nMin= size(vDtbMinHkf)
  allocate(tS(1:nMin,1:nTP))
  !
  !-------------------- compute Thermodyn. properties for all species --
  do iMin=1,nMin!____________
    do jTP=1,nTP
      TdgK= vTPCond(jTP)%TdgC+T_CK
      Pbar= vTPCond(jTP)%Pbar
      !
      if(sCode=="THR") call DtbMinThr_Calc( &
      & vDtbMinThr(iMin), TdgK, Pbar, &
      & S)
      !
      if(sCode=="HKF") call DtbMinHkf_Calc( &
      & vDtbMinHkf(iMin), TdgK, Pbar, &
      & S)
      !
      S%G0rt= S%G0rt*TdgK*R_jK !! -Tref*S%S0Ele
      !
      tS(iMin,jTP)= S
    enddo
  enddo
  !-------------------/ compute Thermodyn. properties for all species --
  !
  do iMin=1,nMin
    if(sCode=="THR") then 
      Mt= vDtbMinThr(iMin)
      write(f1, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f1b,'(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f2, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f3, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f4, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f5, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
    end if
    if(sCode=="HKF") then
      Mh= vDtbMinHkf(iMin)
      write(f1, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f1b,'(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f2, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f3, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f4, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f5, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
    end if
    do jTP=1,nTP
      S= tS(iMin,jTP)
      write(f1, '(G15.6,A1)',advance="no") S%WeitKg / S%V0, T_ 
      write(f1b,'(G15.6,A1)',advance="no") S%V0*1.E6, T_ 
      write(f2, '(G15.6,A1)',advance="no") S%G0rt/1.E3, T_
      write(f3, '(G15.6,A1)',advance="no") S%CP0, T_ !/nOx -> "Sdatcaled" to NrOfOxygenSdat
      write(f4, '(G15.6,A1)',advance="no") S%H0/1.E3, T_
      write(f5, '(G15.6,A1)',advance="no") S%S0/1.E3, T_
    enddo
    write(f1, *)
    write(f1b,*)
    write(f2, *)
    write(f3, *)
    write(f4, *)
    write(f5, *)
    !! write(fLogK, *)
  enddo
  close(f1); close(f1b); close(f2)
  close(f3); close(f4);  close(f5)
  !
  call GetUnit(f)
  open(f,file=trim(DirDtbOut)//trim(sCode)//"_min_all.tab")
  !
  do iMin=1,nMin
    !
    if(sCode=="THR") Str= trim(vDtbMinThr(iMin)%Name)
    if(sCode=="HKF") Str= trim(vDtbMinHkf(iMin)%Name)
    !
    write(f,'(2(A,A1))',advance="no") "VOL",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f,'(G15.6,A1)',advance="no") tS(iMin,jTP)%V0*1.E6, T_ 
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "Gibbs",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%G0rt/1.E3, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "Cp",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%CP0, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "H",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%H0/1.E3, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "S",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%S0/1.E3, T_
    enddo
    write(f, *)
    !
  enddo
  !
  close(f)
  !
  deallocate(tS)
  
  return
end subroutine DtbMin_Tabulate

subroutine DtbSys_Tabulate(sCode,vTPCond) !
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_Numeric_Const, only: Ln10
  use M_Dtb_Const,  only: T_CK, Tref, R_jK
  use M_T_Tpcond,   only: T_TPCond
  use M_T_Species,  only: T_Species
  use M_T_DtbMinHkf,only: T_DtbMinHkf,DtbMinHkf_Calc
  !
  use M_T_DtbMinThr
  !
  use M_Dtb_Vars,   only: vDtbMinThr, vDtbMinHkf
  !
  character(len=*),intent(in):: sCode
  type(T_TPCond),  intent(in):: vTPCond(:)
  !
  type(T_DtbMinThr):: Mt
  type(T_DtbMinHkf):: Mh
  type(T_Species),allocatable:: tS(:,:)
  type(T_Species):: S
  integer :: iMin,jTP,nTP,nMin
  integer :: f1,f1b,f2,f3,f4,f5,f
  !! integer :: ff(1:5)
  real(dp):: TdgK, Pbar
  character(len=63):: Str
  !
  ! must have called InitTPsequence beforehand ...
  !
  call GetUnit(f1)
  open(f1,file=trim(DirDtbOut)//trim(sCode)//"_min_density.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_density.tab",&
  & "DTB: minerals,Thr,Density,kg/m3")
  call WriteTPsequence(f1,vTPCond)
  !
  call GetUnit(f1b)
  open(f1b,file=trim(DirDtbOut)//trim(sCode)//"_min_volum.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_volum.tab",&
  & "DTB: minerals,Thr,Volume,cm3")
  call WriteTPsequence(f1b,vTPCond)
  !
  call GetUnit(f2)
  open(f2,file=trim(DirDtbOut)//trim(sCode)//"_min_gibbs.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_gibbs.tab",&
  & "DTB: minerals, Thr, Gibbs, kiloJoules, Benson-Helgeson Conv.")
  call WriteTPsequence(f2,vTPCond)  
  !
  call GetUnit(f3)
  open(f3,file=trim(DirDtbOut)//trim(sCode)//"_min_cp.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_cp.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,Cp")
  call WriteTPsequence(f3,vTPCond)
  !
  call GetUnit(f4)
  open(f4,file=trim(DirDtbOut)//trim(sCode)//"_min_enthalpy.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_enthalpy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,enthalpy")
  call WriteTPsequence(f4,vTPCond)
  !
  call GetUnit(f5)
  open(f5,file=trim(DirDtbOut)//trim(sCode)//"_min_entropy.tab")
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//trim(sCode)//"_min_entropy.tab",&
  !! & "DTB: minerals,Thr,Cp,SCALED to Oxygen number")
  & "DTB: minerals,Thr,entropy")
  call WriteTPsequence(f5,vTPCond)
  !
  nTP=  size(vTPCond)
  if(sCode=="THR")  nMin= size(vDtbMinThr)
  if(sCode=="HKF")  nMin= size(vDtbMinHkf)
  allocate(tS(1:nMin,1:nTP))
  !
  !-------------------- compute Thermodyn. properties for all species --
  do iMin=1,nMin!____________
    do jTP=1,nTP
      TdgK= vTPCond(jTP)%TdgC+T_CK
      Pbar= vTPCond(jTP)%Pbar
      !
      if(sCode=="THR") call DtbMinThr_Calc( &
      & vDtbMinThr(iMin), TdgK, Pbar, &
      & S)
      !
      if(sCode=="HKF") call DtbMinHkf_Calc( &
      & vDtbMinHkf(iMin), TdgK, Pbar, &
      & S)
      !
      S%G0rt= S%G0rt*TdgK*R_jK !! -Tref*S%S0Ele
      !
      tS(iMin,jTP)= S
    enddo
  enddo
  !-------------------/ compute Thermodyn. properties for all species --
  !
  do iMin=1,nMin
    if(sCode=="THR") then 
      Mt= vDtbMinThr(iMin)
      write(f1, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f1b,'(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f2, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f3, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f4, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
      write(f5, '(3(A,A1),I3,A1)',advance="no") Mt%Typ,T_,Mt%Num,T_,Mt%name,T_,Mt%Div,T_
    end if
    if(sCode=="HKF") then
      Mh= vDtbMinHkf(iMin)
      write(f1, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f1b,'(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f2, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f3, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f4, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
      write(f5, '(3(A,A1),I3,A1)',advance="no") Mh%Typ,T_,Mh%Num,T_,Mh%name,T_,Mh%Div,T_
    end if
    do jTP=1,nTP
      S= tS(iMin,jTP)
      write(f1, '(G15.6,A1)',advance="no") S%WeitKg / S%V0, T_ 
      write(f1b,'(G15.6,A1)',advance="no") S%V0*1.E6, T_ 
      write(f2, '(G15.6,A1)',advance="no") S%G0rt/1.E3, T_
      write(f3, '(G15.6,A1)',advance="no") S%CP0, T_ !/nOx -> "Sdatcaled" to NrOfOxygenSdat
      write(f4, '(G15.6,A1)',advance="no") S%H0/1.E3, T_
      write(f5, '(G15.6,A1)',advance="no") S%S0/1.E3, T_
    enddo
    write(f1, *)
    write(f1b,*)
    write(f2, *)
    write(f3, *)
    write(f4, *)
    write(f5, *)
    !! write(fLogK, *)
  enddo
  close(f1); close(f1b); close(f2)
  close(f3); close(f4);  close(f5)
  !
  call GetUnit(f)
  open(f,file=trim(DirDtbOut)//trim(sCode)//"_min_all.tab")
  !
  do iMin=1,nMin
    !
    if(sCode=="THR") Str= trim(vDtbMinThr(iMin)%Name)
    if(sCode=="HKF") Str= trim(vDtbMinHkf(iMin)%Name)
    !
    write(f,'(2(A,A1))',advance="no") "VOL",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f,'(G15.6,A1)',advance="no") tS(iMin,jTP)%V0*1.E6, T_ 
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "Gibbs",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%G0rt/1.E3, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "Cp",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%CP0, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "H",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%H0/1.E3, T_
    enddo
    write(f, *)
    !
    write(f,'(2(A,A1))',advance="no") "S",T_,Trim(Str),T_
    do jTP=1,nTP
      write(f, '(G15.6,A1)',advance="no") tS(iMin,jTP)%S0/1.E3, T_
    enddo
    write(f, *)
    !
  enddo
  !
  close(f)
  !
  deallocate(tS)
  
  return
end subroutine DtbSys_Tabulate

subroutine DtbH2OHkf_Tabulate(vTPCond)
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_Dtb_Const,  only: T_CK
  use M_T_Tpcond,   only: T_TPCond
  use M_T_DtbH2OHkf,only: T_DtbH2OHkf,DtbH2OHkf_Calc
  use M_T_DtbAquHkf
  use M_IoTools
  !
  type(T_TPCond), dimension(:),intent(in):: vTPCond
  !
  integer::i,f1,NN
  type(T_DtbH2OHkf):: p
  type(T_DtbH2OHkf),allocatable:: vDtbH2OHkf(:)
  
  !----------------------------------------- write Solvent properties --
  !
  !------------------------------------- results of HKF-Haar routines --
  NN= size(vTPCond)
  allocate(vDtbH2OHkf(1:NN))
  
  call GetUnit(f1)
  open(f1,file=trim(DirDtbOut)//"propsh2o.restab")
  
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"propsh2o.restab",&
  & "DTB: Solvent properties, results of HKF-Haar routines")
  
  write(f1,'(3(A,A1))',advance="no")  "TdgC",T_,"Pbar", T_,"RH2O", T_
  write(f1,'(4(A,A1))',advance="no")  "G",   T_,"dG_dP",T_,"dG_dT",T_,"d2G_dT2",T_
  !write(f1,'(4(A,A1))',advance="no") "Z",   T_,"Q",    T_,"Y",    T_,"X",T_
  write(f1,'(2(A,A1))')               "dhA", T_,"dhB"
  
  do I=1,NN
    call DtbH2OHkf_Calc( &
    & vTPCond(i)%TdgC+T_CK,vTPCond(i)%Pbar, &
    & vDtbH2OHkf(i))
  enddo
  
  do I=1,NN
    p= vDtbH2OHkf(i)
    write(f1,'(F7.3,A1,F7.3,A1,G15.8,A1)',advance="no") &
    & p%TdgK-T_CK,T_, &
    & p%Pbar,     T_, &
    & p%Rho,      T_
    write(f1,'(4(G15.8,A1))',advance="no") &
    & p%Gshok,T_,p%dGShdP,T_,p%dGShdT,T_,p%d2GShdT2,T_
    write(f1,*)
  enddo !I
  
  close(f1)
  !----------------------------------------/ write Solvent properties --
  !
  write(fLogK,'(/,A)') "TP.TABLE"
  !
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%TdgK-T_CK,Opt_S="TdgC")
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Pbar,     Opt_S="Pbar")
  !
  write(fLogK,'(A,/)') "END TP.TABLE"
  !
  write(fLogK,'(/,A)') "SOLVENT"
  !
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Rho, Opt_S="Rho")
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%Eps, Opt_S="Eps")
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%dhA, Opt_S="DHA")
  call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%dhB, Opt_S="DHB")
  !! call OutStrVec(fLogK,vDtbH2OHkf(1:NN)%bdot,Opt_S="BdoT")
  !
  write(fLogK,'(A,/)') "END"
  !
  !call OutStrVec(f1,S="Gshok__",vDtbH2OHkf(1:NN)%Gshok)
  !call OutStrVec(f1,S="dGShdP_",vDtbH2OHkf(1:NN)%dGShdP)
  !call OutStrVec(f1,S="dGShdT_",vDtbH2OHkf(1:NN)%dGShdT)
  !call OutStrVec(f1,S="d2GSdT2",vDtbH2OHkf(1:NN)%d2GShdT2)
  !call OutStrVec(f1,S="dEpsdP_",vDtbH2OHkf(1:NN)%dEpsdP)
  !call OutStrVec(f1,S="dEpsdT_",vDtbH2OHkf(1:NN)%dEpsdT)
  !call OutStrVec(f1,S="d2EpsdT",vDtbH2OHkf(1:NN)%d2EpsdT2)
  !call OutStrVec(f1,S="Qeps___",vDtbH2OHkf(1:NN)%Qeps)
  !call OutStrVec(f1,S="Xeps___",vDtbH2OHkf(1:NN)%Xeps)
  !call OutStrVec(f1,S="Yeps___",vDtbH2OHkf(1:NN)%Yeps)
  !call OutStrVec(f1,S="Zeps___",vDtbH2OHkf(1:NN)%Zeps)
  !
end subroutine DtbH2OHkf_Tabulate

subroutine WriteTPsequence(f,vTPCond)
  
  use M_T_Tpcond,only: T_TPCond
  !
  type(T_TPCond), dimension(:),intent(in):: vTPCond
  !
  integer:: f
  integer:: i,NN
  
  NN= size(vTPCond)
  
  write(f,'(4(A6,A1))',advance="no") &
  & ".TYP",T_,".SOURCE",T_,".SPECIES",T_,".TDGC",T_
  do I=1,NN; write(f,'(G15.4, A1)',advance="no") vTPCond(I)%TdgC,T_; enddo; write(f,*)
  
  write(f,'(4(A6,A1))',advance="no") &
  & ".TYP",T_,".SOURCE",T_,".SPECIES",T_,".PBAR",T_
  do I=1,NN; write(f,'(G15.4, A1)',advance="no") vTPCond(I)%Pbar,T_; enddo; write(f,*)
  
  return
end subroutine WriteTPsequence

subroutine DtbMinHkf_Write
!--
!-- check result DtbMinHKF_Read
!--
  use M_Files,    only: DirDtbOut,Files_Index_Write
  use M_Dtb_Const,only: Tref
  use M_T_DtbMinHkf
  use M_Dtb_Vars, only: vDtbMinHkf
  !
  type(T_DtbMinHkf):: M
  integer       :: I,f
  real(dp)      :: X
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbMinHkf_Write"
  !
  call GetUnit(f)
  open(f,file=trim(DirDtbOut)//"hkf_min.restab")
  !do I=1,nEleDtb; write(f,'(G15.8,A1)',advance="no") vS0Ele(I),T_; enddo; write(f,*)
  
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"hkf_min.restab",&
  & "DTB: table of data for Min.species, HKF")
  
  write(f,'(3(A,A1))',advance="no") &
  & "Num",     T_,&
  & "Name",    T_,&
  & "Formula", T_
  write(f,'(6(A,A1))',advance="no") &
  & "(G0R-H0R)/Tr+S0",T_, &
  & "S0Ele",          T_, &
  & "G0R",            T_, &
  & "H0R",            T_, &
  & "S0",             T_, &
  & "V0R",            T_
  write(f,'(4(A,A1))') &
  & "MK1_1__",T_,&
  & "MK1_2__",T_,&
  & "MK1_3__",T_,&
  & "MK1_4__",T_
  
  do I=1,size(vDtbMinHkf)
  
    M= vDtbMinHkf(I)
    X= (M%G0R -M%H0R)/Tref +M%S0_ 
    
    write(f,'(3(A,A1))',advance="no") &
    & M%Num,    T_, &
    & M%Name,   T_, &
    & M%Formula,T_
    write(f,'(6(G15.8,A1))',advance="no") &
    & X,       T_,  &
    & M%S0Ele, T_,  &
    & M%G0R,   T_,  &
    & M%H0R,   T_,  &
    & M%S0_,   T_,  &
    & M%V0R,   T_
    write(f,'(3(G15.8,A1))') &
    & M%MK1(1),T_, &
    & M%MK1(2),T_, &
    & M%MK1(3),T_
  
  enddo
  
  close(f)
  
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbMinHkf_Write"
end subroutine DtbMinHkf_Write

subroutine DtbAquHkf_Write
!--
!-- check result DtbAquHKF_Read
!-- make a new version of Hkf base, one species per line
!--

  use M_Dtb_Const,only: Tref
  use M_Files,    only: DirDtbOut,Files_Index_Write
  use M_T_DtbAquHkf
  use M_Dtb_Vars, only: vDtbAquHkf
  !
  type(T_DtbAquHkf):: M
  integer:: I,f
  
  
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbAquHkf_Write"
  !
  call GetUnit(f)
  open(f,file=trim(DirDtbOut)//"_aqu_hkf.restab")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"aqu_hkf.restab",&
  & "DTB: table of data for Aqu.species, HKF")
  !
  write(f,'(14(A,A1))') &
  & "Num",T_,"Name",T_,"Formula",T_, &
  !!& "(G0R-H0R)/Tr+S0",T_, "S0Ele",T_,&
  & "GfAq",T_, "HfAq",T_, "SPrTrAq",T_, &
  & "A1",T_,"A2",T_,"A3",T_,     "A4",T_, &
  & "C1",T_,"C2",T_,"Omega__",T_,"Chg",T_
  !
  do I=1,size(vDtbAquHkf)
    M=vDtbAquHkf(I)
    write(f,'(3(A,A1))',advance="no") &
    & M%Num,T_,&
    & M%Name,T_,&
    & M%Formula,T_
    !!Y=dot_product(M%Stoik(1:nEl),vEle(1:nEl)%S0)
    !!write(f,'(5(G15.8,A1))',advance="no") X,   T_,Y,   T_,M%G0R, T_,M%H0R,T_,M%S0_,T_
    write(f,'(3(G15.8,A1))',advance="no") M%G0R, T_,M%H0R,T_,M%S0_,T_
    write(f,'(5(G15.8,A1))',advance="no") M%A1,T_,M%A2,T_,M%A3,  T_,M%A4, T_
    write(f,'(3(G15.8,A1),I3,A1)')        M%C1,T_,M%C2,T_,M%wref,T_,M%Chg,T_
  enddo
  !
  close(f)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbAquHkf_Write"
  
  return
end subroutine DtbAquHkf_Write

subroutine DtbMinThr_Write
!--
!--check result DtbMinHKF_Read
!--
  
  use M_Files,      only: DirDtbOut,Files_Index_Write
  use M_Dtb_Vars,   only: vDtbMinThr
  use M_T_DtbMinThr,only: T_DtbMinThr
  
  type(T_DtbMinThr):: M
  integer          :: I,f
  
  !!=parameters=========================================================
  !  integer :: nLanda !number of landau transitions
  !  integer :: codVol !code for volume function
  !  real(dp)::& !for landau transitions
  !  & TQ1B(1:3),TRE(1:3), ASPK(1:3),BSPK(1:3),DHTR(1:3),&
  !  & TEQ(1:3), DVTR(1:3),DVDT(1:3),DVDP(1:3)
  !  real(dp):: G0R,H0R,S0_,V0R
  !  real(dp):: S0Ele,WeitKg
  !  real(dp):: AA0,AAT,BB0,BBT
  !  !          !for VDW/RDK/... cubic EoS of Gases (VanDerWaals,RedlichKwong,...)
  !  real(dp):: Tcrit,Pcrit,Acentric !for SRK/PRS/... cubic EoS of gases
  !  real(dp):: K1,K2,K3,K4,K5,K6,K7,K8,K9  !coeff's for Cp formulas
  !  real(dp):: D1,D2,D3,D4,D5,D6,D7,D8,D9  !disorder /Berman (HP:D1,D2,D3)
  !  real(dp):: TD0,TDMax,VAdj              !disorder /Berman 
  !  real(dp):: VTA,VTB,VPA,VPB             !volume function
  !  !__________________________a la Berman: use VTA,VTB,VPA,VPB
  !  !__________________________a la HP:     use VTA,VTB,VPA,TKRI,SMA
  !  !VL0,VLA,VLN,VL2,VAA,VAB,VB, & !older formats ??
  !  real(dp):: TKri,SMA
  !  logical::Dis !,Min,VdW,RdK !,TL1
  !!--------------------------------------------------------------------
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DtbMinThr_Write"
  !
  call GetUnit(f)
  open(f,file=trim(DirDtbOut)//"_min_thr.restab")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirDtbOut)//"_min_thr.restab",&
  & "DTB: data for Min.species, Theriak, in spreadsheet form")
  !
  write(f,'(2(A,A1))',advance="no") &
  & "Name                ",T_,& !len=20
  & "Formula                                ",T_
  write(f,'(4(A,A1))',advance="no") &
  & "G0R(J)         ",T_,&
  & "H0R(J)         ",T_,&
  & "S0(J/K)        ",T_,&
  & "V0R(J/bar)     ",T_
  write(f,'(9(A,A1))',advance="no") &
  & "K1             ",T_,&
  & "K4             ",T_,&
  & "K3             ",T_,&
  & "K8             ",T_,&
  & "K6             ",T_,&
  & "K2             ",T_,&
  & "K5             ",T_,&
  & "K7             ",T_,&
  & "K9             ",T_
  write(f,'(6(A,A1))',advance="no") &
  & "VTA*1E5         ",T_,&
  & "VTB*1E5         ",T_,&
  & "VPA*1E5         ",T_,&
  & "VPB*1E8         ",T_,&
  & "Tkri           ",T_,&
  & "SMA            ",T_
  write(f,'(9(A,A1))',advance="no") &
  & "D1             ",T_,&
  & "D4             ",T_,&
  & "D3             ",T_,&
  & "D8             ",T_,&
  & "D6             ",T_,&
  & "D2             ",T_,&
  & "D5             ",T_,&
  & "D7             ",T_,&
  & "D9             ",T_
  write(f,*)
  !
  do I=1,size(vDtbMinThr)
    M=vDtbMinThr(I)
    !
    write(f,'(2(A,A1))',advance="no") M%Name,T_,M%Formula,T_
    !
    write(f,'(4(G15.8,A1))',advance="no")  &
    & M%G0R,T_,M%H0R,T_,M%S0_,T_,M%V0R,T_
    !
    write(f,'(9(G15.8,A1))',advance="no")  &
    & M%K1, T_,M%K4, T_,M%K3, T_,M%K8, T_, &
    & M%K6, T_,M%K2, T_,M%K5, T_,M%K7, T_,M%K9, T_
    !
    write(f,'(6(G15.8,A1))',advance="no")  &
    & M%VTA*1.0D5,T_,M%VTB*1.0D5,T_, &
    & M%VPA*1.0D5,T_,M%VPB*1.0D8,T_, &
    & M%Tkri,T_,     M%SMA,T_
    !
    write(f,'(9(G15.8,A1))',advance="no")  &
    & M%D1,T_,M%D2,T_,M%D3,T_,M%D4,T_, &
    & M%D5,T_,M%D6,T_,M%D7,T_,M%D8,T_,M%D9,T_ 
    !
    write(f,*)
  enddo
  !
  close(f)
  !
  call Warning_("Data written to "//trim(DirDtbOut)//"_min_thr.restab")
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DtbMinThr_Write"
end subroutine DtbMinThr_Write

end module M_Dtb_Write
