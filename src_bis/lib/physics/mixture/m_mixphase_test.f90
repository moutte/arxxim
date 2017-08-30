module M_MixPhase_Test
!--
!--- tools for testing mixture phases
!--- (mixture- phase of variable composition, following a symmetric model)
!--
  use M_Kinds
  use M_Trace
  !
  implicit none
  !
  private
  !
  public:: MixPhase_Minimize
  public:: MixPhase_Test
  public:: MixPhase_Show
  
  real(dp),dimension(:,:),allocatable:: tGMix,tGtot,tGMarg
  real(dp),dimension(:,:),allocatable:: tActiv1,tActiv2
  real(dp),dimension(:,:),allocatable:: tGamma1,tGamma2
  real(dp),dimension(:,:),allocatable:: tAtom
  !
  integer,parameter:: nPt=99

contains

subroutine MixPhase_Minimize
  use M_Dtb_Const,  only: T_CK
  use M_T_MixModel, only: T_MixModel
  use M_IoTools,    only: GetUnit
  !
  use M_Global_Vars,only: vSpc,vMixModel
  !
  use M_Optimsolver_Theriak
  use M_MixModel_Optim
  !
  type(T_MixModel):: MM
  real(dp),allocatable:: vX(:)
  real(dp),allocatable:: vMu(:) !
  real(dp):: G
  integer :: N
  integer :: i,j,k
  real(dp):: TdgK,Pbar
  !
  real(dp):: TolX,DeltaInit
  integer :: its,nCallG
  integer :: f,ff
  logical :: Converge
  character(len=80):: sFMT
  !
  real(dp),parameter:: &
  & Tmin= 600.0D0,   &
  & Tmax= 1500.0D0,  &
  & Tstp= 50.0D0
  !
  if(size(vMixModel)<1) then
    call Warning_("NO MIXING MODEL TO BE TESTED")
    return
  end if
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< MixPhase_Minimize"
  !
  if(iDebug>2) then
    call GetUnit(ff)
    open(ff,file= 'optimsolver_theriak.log')
  end if
  !
  call GetUnit(f)
  open(f,file= 'out_mixminim.log')
  !
  TolX= 1.0D-3
  DeltaInit= 0.05D0
  !
  MM= vMixModel(1)
  N= MM%NPole
  !
  write(sFMT,'(a,i3,a,i3,a)') &
  & '(3(G15.6,1X),', N,'(G15.6,1X),', N,'(G15.6,1X))'
  !write(sFMT,'(a,i3,a)') '(2G15.6,',N,'(G15.6,1X))'
  !print *,sFMT  ;  pause
  !
  allocate(vX(1:N))
  allocate(vMu(1:N))
  !
  TdgK= Tmin +T_CK
  Pbar= 1.0D3
  !
  allocate(Mixmodel_Optim_vMu0rt(N))
  allocate(Mixmodel_Optim_vLPole(N))
  
  DoTP: do
  
    call MixModel_Optim_SetParams(TdgK,Pbar,MM)
    do i=1,N
      Mixmodel_Optim_vMu0rt(i)= vSpc(MM%vIPole(i))%G0rt
      Mixmodel_Optim_vLPole(i)= .true.  !! MM%vHasPole(i)   !!
    end do
    
    do i=1,N
      !
      vX(i)= One - 1.0D-3
      do j=1,N
        if(j/=i) vX(j)= 1.0D-3/real(N-1)
      enddo
      !
      call Optimsolver_Theriak( & !
      !& MixModel_Optim_GetGibbs,      & !
      !& MixModel_Optim_GetPotentials, & !
      & MixModel_Optim_GetMu, & !
      & ff,            & ! IN
      & DeltaInit,     & ! IN
      & TolX,          & ! IN
      & vX,            & ! INOUT
      & vMu,           & ! OUT
      & G,             & ! OUT
      & Converge,      & ! OUT
      & its,nCallG)      ! OUT
      !
      write(f,sFMT) TdgK-T_CK,Pbar,G,(vX(k),k=1,N),(vMu(k),k=1,N)
      !write(f,sFMT) TdgK-T_CK,Pbar,(vX(k),k=1,N)
      !
    enddo
    
    TdgK= TdgK +Tstp
    if(TdgK -T_CK > Tmax) exit DoTP
    
    !~ pause
  
  enddo DoTP
  !
  deallocate(vX)
  deallocate(vMu)
  deallocate(Mixmodel_Optim_vMu0rt)
  deallocate(Mixmodel_Optim_vLPole)
  !
  if(iDebug>2) close(ff)
  close(f)
  !
  write(6,'(A)') 'resultS IN out_mixminim.log'
  call Pause_
  !
  if(iDebug>0) write(fTrc,'(/,A)') "</ MixPhase_Minimize"
  !
end subroutine MixPhase_Minimize

subroutine MixPhase_Show
!--
!-- check the mixture model for the phase vMixFas(1) of the current system
!--
  use M_IoTools,only: GetUnit,OutStrVec
  use M_Files,  only: DirOut
  !
  use M_T_MixModel,only: T_MixModel,T_Margul,Mix_Site
  use M_T_MixPhase !, only: T_MixPhase, MixPhase_NormCompo
  !
  !--- global variables
  use M_Global_Vars,only: vMixFas,vMixModel
  !---/
  !
  type(T_MixModel):: S
  type(T_MixPhase):: P
  type(T_Margul)  :: M
  integer:: I,J,K
  integer:: f1
  integer:: iP,iEl
  !
  if(size(vMixFas)<1) then
    call Warning_("NO MIXTURE PHASE TO BE TESTED")
    return
  end if
  !
  P= vMixFas(1)
  S= vMixModel(P%iModel)
  !
  call GetUnit(f1)
  open(f1,file=trim(DirOut)//"_check.res")
  !
  write(f1,'(2A,/)') "Model= ",trim(S%Name)
  !
  if(S%Model==Mix_Site) then
    !
    write(f1,'(A)') "==tPoleAtom=="
    !
    write(f1,'(A12,1X)',advance="NO") "_"
    do iEl=1,S%NAtom
      write(f1,'(1X,A6,1X)',advance="NO") S%vNamAtom(iEl)
    enddo
    write(f1,*)
    !
    do iP=1,S%NPole
      write(f1,'(A12,1X)',advance="NO") S%vNamPole(iP)
      do iEl=1,S%NAtom
        write(f1,'(I7,1X)',advance="NO") S%tPoleAtom(iP,iEl)
      enddo
      write(f1,*)
    enddo
    !
    write(f1,'(A)')  "==vIAtomSite=="
    write(f1,'(A12,1X)',advance="NO") "_"
    do iEl=1,S%NAtom
      write(f1,'(I7,1X)',advance="NO") S%vIAtomSite(iEl)
    enddo
    write(f1,*)
    !
    write(f1,'(A)')  "==vAtomMulti=="
    write(f1,'(A12,1X)',advance="NO") "_"
    do iEl=1,S%NAtom
      write(f1,'(I7,1X)',advance="NO") S%vAtomMulti(iEl)
    enddo
    write(f1,*)
    !
    write(f1,'(A)') "==tPoleAtom (fractions)=="
    do iP=1,S%NPole
      write(f1,'(A12,1X)',advance="NO") S%vNamPole(iP)
      do iEl=1,S%NAtom
        write(f1,'(F7.2,1X)',advance="NO") &
        & S%tPoleAtom(iP,iEl) /real(S%vAtomMulti(iEl))
      enddo
      write(f1,*)
    enddo
    !
    write(f1,'(/,A)') "==Normalization=="
    do iP=1,S%NPole
      write(f1,'(A12,1X,F7.2)') S%vNamPole(iP), S%vPoleCoeff(iP)
    enddo
    write(f1,*)
    !
  end if
  !
  if(S%NMarg>0) then
    write(f1,'(A)') "==Margules=="
    do I=1,S%NMarg
      M= S%vMarg(I)
      J= S%vIMargSite(I)
      write(f1,'(A,I1,2A)',advance="NO") "Site= ",J,"= ",S%vNamSite(J)
      write(f1,'(A)',advance="NO") ",  degrees= "
      do K=1,S%NAtom
        write(f1,'(1X,I2)',advance="NO") M%vDegree(K)
      enddo
      write(f1,*)
    enddo
  end if
  !
end subroutine MixPhase_Show

subroutine MixPhase_Test(Cod)
!--
!-- test the mixture model for the phase vMixFas(1) of the current system
!-- using data from theriak database for end members
!-- !!! must have initialized vDtbMinThr before calling MixPhase_Test_1 !!!
!--
  !
  use M_Dtb_Const,  only: R_jk,T_CK,Tref,Pref
  use M_Dtb_Calc,   only: SpeciesMin_TP_Update_fromDtb
  use M_T_Species,  only: Species_Index
  use M_TPcond_Read,only: TPpath_Read
  use M_T_MixModel, only: T_MixModel,Mix_Site
  use M_T_MixModel, only: MixModel_Margul_Wg_Calc,MixModel_XPoleToXSite
  use M_T_MixPhase, only: T_MixPhase, MixPhase_NormCompo
  use M_T_MixPhase, only: MixPhase_CalcActivs,MixPhase_CalcMixing
  !
  !--- global variables --
  use M_Global_Vars,only: vEle,vSpcDtb,vSpc,vMixFas,vMixModel
  use M_Path_Vars,  only: vTPpath
  !---/
  !
  character(len=*):: Cod
  !
  type(T_MixModel):: S
  type(T_MixPhase):: P
  real(dp):: TdgK,Pbar
  real(dp):: GTotRT,GMecaRT,GIdMix,GMix,GXsMix
  real(dp):: GMix_FromAct,GXsMix_FromAct
  integer :: iP,jTP,iX,N
  integer :: I,I1,I2
  !
  real(dp),allocatable:: vXAtom(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< MixPhase_Test "//trim(Cod)
  !
  TdgK= Tref ! default
  Pbar= Pref ! values
  !
  call TPpath_Read(TdgK,Pbar)
  N= size(vTPpath)
  !
  allocate(tGMix(N,nPt))
  allocate(tGTot(N,nPt))
  allocate(tGMarg(N,nPt))
  !
  allocate(tActiv1(N,nPt))
  allocate(tActiv2(N,nPt))
  !  
  allocate(tGamma1(N,nPt))
  allocate(tGamma2(N,nPt))
  !
  if(size(vMixFas)<1) then
    call Warning_("No Mixture phase in current run")
    return
  end if
  !
  P= vMixFas(1)
  !
  P%vLPole(:)= .false.
  P%vXPole(:)= Zero
  !
  S= vMixModel(P%iModel)
  !
  if(S%NPole==1) return
  !
  !--------------------- update the vector S%vIPole according to vSpc --
  do iP=1,S%NPole
    I= Species_Index(S%vNamPole(iP),vSpc)
    if(I<1) then
      call Stop_("Species "//trim(S%vNamPole(iP))//" NOT in this database !!!")
    else
      S%vIPole(iP)= I !J
    end if
  enddo
  vMixModel(P%iModel)= S
  !
  if(S%Model==Mix_Site) allocate(tAtom(nPt,S%NAtom))
  if(S%Model==Mix_Site) allocate(vXAtom(S%NAtom))
  !
  I1=1  ; I2=3    !-> test for feldspar_ss, albite(=1)-anorthite(=3) join
  I1=2  ; I2=1    !-> test for binary solution
  P%vLPole(I1)=.true.
  P%vLPole(I2)=.true. !; S%vLPole(3)=.true.
  !
  LoopTP: do jTP=1,N
    !
    Pbar= vTPpath(jTP)%Pbar
    TdgK= vTPpath(jTP)%TdgC +T_CK
    !
    if(iDebug>0) write(*,'(I3,1X,G12.3,1X,G12.3)') jTP,Pbar,TdgK
    !
    !----------------- update values of GR for the end-members at T,P --
    do I=1,S%NPole
      call SpeciesMin_TP_Update_fromDtb(TdgK,Pbar,vSpcDtb,vSpc(S%vIPole(I)))
    enddo
    !-----------------/
    !
    !----------------- update Margules parameters for the current T,P --
    if(S%NMarg>0) call MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
    !-----------------/
    !
    !-------------------- test calculation of Gibbs energy / Activity --
    !S%vXPole(3)=0.3D0
    !! P%vXPole(1)= One !-S%vXPole(3)
    !! P%vXPole(2)= 0.0D0
    LoopX: do iX=1,nPt !increment 0.01 is for nPt=99 !!! 
      !
      P%vXPole(I1)= One - iX*0.01D0 !P%vXPole(1)-0.01D0
      P%vXPole(I2)=       iX*0.01D0 !P%vXPole(2)+0.01D0
      !
      call MixPhase_NormCompo(P)
      !
      if(S%Model==Mix_Site) then
        call MixModel_XPoleToXSite( &
        & S, &
        & P%vXPole(1:S%NPole), & 
        & vXAtom(1:S%NAtom))
        if(jTP==1) tAtom(iX,:)= vXAtom(:)
        !do I=1,S%NAtom; write(fTrc,'(F12.3,A1)',advance="no") P%vXAtom(I),T_; enddo
        !write(fTrc,*)
      end if
      !
      !
      call MixPhase_CalcActivs(TdgK,Pbar,S,P)
      !
      !--- save in tables for printing
      tActiv1(jTP,iX)= exp(P%vLnAct(I1))
      tActiv2(jTP,iX)= exp(P%vLnAct(I2))
      tGamma1(jTP,iX)= exp(P%vLGam(I1))
      tGamma2(jTP,iX)= exp(P%vLGam(I2))
      !---/
      !
      call MixPhase_CalcMixing( &
      & TdgK,Pbar,S,P, &    !IN
      & GMix,GIdMix,GXsMix) !OUT
      ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
      !
      GMecaRT= Zero
      do I=1,S%NPole
        if(P%vLPole(I)) &
        & GMecaRT= GMecaRT &
        &        + P%vXPole(I) *vSpc(S%vIPole(I))%G0rt
      enddo
      GTotRT= GMecaRT + GMix/TdgK/R_jk
      !
      !--- save in tables for printing
      tGMix(jTP,iX)=  GMix   /TdgK/R_jk
      tGMarg(jTP,iX)= GXsMix /TdgK/R_jk
      tGTot(jTP,iX)=  GTotRT
      !---/
      !
      GMix_FromAct=   P%vXPole(I1)*P%vLnAct(I1) + P%vXPole(I2)*P%vLnAct(I2)
      GXsMix_FromAct= P%vXPole(I1)*P%vLGam(I1)  + P%vXPole(I2)*P%vLGam(I2)
      !print '(A,3G15.6)',"CHECK XsMix=", GXsMix_FromAct, GXsMix/TdgK/R_jk
      !print '(A,3G15.6)',"CHECK Mix=  ", GMix_FromAct,   GMix/TdgK/R_jk
      !
    enddo LoopX
      !pause
    !------------------ / test calculation of Gibbs energy / Activity --
    !
  enddo LoopTP
  !
  call WriteResults(S)
  !
  deallocate(tGMix,tGtot,tGMarg,tActiv1,tActiv2,tGamma1,tGamma2)
  !
  if(S%Model==Mix_Site) deallocate(vXAtom,tAtom)
  !
  deallocate(vTPpath)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ MixPhase_Test_"//trim(Cod)
  !
  return
end subroutine MixPhase_Test

subroutine WriteResults(S)
  use M_IoTools,   only: GetUnit,OutStrVec
  use M_Numeric_Const,only: Ln10
  use M_Files,     only: DirOut
  use M_Path_Vars, only: vTPpath
  use M_T_MixModel,only: T_MixModel,Mix_Site
  !
  type(T_MixModel),intent(in):: S
  !
  integer:: f1,f1x,f21,f22,f31,f32,fAtom
  integer:: N,I,J,iX
  character(len=30):: sFormat
  !
  ! line= TP-path
  ! column= composition (from 1% to 99% end-member 2)
  !
  ! file out_mix_gmix.restab  = G_mixing and G_mixing_xs (on different lines)
  ! file out_mix_activ.restab = activities (e-m 1, then e-m 2, on same mine)
  !
  call GetUnit(f1);  open(f1, file=trim(DirOut)//"_mix_gmix.restab")
  call GetUnit(f1x); open(f1x,file=trim(DirOut)//"_mix_gxs.restab")
  call GetUnit(f21); open(f21,file=trim(DirOut)//"_mix_activ1.restab")
  call GetUnit(f22); open(f22,file=trim(DirOut)//"_mix_activ2.restab")
  call GetUnit(f31); open(f31,file=trim(DirOut)//"_mix_gamma1.restab")
  call GetUnit(f32); open(f32,file=trim(DirOut)//"_mix_gamma2.restab")
  !
  N= size(vTPpath)
  !
  if(S%Model==Mix_Site) then
    call GetUnit(fAtom)
    open(fAtom,file=trim(DirOut)//"_mix_atom.restab")
    do I=1,S%NAtom; write(fAtom,'(A,A1)',advance="no") S%vNamAtom(I),T_; enddo
    write(fAtom,*)
  end if
  !
  write(sFormat,'(A,I3,A)') '(2(A,A1),',N,'(A15,A1))'
  !
  !------------------------------------------------ header Gibbs file --
  write(f1,sFormat) ".WHAT",T_,"index",T_, (vTPpath(J)%Name,T_,J=1,N)
  !---------------------------------------------- header GibbsXs file --
  write(f1x,sFormat) ".WHAT",T_,"index",T_,(vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Activ file --
  write(f21,sFormat) ".WHAT",T_,"iX",T_,   (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Activ file --
  write(f22,sFormat) ".WHAT",T_,"iX",T_,   (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Gamma file --
  write(f31,sFormat) ".Typ",T_,"iX",T_,    (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Gamma file --
  write(f32,sFormat) ".Typ",T_,"iX",T_,    (vTPpath(J)%Name,T_,J=1,N)
  !
  !---------------------------------------------------- write resultS --
  write(sFormat,'(A,I3,A)') '(A,A1,I3,A1,',N,'(G15.8,A1))'
  !print *,sFormat  ;  pause
  !
  LoopX_: do iX=1,nPt
    ! 
    !------------------------------- total Gibbs energy of mixing /RT --
    write(f1,sFormat)  "GMix/T",T_,  iX,T_, (tGMix(J,iX),T_,J=1,N)
    !------------------------------ excess Gibbs energy of mixing /RT --
    write(f1x,sFormat) "GXsMix/T",T_,iX,T_, (tGMarg(J,iX),T_,J=1,N)
    !
    !------------------------------------------------------- activity --
    ! results for component 1
    write(f21,sFormat) "ACTIV1",T_,  iX,T_, (tActiv1(J,iX),T_,J=1,N)
    ! results for component 2
    write(f22,sFormat) "ACTIV2",T_,  iX,T_, (tActiv2(J,iX),T_,J=1,N)
    !
    !---------------------------------------------------------- gamma --
    ! results for component 1
    write(f31,sFormat) "GAMMA1",T_,  iX,T_, (tGamma1(J,iX),T_,J=1,N)
    ! results for component 2
    write(f32,sFormat) "GAMMA2",T_,  iX,T_, (tGamma2(J,iX),T_,J=1,N)
    !
    if(S%Model==Mix_Site) then
      do I=1,S%NAtom; write(fAtom,'(F12.3,A1)',advance="no") tAtom(iX,I),T_; enddo
      write(fAtom,*)
    end if
  enddo LoopX_
  !
  close(f1)   ;  close(f1x)
  close(f21)  ;  close(f22)
  close(f31)  ;  close(f32)
  if(S%Model==Mix_Site) close(fAtom)
  !-------------------------------------------------- / write resultS --
  !
  if(iDebug>0) write(fTrc,'(A)') "results in files "//trim(DirOut)//"_mix_*.restab"
  if(iDebug>0) print      '(A)', "results in files "//trim(DirOut)//"_mix_*.restab"
  if(iDebug>2) call Pause_
  !
end subroutine WriteResults

end module M_MixPhase_Test
