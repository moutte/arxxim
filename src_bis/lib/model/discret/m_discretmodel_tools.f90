module M_DiscretModel_Tools
!--
!-- tools for "discretization" of a mixture phase
!-- to an array of pure phases
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_
  implicit none
  !
  private
  !
  public:: DiscretParam_Init
  public:: DiscretSpecies_TP_Update
  public:: DiscretSpecies_Stoikio_Calc
  !
contains

subroutine DiscretModel_Init( &
& vEle,         & !IN
& vSpc,         & !IN
& MixModel,     & !IN
& DiscretModel, & !IN
& vParam,       & !OUT
& vSpcOut)        !OUT properties of species resulting from discretization
!--
!--!!! called by DiscretParam_Init !!!
!--
!-- discretize one mixture phase, following the mixing model MixModel,
!-- into nDiscret pure species (for a binary solution);
!-- -> for each species generated,
!--    build the stoichiometric formula, e.g. CA(10)MG(01)FE(09)SI(20)O(60)/(10)
!--    build the species name from the names of the end members, e.g. DIO01HED09
!--
  use M_T_Element, only: T_Element
  use M_T_Species, only: T_Species, Species_Zero
  use M_T_MixPhase,only: T_MixPhase, MixPhase_NormCompo
  use M_T_MixModel,only: T_MixModel, MixModel_Index
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  !
  type(T_Element),    intent(in) :: vEle(:)
  type(T_Species),    intent(in) :: vSpc(:)
  type(T_MixModel),   intent(in) :: MixModel
  type(T_DiscretModel),intent(in):: DiscretModel
  !
  type(T_DiscretParam),intent(out):: vParam(:)
  type(T_Species),     intent(out):: vSpcOut(:)
  !
  type(T_MixPhase):: P
  type(T_Species) :: S
  ! character(len=7):: cName
  ! integer      :: nDiscret
  integer :: nEl,N,I,J,K,i0,Dims
  ! integer :: Z,nDiv1,nDiv2 !,nDiv3
  real(dp):: ntot
  logical :: Ternary
  character(len=23):: DiscretName
  character(len=71):: Formul
  integer,allocatable:: vStoik1(:),vStoik2(:),vStoik3(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_Init"
  !
  Ternary= (MixModel%NPole==3)
  !! MixModel= vMixModel(DiscretModel%iMix)
  Dims= DiscretModel%Dim1 +1
  !
  P%Name=      trim(DiscretModel%Name)
  P%iModel=    DiscretModel%iMix !!MixModel_Index(vMixModel,MixModel%Name)
  P%vLPole=    .false.
  P%vXPole=    Zero
  !
  ! compute stoikio vectors of end members
  if(iDebug>0) write(fTrc,'(A)') trim(vSpc(MixModel%vIPole(1))%Formula)
  if(iDebug>0) write(fTrc,'(A)') trim(vSpc(MixModel%vIPole(2))%Formula)
  !
  if(MixModel%NPole==1) return
  !
  nEl= size(vEle)
  !
  allocate(vStoik1(0:nEl))  ;  vStoik1(:)= 0
  allocate(vStoik2(0:nEl))  ;  vStoik2(:)= 0
  allocate(vStoik3(0:nEl))  ;  vStoik3(:)= 0
  !
  vStoik1(0:nEl)= vSpc(MixModel%vIPole(1))%vStoikio(0:nEl)
  vStoik2(0:nEl)= vSpc(MixModel%vIPole(2))%vStoikio(0:nEl)
  if(Ternary) &
  & vStoik3(0:nEl)= vSpc(MixModel%vIPole(3))%vStoikio(0:nEl)
  !
  if(Ternary) then  ;  i0= 2
  else              ;  i0= 1
  end if
  !
  N=0
  i=0
  do
    i= i+1
    k= 0
    do
      if(Ternary) k= k +1
      j= Dims -i -k
      N= N+1
      !----------------------------------------------------------------!
      !
      vParam(N)%I= I
      vParam(N)%J= J
      vParam(N)%K= K
      !
      ntot= real(I+J+K)
      P%vLPole(1)=I>0   ; P%vXPole(1)= I /ntot
      P%vLPole(2)=J>0   ; P%vXPole(2)= J /ntot
      P%vLPole(3)=K>0   ; P%vXPole(3)= K /ntot
      !
      call MixPhase_NormCompo(P)
      !
      call Species_Zero(S)
      !
      S%Typ= vSpc(MixModel%vIPole(1))%Typ
      S%WeitKg= dot_product(P%vXPole(1:MixModel%NPole), &
      &         vSpc(MixModel%vIPole(1:MixModel%NPole))%WeitKg)
      !
      if(Ternary) then
        call MixtureName( &
        & vSpc(MixModel%vIPole(1))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(2))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(3))%NamSp(1:3), &
        & I,J,K, &
        & DiscretName)
      else
        call MixtureName( &
        & vSpc(MixModel%vIPole(1))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(2))%NamSp(1:3), &
        & "xxx", &
        & I,J,K, &
        & DiscretName)
      end if
      !
      S%NamSp= trim(DiscretName)
      !
      S%vStoikio(0:nEl)= I*vStoik1(0:nEl) + J*vStoik2(0:nEl) + K*vStoik3(0:nEl)
      !
      call DiscretModel_Formula_Build( &
      & vEle,S, &
      & Formul)
      S%Formula= trim(Formul)
      !
      vSpcOut(N)= S
      !
      if(iDebug>0) write(fTrc,'(2(A,A1))') &
      & trim(DiscretName),T_,trim(Formul),T_
      !
      !----------------------------------------------------------------!
      if(.not. Ternary) exit
      if(k==Dims-1 -i) exit
    enddo
    !pause_
    if(i==Dims -i0) exit
  enddo
  !
  deallocate(vStoik1,vStoik2,vStoik3)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretModel_Init"
  !
  return
end subroutine DiscretModel_Init

subroutine MixtureName(S1,S2,S3,I,J,K,S)
  character(len=3),intent(in):: S1,S2,S3
  integer,         intent(in) ::I,J,K
  character(len=*),intent(out)::S
  !
  character(len=3):: Str
  !
  S=""
  write(Str,'(I2)') I  ; if(I<10) Str(1:1)='0'
  S=trim(S)//trim(S1)//"-"//trim(adjustl(Str))
  !
  S=trim(S)//"_"
  write(Str,'(I2)') J  ; if(J<10) Str(1:1)='0'
  S=trim(S)//trim(S2)//"-"//trim(adjustl(Str))
  !
  if(K==0) return
  !
  S=trim(S)//"_"
  write(Str,'(I2)') K  ; if(K<10) Str(1:1)='0'
  S=trim(S)//trim(S3)//"-"//trim(adjustl(Str))
  !
  return
end subroutine MixtureName

subroutine DiscretSpecies_Stoikio_Calc( &
& vEle,          & !IN
& vMixModel,     & !IN
& vDiscretModel, & !IN
& vDiscretParam, & !IN
& vSpc)            !INOUT with updated stoichio's for discretized species
!--
!-- update the stoichiometries
!-- of all discretized mixture phases
!-- according to models read from vDiscretModel,
!--
  use M_T_Element, only: T_Element
  use M_T_Species, only: T_Species
  use M_T_MixPhase,only: T_MixPhase
  use M_T_MixModel,only: T_MixModel, MixModel_Index
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  !
  type(T_Element),     intent(in):: vEle(:)
  type(T_MixModel),    intent(in):: vMixModel(:)
  type(T_DiscretModel),intent(in):: vDiscretModel(:)
  type(T_DiscretParam),intent(in):: vDiscretParam(:)
  !
  type(T_Species),  intent(inout):: vSpc(:)
  !
  type(T_DiscretModel):: DsModl
  type(T_MixModel):: MxModl
  ! type(T_MixPhase):: MxPhas
  !
  logical :: Ternary
  integer :: I,J,K
  integer :: iSp,iMix0,iMix,iDis,nEl
  !
  integer,allocatable:: vStoik1(:),vStoik2(:),vStoik3(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretSpecies_Stoikio_Calc"
  !
  nEl= size(vEle)
  !
  iMix0= 0
  !
  allocate(vStoik1(0:nEl))  ;  vStoik1(:)= 0
  allocate(vStoik2(0:nEl))  ;  vStoik2(:)= 0
  allocate(vStoik3(0:nEl))  ;  vStoik3(:)= 0
  !
  do iSp=1,size(vSpc)
    !
    if(vSpc(iSp)%iDiscret==0) cycle
    !
    iDis=    vSpc(iSp)%iDiscret
    DsModl=  vDiscretModel(vDiscretParam(iDis)%iModel)
    iMix=    DsModl%iMix
    MxModl=  vMixModel(iMix)
    !nDiscret= DsModl%Dim1
    !
    if(iMix/=iMix0) then !> change parameters only when needed
      !
      Ternary= (MxModl%NPole==3)
      !
      vStoik1(0:nEl)= vSpc(MxModl%vIPole(1))%vStoikio(0:nEl)
      vStoik2(0:nEl)= vSpc(MxModl%vIPole(2))%vStoikio(0:nEl)
      if(Ternary) &
      & vStoik3(0:nEl)= vSpc(MxModl%vIPole(3))%vStoikio(0:nEl)
      !
      if(iDebug>0) &
      write(fTrc,'(4(2A,/))') &
      & "disc.model=", DsModl%Name, &
      & "mixt.model=", MxModl%Name, &
      & "pole 1=    ", vSpc(MxModl%vIPole(1))%Formula, &
      & "pole 2=    ", vSpc(MxModl%vIPole(2))%Formula
      !
      iMix0= iMix
    end if
    !
    I= vDiscretParam(iDis)%I
    J= vDiscretParam(iDis)%J
    K= vDiscretParam(iDis)%K
    !
    vSpc(iSp)%vStoikio(0:nEl)= &
    & I*vStoik1(0:nEl) + J*vStoik2(0:nEl) + K*vStoik3(0:nEl)
    !
    !~ call DiscretModel_Formula_Build( &
    !~ & vEle,Dims-2,vStoik1(0),vStoik2(0),I,J, &
    !~ & vStoik1(1:nEl),vStoik2(1:nEl), &
    !~ & Formul)
    !~ S%Formula= trim(Formul)
    !
  enddo
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretSpecies_Stoikio_Calc"
  !
  return
end subroutine DiscretSpecies_Stoikio_Calc

subroutine DiscretParam_Init( &
& vEle,vSpc,vMixModel,vDiscretModel, &
& vDiscretParam,vSpcOut)
!--
!--initialize vDiscretParam
!--
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  use M_T_Element, only: T_Element
  use M_T_Species, only: T_Species
  use M_T_MixModel,only: T_MixModel
  !
  type(T_Element),     intent(in):: vEle(:)
  type(T_Species),     intent(in):: vSpc(:)
  type(T_MixModel),    intent(in):: vMixModel(:)
  type(T_DiscretModel),intent(in):: vDiscretModel(:)
  !
  type(T_DiscretParam),intent(out):: vDiscretParam(:)
  type(T_Species),     intent(out):: vSpcOut(:)
  !
  type(T_Species),     allocatable:: vSpcTmp(:)
  type(T_DiscretParam),allocatable:: vParam(:)
  integer:: I, iM, iMix, dimDis, N
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretParam_Init"
  !
  N= 0
  do iM=1,size(vDiscretModel)
    !
    dimDis= vDiscretModel(iM)%Dim1
    dimDis= vDiscretModel(iM)%DimTot
    !
    allocate(vSpcTmp(1:dimDis))
    allocate(vParam(1:dimDis))
    !
    iMix= vDiscretModel(iM)%iMix
    !
    call DiscretModel_Init( &
    & vEle,vSpc,vMixModel(iMix),vDiscretModel(iM), &
    & vParam,vSpcTmp)
    !
    do I=1,dimDis
      N= N+1
      !
      vSpcOut(N)= vSpcTmp(i)
      vSpcOut(N)%iDiscret= N
      !
      vDiscretParam(N)= vParam(i)  !-> for I,J,K
      vDiscretParam(N)%iModel= iM  !must come after
      !
      if(iDebug>0) then
        write(fTrc,'(2A)') vSpcOut(N)%NamSp, vSpcOut(N)%Typ
      end if
      !
    enddo
    !
    deallocate(vSpcTmp)
    deallocate(vParam)
    !
  enddo
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretParam_Init"
  !
end subroutine DiscretParam_Init

subroutine DiscretSpecies_TP_Update( &
!--
!--.update the T,P dependent parameters
!--.of all mixture phases discretized according
!--.to models read from vDiscretModel,
!--
& vMixModel,     & !IN
& TdgK,Pbar,     & !IN
& vDiscretModel, & !IN
& vDiscretParam, & !IN
& vSpc)            !INOUT with updated parameters for discretized species
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  use M_Dtb_Const, only: T_CK,R_jk
  use M_T_Species, only: T_Species
  use M_T_MixPhase,only: T_MixPhase,MixPhase_CalcMixing,MixPhase_NormCompo
  use M_T_MixModel,only: T_MixModel, MixModel_Margul_Wg_Calc, Mix_Site
  use M_T_MixModel,only: MixModel_Index,MixModel_XPoleToXSite
  !
  type(T_MixModel),    intent(in):: vMixModel(:)
  real(dp),            intent(in):: TdgK,Pbar
  type(T_DiscretModel),intent(in):: vDiscretModel(:)
  type(T_DiscretParam),intent(in):: vDiscretParam(:)
  !
  type(T_Species),  intent(inout):: vSpc(:)
  !
  type(T_DiscretModel):: DsModl
  type(T_MixModel):: MxModl
  type(T_MixPhase):: MxPhas
  real(dp),allocatable:: vXAtom(:)
  integer :: nDiscret
  real(dp):: G_IdMix,GMix,G_XsMix,GMecaRT,ntot !,GTotRT !,V0,WeitKg
  integer :: I,J,K
  integer :: iSp,iMix0,iMix,iDis
  !
  !character(len=30):: cFormat
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretSpecies_TP_Update"
  !
  iMix0= 0
  !
  do iSp=1,size(vSpc)
    !
    if(vSpc(iSp)%iDiscret==0) cycle
    !
    iDis=    vSpc(iSp)%iDiscret
    DsModl=  vDiscretModel(vDiscretParam(iDis)%iModel)
    !
    iMix=    DsModl%iMix
    MxModl=  vMixModel(iMix)
    nDiscret= DsModl%Dim1
    !
    !! print *,vSpc(iSp)%NamSp,iDis," ",DsModl%Name,MxModl%Name
    !
    MxPhas%iModel=    iMix
    MxPhas%vLPole=    .false.
    MxPhas%vXPole=    Zero
    !
    !---------------------check that parameters are not already uptodate
    if(iMix/=iMix0) then
      !--------------update Margules parameters for the current T,MxPhas
      if(MxModl%NMarg>0) &
      & call MixModel_Margul_Wg_Calc(TdgK,Pbar,MxModl)
      !----------------------------------------------------------------/
      iMix0= iMix
      !
      if(iDebug>2) then
        write(fTrc,'(2(2A,/))') &
        & "disc.model= ", trim(DsModl%Name), &
        & "mixt.model= ", trim(MxModl%Name)
        write(fTrc,'(2(2A,1X,G15.6,/))') &
        & "pole 1    = ", trim(vSpc(MxModl%vIPole(1))%Formula),vSpc(MxModl%vIPole(1))%G0rt, &
        & "pole 2    = ", trim(vSpc(MxModl%vIPole(2))%Formula),vSpc(MxModl%vIPole(2))%G0rt
      end if
      !
    end if
    !------------------------------------------------------------------/
    !
    I= vDiscretParam(iDis)%I
    J= vDiscretParam(iDis)%J
    K= vDiscretParam(iDis)%K
    !
    ntot= real(I+J+K)
    MxPhas%vLPole(1)=I>0 ; MxPhas%vXPole(1)= I /ntot
    MxPhas%vLPole(2)=J>0 ; MxPhas%vXPole(2)= J /ntot
    MxPhas%vLPole(3)=K>0 ; MxPhas%vXPole(3)= K /ntot
    !
    if(MxModl%Model==Mix_Site) then
      allocate(vXAtom(MxModl%NAtom))
      call MixModel_XPoleToXSite(  &
      & MxModl,MxPhas%vXPole(1:MxModl%NPole), &
      & vXAtom)
    end if
    !
    !write(cFormat,'(A,I3,A)') '(A,1X',MxModl%NPole,'(G12.3,1X))'
    !if(iDebug>2) &
    !& write(fTrc,cFormat) "MxModl%vIPole=",(MxPhas%vXPole(I),I=1,MxModl%NPole)
    !
    !~ M= 4 !size(vSpc(iSp)%vStoikio)
    !~ write(cFormat,'(A,I3,A)') '(A,1X,I6,1X,',M+1,'(I6,1X))'
    !~ if(iDebug>2) then
      !~ write(91,cFormat) "Pol1%vStoik=",&
      !~ & vSpc(MxModl%vIPole(1))%Div, &
      !~ & (vSpc(MxModl%vIPole(1))%vStoikio(I),I=0,M)

      !~ write(91,cFormat) "Pol2%vStoik=",&
      !~ & vSpc(MxModl%vIPole(2))%Div, &
      !~ & (vSpc(MxModl%vIPole(2))%vStoikio(I),I=0,M)

      !~ write(91,cFormat) "Spc%vStoik= ", &
      !~ & vSpc(iSp)%Div,(vSpc(iSp)%vStoikio(I),I=0,M)

      !~ write(91,*)
    !~ end if
    !
    if(iDebug>2) write(fTrc,'(A,3(F7.3,1X))',advance="NO") &
    & "I,J,K=", I/ntot, J/ntot, K/ntot
    !
    !!call MixPhase_NormCompo(MxPhas)
    !
    vSpc(iSp)%V0= dot_product( &
    & MxPhas%vXPole(1:MxModl%NPole), &
    & vSpc(MxModl%vIPole(1:MxModl%NPole))%V0)
    !
    call MixPhase_CalcMixing( &
    & TdgK,Pbar,MxModl,MxPhas, & !IN
    & GMix,G_IdMix,G_XsMix) !OUT
    ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
    !
    GMecaRT= dot_product( &
    & MxPhas%vXPole(1:MxModl%NPole), &
    & vSpc(MxModl%vIPole(1:MxModl%NPole))%G0rt)
    !
    vSpc(iSp)%G0rt= GMecaRT + GMix/R_jk/TdgK
    !
    if(iDebug>2) write(fTrc,'(4(G15.6,1X))') &
    & G_XsMix/R_jk/TdgK, &
    & GMix/R_jk/TdgK, &
    & GMecaRT, &
    & vSpc(iSp)%G0rt
    !
    if(allocated(vXAtom)) deallocate(vXAtom)
    !
  enddo
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretSpecies_TP_Update"
  !
end subroutine DiscretSpecies_TP_Update

subroutine DiscretModel_Formula_Build( &
& vEle,Spc, &
& S)
!--
!-- build the stoikio formula  of "pseudophases" in terms of elements
!--
  use M_T_Element,only: T_Element
  use M_T_Species,only: T_Species
  !
  type(T_Element), intent(in) :: vEle(:)
  type(T_Species), intent(in) :: Spc
  character(len=*),intent(out):: S
  !
  character(len=5):: Str
  integer:: iEl, I, nDiv
  !
  S=""
  nDiv= Spc%vStoikio(0)
  !
  do iEl=1,size(vEle)

    I= Spc%vStoikio(iEl)
    if(I>0) then
      write(Str,'(I5)') I
      S= trim(S)//vEle(iEl)%NamEl(1:2)//"("//trim(adjustl(Str))//")"
    end if

  enddo
  !
  write(Str,'(I3)') nDiv; if(nDiv<10) Str(1:1)='0'
  !
  S= trim(S)//"/("//trim(adjustl(Str))//")"
  !
end subroutine DiscretModel_Formula_Build

subroutine FelsparName(I,J,K,S)
  integer,         intent(in) ::I,J,K
  character(len=*),intent(out)::S
  !
  integer         ::I_,J_,K_
  character(len=2)::Str
  !
  S=""
  I_=I; J_=J; K_=K
  write(Str,'(I2)') I_; if(I_<10) Str(1:1)='0'; S=trim(S)//"AB"//Str//"_" !if(I_>0) I_=3*I_+1;
  write(Str,'(I2)') J_; if(J_<10) Str(1:1)='0'; S=trim(S)//"OR"//Str//"_" !if(J_>0) J_=3*J_+1;
  write(Str,'(I2)') K_; if(K_<10) Str(1:1)='0'; S=trim(S)//"AN"//Str//"_" !if(K_>0) K_=3*K_+1;
  !
end subroutine FelsparName

end module M_DiscretModel_Tools

!~ subroutine DiscretModel_Formula_Build_( &
!~ & vEle,nDiscret,nDiv1,nDiv2,I,J, &
!~ & vStoik1,vStoik2, &
!~ & S)
!~ !--
!~ !-- build the stoikio formula  of "pseudophases" in terms of elements
!~ !--
  !~ use M_T_Element,only: T_Element
  !~ !
  !~ type(T_Element), intent(in) :: vEle(:)
  !~ integer,         intent(in) :: nDiscret,nDiv1,nDiv2,I,J
  !~ integer,         intent(in) :: vStoik1(:),vStoik2(:)
  !~ character(len=*),intent(out):: S
  !~ !
  !~ integer,allocatable:: vStoik(:)
  !~ character(len=3):: Str
  !~ integer:: iEl
  !~ !
  !~ allocate(vStoik(1:size(vEle)))
  !~ !
  !~ S=""
  !~ do iEl=1,size(vEle)
    !~ !
    !~ vStoik(iEl)= I*vStoik1(iEl)*nDiv1 + J*vStoik2(iEl)*nDiv2
    !~ if(vStoik(iEl)>0) then
      !~ write(Str,'(I3)') vStoik(iEl)
      !~ S= trim(S)//vEle(iEl)%NamEl(1:2)//"("//trim(adjustl(Str))//")"
      !~ !!if(iDebug>0) write(fTrc,'(A)') trim(S)
    !~ end if
    !~ !
  !~ enddo
  !~ !
  !~ write(Str,'(I2)') nDiscret+2; if(nDiscret<10) Str(1:1)='0'
  !~ S= trim(S)//"/("//trim(adjustl(Str))//")"
  !~ !
  !~ deallocate(vStoik)
  !~ !
!~ end subroutine DiscretModel_Formula_Build_

!~ subroutine FelsparDiscretize
  !~ use M_Dtb_Const,   only: T_CK,R_jk
  !~ use M_T_Species,   only: T_Species
  !~ use M_Global_Vars, only: vEle,vSpc,vMixModel,tFormula
  !~ use M_Path_Vars,   only: tGrt,vTPpath
  !~ use M_T_MixModel,  only: T_MixModel,MixModel_Margul_Wg_Calc
  !~ use M_T_MixPhase
  !~ !
  !~ integer::I,J,K,N
  !~ !! type(T_Species),dimension(:),allocatable::vSpc0 !store vSpc
  !~ !
  !~ type(T_MixModel)::S
  !~ type(T_MixPhase)::P
  !~ real(dp):: Pbar,TdgK
  !~ real(dp):: GIdMix,GMix,GXsMix,GMeca,GTot
  !~ integer::nEl,iTP,iEl
  !~ character(len=15)::SolName
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A)') "FelsparDiscretize"
  !~ !
  !~ N=0
  !~ do I=0,Discret
    !~ do J=0,Discret-I
      !~ N=N+1
      !~ K=Discret-I-J
    !~ enddo
  !~ enddo
  !~ write(*,*) "N_tot=",N  !!; call Pause_
  !~ !
  !~ iTP=1
  !~ nEl=size(vEle)
  !~ P%Name= "FELSPAR"
  !~ P%iModel= 1
  !~ !
  !~ S=vMixModel(1)
  !~ P%vLPole=.false.; P%vXPole=0.0D0
  !~ !
  !~ do iTP=1,size(vTPpath)
    !~ Pbar= vTPpath(iTP)%Pbar
    !~ TdgK= vTPpath(iTP)%TdgC+T_CK
    !~ if(iDebug>0) write(fTrc,'(I3,2(1X,G12.3))') iTP,Pbar,TdgK
    !~ !
    !~ !update Margules parameters for the current T,P
    !~ if(S%NMarg>0) call MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
    !~ vMixModel(1)=S
    !~ !
    !~ N=0
    !~ do I=0,Discret
      !~ do J=0,Discret-I
        !~ N=N+1
        !~ K=Discret-I-J
        !~ P%vLPole(1)=I>0; P%vXPole(1)= real(I)/real(Discret)
        !~ P%vLPole(2)=J>0; P%vXPole(2)= real(J)/real(Discret)
        !~ P%vLPole(3)=K>0; P%vXPole(3)= real(K)/real(Discret)
        !~ !
        !~ call MixPhase_NormCompo(P)
        !~ call MixPhase_CalcMixing(TdgK,Pbar,S,P,GMix,GIdMix,GXsMix) !all OUT in joules
        !~ !GMeca=dot_product(P%vXPole(1:S%NPole),vSpc(S%vIPole(1:S%NPole))%tG0(iTP))*TdgK*R_jk
        !~ GMeca= dot_product(P%vXPole(1:S%NPole),tGrt(S%vIPole(1:S%NPole),iTP)) &
        !~ &    * TdgK*R_jk
        !~ GTot=GMeca+GMix
        !~ !
        !~ !! vGibFas(nFasPur+N,iTP)=GTot*real(Discret)
        !~ !
        !~ if(iTP==1) then !do only for first TP value
          !~ if(iDebug>0) write(fTrc,'(A,A1,I3,A1,4(G15.6,A1))') &
          !~ & Solname,T_,N,T_,P%vXPole(1),T_,P%vXPole(2),T_,P%vXPole(3),T_,GTot,T_
          !~ !
          !~ call FelsparName(I,J,K,SolName)
          !~ !! v%namSol(nFasPur+N)=trim(SolName)
          !~ !
          !~ !calculate stoikio of "pseudophases" in terms of elements
          !~ do iEl=1,nEl
            !~ tFormula(nFasPur+N,iEl)= &
            !~ & I*tFormula(S%vIPole(1),iEl) &
            !~ + J*tFormula(S%vIPole(2),iEl) &
            !~ + K*tFormula(S%vIPole(3),iEl)
            !~ !! write(*,'(A3,I6)') vNamEle(iEl),tFormula(nFasPur+N,iEl)
          !~ enddo
          !~ !call Pause_
        !~ end if
      !~ enddo
    !~ enddo
  !~ enddo
  !~ write(fTrc,'(/,A,I3,/)') "N_tot=",N
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A)') "FelsparDiscretize"
  !~ !! nFasSol=N
!~ end subroutine FelsparDiscretize

!~ subroutine DiscretModel_TP_Update( &
!~ !--
!~ !-- §§§!!! OBSOLETE !!!§§§!!! OBSOLETE !!!§§§!!! OBSOLETE !!!§§§!!!
!~ !-- for a mixture phase discretized according to model DiscretModel,
!~ !-- update the T,P dependent parameters
!~ !--
!~ & vSpc,        & !IN
!~ & vMixModel,   & !IN
!~ & TdgK,Pbar,   & !IN
!~ & DiscretModel, & !IN
!~ & vSpcDiscret)    !INOUT discretized species with updated parameters
  !~ use M_T_DiscretModel,only: T_DiscretModel !,T_DiscretParam
  !~ use M_Dtb_Const, only: T_CK,R_jk
  !~ use M_T_Species, only: T_Species
  !~ use M_T_MixPhase,only: T_MixPhase,MixPhase_CalcMixing,MixPhase_NormCompo
  !~ use M_T_MixModel,only: MixModel_Margul_Wg_Calc, T_MixModel
  !~ use M_T_MixModel,only: MixModel_Index,MixModel_XPoleToXSite
  !~ !
  !~ type(T_Species),     intent(in):: vSpc(:)
  !~ type(T_MixModel),    intent(in):: vMixModel(:)
  !~ real(dp),            intent(in):: TdgK,Pbar
  !~ type(T_DiscretModel),intent(in):: DiscretModel
  !~ !
  !~ type(T_Species),  intent(inout):: vSpcDiscret(:)
  !~ !
  !~ type(T_MixModel):: MxModl
  !~ type(T_MixPhase):: MxPhas
  !~ integer :: nDiscret
  !~ real(dp):: G_IdMix,GMix,G_XsMix,GMecaRT !,GTotRT !,V0,WeitKg
  !~ integer :: N,I,J,K
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_TP_Update"
  !~ !
  !~ MxModl=   vMixModel(DiscretModel%iMix)
  !~ nDiscret= DiscretModel%DimTot
  !~ !
  !~ MxPhas%iModel=    DiscretModel%iMix
  !~ MxPhas%vLPole=    .false.
  !~ MxPhas%vXPole=    Zero
  !~ !
  !~ !--- compute stoikio vectors of end members
  !~ if(iDebug>0) write(fTrc,'(A)') vSpc(MxModl%vIPole(1))%Formula !debug_
  !~ if(iDebug>0) write(fTrc,'(A)') vSpc(MxModl%vIPole(2))%Formula
  !~ !
  !~ !--- update Margules parameters for the current T,MxPhas
  !~ if(MxModl%NMarg>0) &
  !~ & call MixModel_Margul_Wg_Calc(TdgK,Pbar,MxModl)
  !~ !
  !~ N= 0
  !~ K= 0
  !~ !!do I=0,nDiscret-1
  !~ do I=1,nDiscret
    !~ N=N+1 !-> for binary solution
    !~ !do J=0,Discret-I !-> for ternary solution
      !~ !!J=nDiscret-1 -I
      !~ J= nDiscret+1 -I
      !~ !!N=N+1 !-> for ternary solution
      !~ !!K=Discret-I-J !-> for ternary solution
      !~ !
      !~ MxPhas%vLPole(1)= I>0  ; MxPhas%vXPole(1)= real(I) /real(nDiscret+1)
      !~ MxPhas%vLPole(2)= J>0  ; MxPhas%vXPole(2)= real(J) /real(nDiscret+1)
      !~ MxPhas%vLPole(3)= K>0  ; MxPhas%vXPole(3)= real(K) /real(nDiscret+1)
      !~ !
      !~ call MixPhase_NormCompo(MxPhas)
      !~ !
      !~ !if(trim(MxModl%Model)=="SITE") then
      !~ !  allocate(vXAtom(MxModl%NAtom))
      !~ !  call MixModel_XPoleToXSite(  &
      !~ !  & MxModl,MxPhas%vXPole(1:MxModl%NPole), &
      !~ !  & vXAtom)
      !~ !end if
      !~ !
      !~ vSpcDiscret(N)%V0= &
      !~ & dot_product(MxPhas%vXPole(1:MxModl%NPole),vSpc(MxModl%vIPole(1:MxModl%NPole))%V0)
      !~ !
      !~ call MixPhase_CalcMixing( &
      !~ & TdgK,Pbar,MxModl,MxPhas, &  !IN
      !~ & GMix,G_IdMix,G_XsMix) !OUT
      !~ ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
      !~ !
      !~ GMecaRT= dot_product( &
      !~ & MxPhas%vXPole(1:MxModl%NPole), &
      !~ & vSpc(MxModl%vIPole(1:MxModl%NPole))%G0rt)
      !~ !
      !~ vSpcDiscret(N)%G0rt= GMecaRT + GMix/R_jk/TdgK
      !~ !
    !~ !enddo
  !~ enddo
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretModel_TP_Update"
  !~ !
!~ end subroutine DiscretModel_TP_Update
