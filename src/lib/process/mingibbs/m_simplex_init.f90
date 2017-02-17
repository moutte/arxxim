module M_Simplex_Build
!--
!-- io module for simplex calculation
!-- independant (provisionally) from other io modules
!--
  use M_Kinds
  use M_Trace,      only: fTrc,iDebug,T_,Stop_,Pause_,Warning_
  use M_T_Component,only: T_Component, Component_Zero
  !
  implicit none
  !
  private
  !
  public:: Simplex_Build
  public:: Simplex_CloseFiles
  !
  integer:: nFasPur,nFasSol
  integer:: Discret
  != dimension of discretization (=number of points on each parameter)
  !
  !! logical:: Redox
  !
  integer,public:: fSpl1=0, fSpl2=0
  ! files that should remain accessible
  ! by different subroutines throughout program
  !
contains

subroutine Simplex_Build
  !
  use M_T_Element,   only: Element_Index
  use M_T_Species,   only: T_Species,Species_Stoikio_Calc,Species_Append
  use M_Dtb_Const,   only: Tref,Pref,T_CK
  !
  use M_Global_Alloc,only: MixModels_Alloc,Phases_Alloc,MixPhases_Alloc
  use M_Global_Alloc,only: DiscretParam_Alloc
  use M_Global_Tools,only: Global_TP_Update
  !
  use M_DiscretModel_Read
  use M_DiscretModel_Tools !DiscretParam_Init
  !
  use M_Global_Vars, only: vEle,vSpcDtb,vSpc,vMixModel,vMixFas,vFas
  use M_Global_Vars, only: vDiscretModel,vDiscretParam
  use M_Global_Vars, only: tFormula
  !
  use M_GEM_Vars,    only: TdgK,Pbar,vCpnGEM,tStoikioGEM
  !---------------------------------------------------------------------
  logical:: Ok
  integer:: N
  type(T_Species),allocatable:: vSpcTmp(:)
  !---------------------------------------------------------------------
  if(iDebug>2) write(fTrc,'(/,A)') "<---------------------Simplex_Build"
  !
  nFasPur=0
  nFasSol=0
  !
  TdgK= Tref
  Pbar= Pref
  !
  ! if(Element_Index("O__",vEle)==0) &
  ! & call Stop_("Element Oxygen not found in vEle")
  !
  !--build vCpnGEM
  call Component_Read(vEle)
  !
  !--rebuild species array vSpc consistent with vCpnGEM
  call Species_Select(size(vEle),vCpnGEM)
  !
  !--build new vEle
  call Element_Select(vCpnGEM)
  !
  call Element_Sort(vEle)
  !
  call Component_Stoikio_Calc(vEle,vCpnGEM,Ok)
  !
  !--build species array vSpc consistent with vEle
  !--OBSOLETE
  ! call Spl_Species_Alloc(vEle)
  !
  call Species_Stoikio_Calc(vEle,vSpc,Ok)
  !
  !--build MixModel database -> vMixModel
  call MixModels_Alloc(vSpc, vMixModel)
  !
  !-------------------------------------------append discretized species
  !--read discretization models, build vDiscretModel
  call DiscretModel_Read(vMixModel)
  !
  N= size(vDiscretModel)
  !
  if(N>0) then
    !
    call DiscretParam_Alloc(vDiscretModel,  vDiscretParam)
    !
    allocate(vSpcTmp(size(vDiscretParam)))
    !
    call DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    call Species_Append(vSpcTmp, vSpc) !-> new vSpc !!!
    !
    deallocate(vSpcTmp)
    !
    call DiscretSpecies_Stoikio_Calc( & !
    & vEle,          & !IN
    & vMixModel,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
    !
  end if
  !------------------------------------------/append discretized species
  !
  !--- read phase compositions, build vMixFas
  call MixPhases_Alloc(vSpc,vMixModel,   vMixFas)
  !
  !--- build vFas
  call Phases_Alloc(vSpc,vMixFas,    vFas)
  !
  if(iDebug>2) print *,"TdgC/Pbar=", TdgK-T_CK,Pbar  !;  pause
  !
  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  !
  !--- compute tFormula
  call Spl_Tables_Build(vEle,vSpc,vFas)
  !
  !--- allocate, compute tStoikioGEM
  if(allocated(tStoikioGEM)) deallocate(tStoikioGEM)
  allocate(tStoikioGEM(size(vFas),size(vCpnGEM)))
  !
  call Spl_Stoikio_Calc( &
  & vEle, vFas,          &
  & tFormula,            &
  & vCpnGEM,             &
  & tStoikioGEM)
  !---/
  !
  !--- write system stoikiometry
  if(iDebug>2) call Spl_WriteSystem(vFas,tFormula,vCpnGEM,tStoikioGEM)
  !
  if(iDebug>2) write(fTrc,'(A,/)') "</--------------------Simplex_Build"
  !
  return
end subroutine Simplex_Build

subroutine Species_Select(nEl,vCpn)
!--
!-- rebuild species array vSpc consistent with vCpnGEM
!--
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  !
  use M_Global_Vars, only: vSpc
  !
  integer, intent(in):: nEl
  type(T_Component),intent(inout):: vCpn(:)
  !
  type(T_Species),allocatable:: vTmp(:)
  real(dp),allocatable:: tStoikio(:,:) !tStoikio(:,:),
  integer :: I,iS,nS,nC,N
  logical :: Ok

  if(iDebug>2) write(fTrc,'(/,A)') "< Species_Select"
  
  nC= size(vCpn)
  nS= size(vSpc)
  allocate(vTmp(nS))
  vTmp(:)= vSpc(:)

  allocate(tStoikio(nC+1,nEl))
  do I=1,nC
    tStoikio(I,1:nEl)= vCpn(I)%vStoikCp(1:nEl)  !tStoikio(I,:)
  end do

  n= 0
  do iS=1,nS
  
    if(vSpc(iS)%Typ /= "AQU" .or. trim(vSpc(iS)%NamSp)=="H2O") then
    
      tStoikio(nC+1,1:nEl)= vSpc(iS)%vStoikio(1:nEl)
      
      call Check_Independent(tStoikio,Ok)
      
      if(.not. Ok) then
        ! if species is dependent on basis vCpn(:),
        ! then add to vSpc
        N= N+1
        vTmp(N)= vSpc(iS)
        !~ write(6,'(2A)') "OK SPECIES= ", trim(vSpc(is)%NamSp)
        if(iDebug>2) write(fTrc,'(A)') trim(vSpc(iS)%NamSp)
      end if
      
    end if
  
  end do

  deallocate(vSpc)
  allocate(vSpc(N))
  vSpc(1:N)= vTmp(1:N)

  deallocate(vTmp)
  deallocate(tStoikio)

  if(iDebug>2) write(fTrc,'(A,/)') "</ Species_Select"
  
  return
end subroutine Species_Select

subroutine Check_Independent(tStoikio,Ok)
  use M_Numeric_Tools,only: iMinLoc_R,iFirstLoc,iMaxLoc_R
  !
  real(dp),intent(in) :: tStoikio(:,:)
  logical, intent(out):: Ok
  !
  real(dp),allocatable:: tStoikTmp(:,:) !tStoikio(:,:),
  integer, allocatable:: vIPivot(:)
  !
  integer :: I,J,K
  integer :: N,nEl
  real(dp):: Y,Pivot
  
  N=   size(tStoikio,1)
  nEl= size(tStoikio,2)
  allocate(tStoikTmp(N,nEl))     ;  tStoikTmp(:,:)= Zero
  allocate(vIPivot(nEl))         ;  vIPivot(:)= 0
  !
  Ok= .true.
  K= 0
  do I=1,N
    !
    K= K+1
    tStoikTmp(K,:)= tStoikio(I,:)
    do J=1,K-1
      Y= tStoikTmp(K,vIPivot(J))
      tStoikTmp(K,:)= tStoikTmp(K,:)- Y*tStoikTmp(J,:)
    end do
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(K,:)))
    Pivot= tStoikTmp(K,vIPivot(K))
    if(ABS(Pivot)<1.D-9) then
      != all coeff's are 0 -> this species is not independent -> exit ==
      Ok= .false.
      deallocate(tStoikTmp,vIPivot) !tStoikio,
      return
    end if
    !
    tStoikTmp(K,:)= tStoikTmp(K,:) /Pivot
    !
  end do
  !
  ! write(sFMT,'(a,i3,a)') '(',nEl,'(G12.3,1X))'
  ! do I=1,N
  !   write(12,sFMT) (tStoikTmp(I,K),K=1,nEl)
  ! end do
  !
  deallocate(tStoikTmp,vIPivot) !tStoikio,
  !
end subroutine Check_Independent

subroutine Simplex_CloseFiles
  if(fSpl1>0) close(fSpl1) ; fSpl1= 0
  if(fSpl2>0) close(fSpl2) ; fSpl2= 0
end subroutine Simplex_CloseFiles

subroutine Component_Read(vEle)
!--
!-- read components from SYSTEM.GEM block -> build vCpnGEM
!--
  use M_IOTools,  only: Str_Append
  use M_Files,    only: NamFInn
  use M_Dtb_Const,only: T_CK
  use M_IoTools,  only: GetUnit,LinToWrd,AppendToEnd,WrdToReal,Str_Upper
  use M_T_Element,only: T_Element,Formula_Read,Element_Index,Formula_Build
  use M_Dtb_Read_Tools
  !
  use M_GEM_Vars,only: vCpnGEM,TdgK,Pbar
  !
  type(T_Element),intent(in):: vEle(:)
  !
  character(len=255)    :: L,W
  type(T_Component)     :: Cpn
  logical :: sOk,fOk,CpnOk
  integer :: nEl,N,i,iOx
  integer :: ZSp,nDiv,Ztest
  integer :: f,ios
  logical :: UnitIsMole,UnitIsGram
  real(dp):: X
  character(len=4):: Str
  !
  real(dp),allocatable:: tStoikio(:,:)
  integer, allocatable:: vStoik(:)
  type(T_Component),allocatable:: vC(:)
  !
  logical:: EcformIsOk
  integer:: fFormula
  character(len=2),dimension(:),allocatable:: vElement
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Component_Read"
  !
  !CodFormula= "ECFORM"
  fFormula= 0
  allocate(vElement(size(vEle)+3))
  call DtbRead_Build_vElement(vEle,vElement)
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  iOx= Element_Index("OX_",vEle)
  !
  nEl=size(vEle)
  allocate(vStoik(1:nEl))
  allocate(vC(1:nEl))
  !
  if(iDebug>2) print *,"Component_Build"
  !
  N=0
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,sOk)
    if(W(1:1)=='!') cycle DoFile
    call AppendToEnd(L,W,sOk)
    !
    select case(W)
    !
    case("SYSTEM.SIMPLEX","SYSTEM.GEM")
      !
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,sOk)
        if(W(1:1)=='!') cycle DoBlock
        call AppendToEnd(L,W,sOk)
        !
        UnitIsMole= .true.
        UnitIsGram= .false.
        !
        select case(W)
        
        case("ENDINPUT")                 ;  exit DoFile
        case("END","ENDSYSTEM.SIMPLEX")  ;  exit DoBlock
        case("ENDSYSTEM.GEM")            ;  exit DoBlock
        
        case("TDGC")
          call LinToWrd(L,W,sOk)
          call WrdToReal(W,X)
          TdgK= X +T_CK
          cycle DoBlock !----------------------------------cycle DoBlock
        
        case("TDGK")
          call LinToWrd(L,W,sOk)
          call WrdToReal(W,X)
          TdgK= X
          cycle DoBlock !----------------------------------cycle DoBlock
        
        case("PBAR")
          call LinToWrd(L,W,sOk)
          call WrdToReal(W,X)
          Pbar= X
          cycle DoBlock !----------------------------------cycle DoBlock
        
        case("MOLE")
          UnitIsMole= .true.
          call LinToWrd(L,W,sOk)
        
        case("GRAM")
          UnitIsGram= .true.
          call LinToWrd(L,W,sOk)
          
        end select
        !
        call Component_Zero(Cpn)
        !
        call Str_Append(W,3)
        Cpn%NamCp= trim(W)
        !
        call LinToWrd(L,W,sOk,"NO")
        !-> compact formula, character case should be conserved :!!
        !
        !if(CodFormula=="SCFORM") then
          call DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
          if(.not.EcformIsOk) then
            if(iDebug>0) print '(3A)',trim(W)," =ERROR IN FORMULA, SKIPPED !!"
            cycle DoBlock
          end if
        !end if
        !
        call Str_Upper(W)
        !
        Cpn%Formula= trim(W)
        !
        !----------------------------------------------- process formula 
        !--------------------- if formula is Ok, save component in vC(:) 
        call Formula_Read(Cpn%Formula,vEle,ZSp,nDiv,fOk,vStoik)
        !
        Ztest= dot_product(vStoik(1:nEl),vEle(1:nEl)%Z)
        ! if(fOk .and. ZSp==Ztest) then
        if(fOk) then
          !
          ! vStoikio has globally fixed size, nElMax,
          ! vStoik has local size nEl !!
          !
          Cpn%Statut= "INERT"
          Cpn%vStoikCp(1:nEl)= vStoik(1:nEl)
          Cpn%vStoikCp(0)=     nDiv
          Cpn%vStoikCp(nEl+1)= ZSp
          !
          if(Zsp/=Ztest) then
            if(iOx/=0) then
              Cpn%vStoikCp(iOx)= Zsp - Ztest
              write(Str,'(I4)') Cpn%vStoikCp(iOx)
              Cpn%Formula= trim(Cpn%Formula)//"OX("//trim(adjustl(Str))//")"
            else
              print *,"ZSp,Stoikio=",ZSp,Ztest
              if(iDebug>0) write(fTrc,'(A,A1,A15,A1,A39)') &
              & "REJECT",T_,Cpn%NamCp,T_,Cpn%Formula
              if(iDebug>0) print '(3A)',"rejected: ",Cpn%NamCp,Cpn%Formula
              cycle
            end if
          end if
          !
          !call Formula_Build(vEle,Cpn%vStoikCp,Zsp,nDiv,Cpn%Formula)
          !
          Cpn%iMix= 0
          !
          if(iDebug>0) write(fTrc,'(A,A1,A15,A1,A39,A1,I3)') &
          & "ACCEPT",T_,Cpn%NamCp,T_,Cpn%Formula,T_,N
          if(iDebug>2) print '(3A)',"accepted: ",Cpn%NamCp,Cpn%Formula
          !
          call LinToWrd(L,W,sOk)
          call WrdToReal(W,Cpn%Mole)
          if(UnitIsGram) then
            X= dot_product( &
            & Cpn%vStoikCp(1:nEl), &
            & vEle(1:nEl)%WeitKg)/Cpn%vStoikCp(0)
            Cpn%Factor= X *1.0D3
            Cpn%Mole= Cpn%Mole /Cpn%Factor
            !! print *,"Cpn%Factor",Cpn%Factor
          end if
          !
          N=N+1
          vC(N)= Cpn
          !
        else
          !
          print *,"ZSp,Stoikio=",ZSp,Ztest
          if(iDebug>0) write(fTrc,'(A,A1,A15,A1,A39)') &
          & "REJECT",T_,Cpn%NamCp,T_,Cpn%Formula
          if(iDebug>0) print '(3A)',"rejected: ",Cpn%NamCp,Cpn%Formula
          !
        end if
        !-------------------------------------------/ process formula --
        !
      end do DoBlock
      !
      if(N==0) call Stop_("< Found No Components...Stop") !--=== stop ==
      !
    end select
    !
  end do DoFile
  !
  !! if(iDebug>2) call Pause_
  close(f)
  !
  if(allocated(vCpnGEM)) deallocate(vCpnGEM)
  allocate(vCpnGEM(1:N))
  !
  vCpnGEM(1:N)= vC(1:N)
  !
  deallocate(vStoik)
  deallocate(vC)
  !
  !----------------------------------------------- check independency --
  allocate(tStoikio(N,nEl))
  do I=1,N
    tStoikio(I,1:nEl)= vCpnGEM(I)%vStoikCp(1:nEl)  !tStoikio(I,:)
  end do
  call Check_Independent(tStoikio,CpnOk)
  if(.not. CpnOk) call Stop_("Components Not Independent")
  deallocate(tStoikio)
  !-----------------------------------------------/check independency --
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Component_Read"
  !
  return
end subroutine Component_Read

subroutine Element_Select(vCpn)
!--
!-- reduce element list to those involved in vCpn,
!-- then re-build vEle
!--
  use M_T_Component,only: T_Component
  use M_T_Element,  only: T_Element
  !
  use M_Global_Vars,only: vEle
  !
  type(T_Component),intent(inout):: vCpn(:)
  !
  integer,allocatable:: tStoik(:,:)
  integer:: nEl, i, N
  type(T_Element),allocatable:: vE(:)
  !
  if(iDebug>2) write(fTrc,'(A)') "<----------------------Element_Select"
  if(iDebug>2) write(fTrc,'(A)') "-> New Element List"
  !
  nEl= size(vEle)
  !
  allocate(tStoik(size(vCpn),nEl))
  do i=1,size(vCpn)
    tStoik(i,1:nEl)= vCpn(i)%vStoikCp(1:nEl)
  end do
  !
  if(iDebug>2) then
    write(fTrc,'(*(A,A1))') (vEle(i)%NamEl,t_,i=1,nEl)
    call WriteMat_I(fTrc,tStoik)
  end if
  !
  allocate(vE(1:nEl))
  N=0
  !
  !---------------------------select elements involved in the components
  do I=1,nEl
    if(ANY(tStoik(:,i)/=0)) then
      N=N+1
      vE(N)= vEle(I) !-> new element list
      if(iDebug>2) write(fTrc,'(I3,1X,A)') N,vEle(I)%NamEl
    end if
  end do
  !---/
  !
  deallocate(vEle); allocate(vEle(1:N))
  vEle(1:N)= vE(1:N)
  !
  deallocate(tStoik)
  deallocate(vE)
  !
  if(iDebug>2) write(fTrc,'(A,/)') "</-------------------Element_Select"
  !
  return
end subroutine Element_Select

subroutine Element_Sort(vEle)
  use M_T_Element,only: T_Element,Element_Index
  !
  type(T_Element),intent(inout):: vEle(:)
  !
  type(T_Element):: E
  integer:: i
  !
  i= Element_Index("O__",vEle)
  if(i/=1) then
    E= vEle(1)
    vEle(1)= vEle(i)
    vEle(i)= E
  end if
  !
end subroutine Element_Sort

subroutine Component_Stoikio_Calc(vEle,vCpn,Ok)
!--
!-- update vCpn(:)%vStoikCp
!--
  use M_T_Element,  only: T_Element,Element_Index
  use M_T_Component,only: T_Component,Component_Stoikio
  !
  type(T_Element),  intent(in)   :: vEle(:)
  type(T_Component),intent(inout):: vCpn(:) !modifs in nDiv,Z,Oxy
  logical,          intent(out)  :: Ok
  !
  integer:: I,J,ieOx
  logical:: fOk
  !
  if(iDebug>2)  write(fTrc,'(A)') "< Component_Stoikio_Calc"
  !
  ieOx= Element_Index("OX_",vEle)
  !
  Ok= .true.
  do I=1,size(vCpn)
    call Component_Stoikio(vEle,ieOx,vCpn(I),fOk)
    if(.not. fOk) then
      Ok= .false.
      !! Msg= "Stoikiometry problem in species "//trim(vSpc(I)%NamSp)
      return
    end if
  end do
  !
  !------------------------------------------------------------ trace --
  if(iDebug>2) then
    do i=1,size(vEle)
      write(fTrc,'(A,A1)',advance="no") vEle(i)%NamEl,t_
    end do
    write(fTrc,*)
    do i=1,size(vCpn)
      do j=1,size(vEle)
        write(fTrc,'(I3,A1)',advance="no") vCpn(i)%vStoikCp(j),t_
      end do
      write(fTrc,*)
    end do
  end if
  !-----------------------------------------------------------/ trace --
  if(iDebug>2)  write(fTrc,'(A)') "</ Component_Stoikio_Calc"
  !
  return
end subroutine Component_Stoikio_Calc

integer function iMaxLoc_I(ARR) !"NR"
  real(dp),intent(in) :: ARR(:)
  integer:: IMAX(1)
  IMAX=MAXLOC(ARR(:))
  iMaxLoc_I=IMAX(1)
end function iMaxLoc_I

!old!!-----------------------------------------------------------------------
!old!!-- retrieve all species consistent with current vEle
!old!!-- OBSOLETE: replaced by Species_Select
!old!!-----------------------------------------------------------------------
!old!subroutine Spl_Species_Alloc(vEle)
!old!  use M_T_Element,   only: T_Element,Element_Index
!old!  use M_T_Species,   only: T_Species,Species_Stoikio
!old!  use M_Global_Vars, only: vSpc !,vSpcDtb
!old!  use M_Dtb_Calc,    only: DtbSpc_Table_Build
!old!  !
!old!  type(T_Element),intent(in):: vEle(:)
!old!  !
!old!  type(T_Species),allocatable:: vSpc_(:)
!old!  integer:: nEl,nSp,ieOx,I,N
!old!  logical:: fOk
!old!  !
!old!  if(iDebug>0) write(fTrc,'(/,A)') "< SplSpecies_Alloc"
!old!  !
!old!  nSp=  size(vSpc)
!old!  nEl=  size(vEle)
!old!  ieOx= Element_Index("OX_",vEle)
!old!  !
!old!  allocate(vSpc_(1:nSp))
!old!  !
!old!  !--- select species consistent with vEle
!old!  N=0
!old!  do I=1,nSp
!old!    if(vSpc(I)%Typ /= "AQU" .or. trim(vSpc(I)%NamSp)=="H2O") then
!old!      call Species_Stoikio(vEle,ieOx,vSpc(I),fOk)
!old!      if(fOk) then
!old!        N=N+1
!old!        vSpc_(N)= vSpc(I)
!old!      end if
!old!    end if
!old!  end do
!old!  !
!old!  deallocate(vSpc); allocate(vSpc(1:N))
!old!  vSpc(1:N)=   vSpc_(1:N)
!old!  deallocate(vSpc_)
!old!  !
!old!  !~ nSp=  size(vSpc)
!old!  !~ call Warning_("vTPpath must be allocated before Spl_Species_Alloc")
!old!  !~ nTP=  size(vTPpath)
!old!  !~ if(allocated(tGrt)) deallocate(tGrt); allocate(tGrt(1:nSp,1:nTP))
!old!  !~ call DtbSpc_Table_Write(vTPpath,vSpcDtb,vSpc,tGrt)
!old!  !~ !
!old!  if(iDebug>0) then
!old!    write(fTrc,'(A)') "-> New Species List"
!old!    do i=1,size(vSpc)
!old!      write(fTrc,'(I3,1X,A)') I,vSpc(I)%NamSp
!old!    end do
!old!  end if
!old!  !
!old!  if(iDebug>0) write(fTrc,'(A,/)') "</ SplSpecies_Alloc"
!old!  !
!old!end subroutine Spl_Species_Alloc

subroutine Spl_Tables_Build(vEle,vSpc,vFas)
!--
!-- build tFormula(1:nFs,1:nEl): table of species stoikios in lines
!--
  use M_T_Element,   only: T_Element,Formula_Read
  use M_T_Species,   only: T_Species
  use M_T_Phase,     only: T_Phase
  use M_Global_Vars, only: tFormula
  !
  type(T_Element), intent(in):: vEle(:)
  type(T_Species), intent(in):: vSpc(:)
  type(T_Phase),   intent(in):: vFas(:)
  !
  integer:: nEl,nFs,I
  integer,dimension(:),allocatable:: vStoik
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Spl_BuildTables"
  !
  nEl=size(vEle)
  nFs=size(vFas)
  !
  deallocate(tFormula)
  allocate(tFormula(1:nFs,1:nEl))
  !
  allocate(vStoik(0:nEl))
  !
  do I=1,nFs
    vStoik(0:nEl)= vSpc(vFas(I)%iSpc)%vStoikio(0:nEl)
    tFormula(I,1:nEl)= vStoik(1:nEl) /real(vStoik(0))
  end do
  !
  deallocate(vStoik)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Spl_BuildTables"
end subroutine Spl_Tables_Build

subroutine Spl_Stoikio_Calc( &
& vEle, vFas,                &
& tFormula,                  &
& vCpn,                      &
& tStoikio)
!--
!-- Compute tStoikioCpn from tStoikEle
!-- tStoikio(1:nFs,1:nCp)- stoikio of phases (line) in terms of components (col)
!--
  !---------------------------------------------------------------------
  use M_T_Component, only: T_Component
  use M_T_Element,   only: T_Element
  use M_T_Phase,     only: T_Phase
  use M_T_Element,   only: Formula_Read
  !
  use M_Numeric_Mat, only: LU_Decomp, LU_BakSub
  !---------------------------------------------------------------------
  type(T_Element),  intent(in) :: vEle(:)
  type(T_Phase),    intent(in) :: vFas(:)
  real(dp),         intent(in) :: tFormula(:,:)
  type(T_Component),intent(in) :: vCpn(:)
  real(dp),         intent(out):: tStoikio(:,:)
  !---------------------------------------------------------------------
  real(dp),allocatable:: A(:,:)
  integer, allocatable:: Indx(:)
  real(dp),allocatable:: Y(:)
  !
  ! stoikiometry table of the components in terms of the elements
  real(dp),allocatable:: tStoik(:,:)
  !
  logical :: bSingul
  real(dp):: D
  integer :: I,nEl,nCp,nFs
  !---------------------------------------------------------------------
  if(iDebug>0) write(fTrc,'(/,A)') "< Spl_Stoikio_Calc"
  !
  nFs=size(vFas)
  nCp=size(vCpn)
  nEl=size(vEle)
  !
  allocate(tStoik(nCp,nEl))
  do i=1,nCp
    tStoik(i,1:nEl)= vCpn(i)%vStoikCp(1:nEl) /real(vCpn(i)%vStoikCp(0))
  end do
  !
  allocate(A(1:nCp,1:nCp))
  A(1:nCp,1:nCp)= transpose(tStoik(1:nCp,2:nEl))
  ! column_1=Oxygen, -> already taken in account by Spl_ReadFormula
  ! -> columns 2:nEl of ElementFormula Matrix gives ElementStoikio of components 1:nCp
  !
  allocate(Indx(1:nCp))
  call LU_Decomp(A,Indx,D,bSingul)
  if(bSingul) call Stop_("SINGUL IN Spl_Stoikio_Calc")
  !-> inverse of A gives Elements 2:nEl in terms of components 1:nCp
  !
  ! obtain Phases 1:nFs expressed in terms of Components 1:nCp
  allocate(Y(1:nCp))
  do I=1,nFs
    Y= tFormula(I,2:nEl)
    call LU_BakSub(A,Indx,Y)
    tStoikio(I,1:nCp)=Y
  end do
  !
  deallocate(A)
  deallocate(Indx)
  deallocate(Y)
  !
  if(iDebug>2) call Spl_Stoikio_Show
  !
  deallocate(tStoik)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Spl_Stoikio_Calc"
  !
contains

subroutine Spl_Stoikio_Show
!--
!-- write stoikio on stokio "trace file"
!--
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut !,NamFInn
  !
  integer:: I,iEl
  integer:: FF
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< ShowStoikio"
  !
  if(iDebug>0) write(fTrc,'(A)') &
  & "Spl_Stoikio_Show -> "//trim(DirOut)//"_spl_stoikio.log"
  !
  call GetUnit(FF)
  !
  open(FF,file=trim(DirOut)//"_spl_stoikio.log")
  write(FF,'(A)') "!! system stoichiometry"
  !
  write(FF,'(/,A,/)') "Components in terms of elements"
  do iEl=1,nEl
    write(FF,'(A1,A2,A1)',advance="no") " ",vEle(iEl)%NamEl,T_
  end do
  write(FF,*)
  do I=1,nCp
    do iEl=1,nEl; write(FF,'(F7.2,A1)',advance="no") tStoik(I,iEl),T_; end do
    write(FF,'(A,A1,F15.3)') trim(vCpn(I)%NamCp),T_,vCpn(I)%Mole
  end do
  !
  write(FF,'(/,A,/)') "Phases in terms of elements"
  do I=1,nFs
    do iEl=1,nEl
       write(FF,'(F7.3,A1)',advance="no") tFormula(I,iEl),T_
    end do
    write(FF,'(A,A1,F15.3)') trim(vFas(I)%NamFs),T_,vFas(I)%Grt
  end do
  !
  close(FF)
  !
  !call GetUnit(FF)
  open(FF,file=trim(DirOut)//"_spl_stoikio2.log")
  write(FF,'(2A)') &
  & "!! Stoikio Table of Fases (line) in terms of Components (column)", &
  & " and vFas(I)%Grt"
  !
  do iEl=1,nCp
    write(FF,'(A,A1)',advance="no") trim(vCpn(iEl)%NamCp),T_
  end do
  write(FF,*)
  
  do I=1,nFs
    do iEl=1,nCp
      write(FF,'(F7.3,A1)',advance="no") tStoikio(I,iEl),T_
    end do
    write(FF,'(A,A1,F15.3)') trim(vFas(I)%NamFs),T_,vFas(I)%Grt
  end do
  !
  close(FF)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ ShowStoikio"
end subroutine Spl_Stoikio_Show

end subroutine Spl_Stoikio_Calc

subroutine Table_Transform(tInput,tTransform,tOutput)
!--
!-- transform, using tTransform, tInput to tOutput
!-- tTransform is invertible, (nC,nC)
!-- tInput and tOutput are matrices of same shape (nC,nS), nS>nC
!-- tOutput - inv(tTransform) * tInput
!--
  use M_Trace,only: Stop_
  use M_Numeric_Mat,only: LU_Decomp, LU_BakSub
  real(dp),intent(inout) :: tTransform(:,:) !square, non singular, (nC,nC)
  real(dp),intent(in)    :: tInput(:,:)  !(nC,nS)
  real(dp),intent(out)   :: tOutput(:,:) !(nC,nS)
  !
  real(dp),allocatable:: vY(:)
  integer, allocatable:: vIndx(:)
  logical :: bSingul
  real(dp):: D
  integer :: I,nC,nS
  !
  nC= size(tInput,1)
  nS= size(tInput,2)
  !
  allocate(vY(1:nC))
  allocate(vIndx(1:nC))
  !
  !the transformation matrix, tTransform, must be invertible !!
  call LU_Decomp(tTransform,vIndx,D,bSingul)
  if(bSingul) call Stop_("SINGUL IN Table_Transform")
  do I= 1,nS
    vY= tInput(:,I)
    call LU_BakSub(tTransform,vIndx,vY)
    tOutput(:,I)= vY(:)
  end do
  !
  deallocate(vY,vIndx)
  !
end subroutine Table_Transform

subroutine WriteRMat(f,Mat) !,N1,N2,M1,M2)
  integer, intent(in)::f !=File Index
  real(dp),intent(in)::Mat(:,:)
  !
  integer::I,J
  !
  do I=1,size(Mat,1)
    do J=1,size(Mat,2)
      write(f,'(F7.1, A1)',advance="no") Mat(I,J), T_
    end do
    write(f,*)
  end do
  write(f,*)
  !
end subroutine WriteRMat

subroutine WriteMat_I(f,Mat)
  integer,intent(in)::f !=File Index
  integer,intent(in)::Mat(:,:)
  !
  integer::I,J
  !
  do I=1,size(Mat,1)
    do J=1,size(Mat,2)
      write(f,'(I3, A1)',advance="no") Mat(I,J), T_
    end do
    write(f,*)
  end do
  write(f,*)
  !
end subroutine WriteMat_I

subroutine Spl_WriteSystem( &
& vFas,                     &
& tFormula,                 &
& vCpn,tStoikio)
!--
!-- write system stoikiometry
!--
  use M_IoTools,    only: GetUnit
  use M_Files,      only: DirOut !,NamFInn
  use M_T_Component,only: T_Component
  use M_T_Phase,    only: T_Phase
  !
  ! use M_Global_Vars,only: tFormula,vFas
  !
  type(T_Phase),    intent(in):: vFas(:)
  real(dp),         intent(in):: tFormula(:,:)
  type(T_Component),intent(in):: vCpn(:)
  real(dp),         intent(in):: tStoikio(:,:)

  integer :: nCp,nFs,iCpn,iFs
  integer :: FF
  real(dp):: X
  !
  nCp=size(vCpn)
  nFs=size(vFas)
  !
  write(fTrc,'(/,A)') "< WriteSystem"
  write(fTrc,'(A)') "Spl_WriteSystem -> "//trim(DirOut)//"_spl_system.log"

  call GetUnit(FF)
  open(FF,file=trim(DirOut)//"_spl_system.log")

  write(FF,'(A,A1)',advance="no") "_", T_
  do iCpn=1,nCp
    write(FF,'(A,A1)',advance="no") trim(vCpn(iCpn)%NamCp), T_
  end do
  write(FF,'(A)') "GIBBS"

  write(FF,'(A,A1)',advance="no") "MOLES", T_
  do iCpn=1,nCp
    write(FF,'(G12.3,A1)',advance="no") vCpn(iCpn)%Mole, T_
  end do
  write(FF,'(A)') "0.000"

  do iFs=1,nFs
    write(FF,'(A,A1)',advance="no") trim(vFas(iFs)%NamFs), T_
    do iCpn=1,nCp
      write(FF,'(F7.2,A1)',advance="no") tStoikio(iFs,iCpn), T_
    end do
    X= vFas(iFs)%Grt
    if(tFormula(iFs,1)/=0) X= vFas(iFs)%Grt /tFormula(iFs,1)
    ! write(FF,'(2F15.3)') vFas(iFs)%Grt,X
    write(FF,'(2F15.3)') vFas(iFs)%Grt
  end do

  write(fTrc,'(/,A,/)') "</ WriteSystem"

  close(FF)

  return
end subroutine Spl_WriteSystem

end module M_Simplex_Build

