module M_Basis_Tools
!--
!-- tools used by M_Basis for operations on basis (stoichiometry)
!--
  use M_Kinds
  use M_Trace,only: fTrc,fHtm,T_,Stop_,iDebug,Pause_
  !
  implicit none
  !
  private
  !
  public:: Basis_Dimensions
  public:: Basis_Alloc
  public:: Basis_Calc_iO_iH_iOx
  public:: Basis_MixPhase_CpnUpdate
  public:: Basis_SpcTypes_Build
  public:: Basis_Cpn_Order
  public:: Basis_SpcOrdPrm_Build
  public:: Basis_AlfaNu_Build
  public:: Basis_FindIBal
  public:: Basis_Cpn_ISpc_Select
  public:: Basis_Species_Select_Test
  public:: Basis_FreeSet_Select
  public:: Basis_StoikTable_Sho

  !~ public:: Basis_NuTable_Change
  !~ public:: Basis_Calc_iW_iH_iOH_iO2
  !~ public:: Basis_RedoxAssign

contains

subroutine Basis_Alloc(nSp,nCp,nFs,nAq,nAs,nMs)
  use M_Basis_Vars,only: tStoikio
  use M_Basis_Vars,only: Basis_SpcTypes_Alloc,Basis_SpcOrdPrm_Alloc,Basis_AlfaNu_Alloc
  !
  integer,intent(in):: nSp,nCp,nFs,nAq,nAs,nMs
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Alloc"
  !
  !---allocate logical tables vLxx>
  call Basis_SpcTypes_Alloc(nSp)
  !
  !---allocate vOrdXx,vPrmXx>
  call Basis_SpcOrdPrm_Alloc(nCp,nSp,nAq,nAs,nMs)
  !
  !---allocate stoichio tables tAlfXx,tNuXx>
  call Basis_AlfaNu_Alloc(nCp,nSp,nFs,nAs,nMs)
  !
  if(allocated(tStoikio)) deallocate(tStoikio)
  allocate(tStoikio(1:nCp,1:nSp))
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Alloc"
  !
end subroutine Basis_Alloc

subroutine Basis_Calc_iO_iH_iOx( & !
& vEle,vCpn, &    !in
& icO_,icH_,icOx) !out
!--
!-- find indexes of elements O, H, OX in vCpn (may change when basis change)
!--
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  !
  type(T_Element),  dimension(:),intent(in) :: vEle
  type(T_Component),dimension(:),intent(in) :: vCpn
  integer,                       intent(out):: icO_,icH_,icOx
  ! ranks of elements O__,H__,OX_ in Cpn of current component list
  integer:: I
  !
  icO_=0; icH_=0; icOx=0
  do I=1,size(vCpn)
    select case(vEle(vCpn(I)%iEle)%NamEl)
    case("O__"); icO_= I
    case("H__"); icH_= I
    case("OX_"); icOx= I
    end select
  end do
  if(iDebug==5) then
    print '(A,3I3)',"icO_icH_icOx= ",icO_,icH_,icOx
    call Pause_
  end if
end subroutine Basis_Calc_iO_iH_iOx

subroutine Basis_Dimensions( & !
& BufferIsExtern,        & !IN
& vCpn,vSpc,nAq,nMn,nGs, & !IN
& nCi,nMx,nAx,nMs,nAs)     !OUT
!--
!-- count number of Inert / Mobile components
!-- assign species' statut
!--
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  use M_T_Phase,    only: T_Phase
  !
  logical,          intent(in) :: BufferIsExtern
  type(T_Component),intent(in) :: vCpn(:)
  type(T_Species),  intent(in) :: vSpc(:)
  integer,          intent(in) :: nAq,nMn,nGs
  integer,          intent(out):: nCi,nMx,nAx,nMs,nAs
  !
  integer:: nCp,nSp,iSp,iCp,N
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Dimensions"
  !
  N=nGs !currently not used ...
  nCp= size(vCpn)
  nSp= size(vSpc)
  nMx= 0
  nAx= 0
  !
  do iCp=1,nCp
    !
    if(vCpn(iCp)%Statut=="INERT") cycle 
    !--species not always assigned at this stage--
    !
    iSp=vCpn(iCp)%iSpc !-> index of Prim'species in vSpc(1:nSp)
    !
    select case(vSpc(iSp)%Typ)

      case("AQU")
        !!!<BUFFER MODIFIED
        if(BufferIsExtern) then
          if(vCpn(iCp)%Statut=="BUFFER") &
          & call Stop_("Buffer species should be NON AQUEOUS !!!...")
        else
          if(vCpn(iCp)%Statut=="BUFFER") nAx=nAx+1
        end if
        !!!</BUFFER MODIFIED
        if(vCpn(iCp)%Statut=="MOBILE") nAx=nAx+1

      case("MIN"); nMx=nMx+1 !nr of Mobile Min.Species

      case("GAS"); nMx=nMx+1 !nr of Mobile Gas.Species (currently not different from Min ...)

    end select
    !
  end do
  !
  nCi= nCp -nAx -nMx
  nMs= nMn -nMx
  nAs= nAq -nCi -nAx
  !
  if(iDebug>2) write(fTrc,'(7(A,I4))') &
  & "nCp=",  nCp,"/ nCi=",nCi,"/ nAx=",nAx,"/ nAs=",nAs, &
  & "/ nMx=",nMx,"/ nMs=",nMs,"/ nMn=",nMn
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Dimensions"

  return
end subroutine Basis_Dimensions

subroutine Basis_MixPhase_CpnUpdate( &
& vMixFas,   & !in
& vSpc,vCpn)   !inout
!--
!-- update activities of end-members in non-aqueous solutions
!-- update activities of mobile components
!--
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component
  use M_T_MixPhase, only: T_MixPhase
  !
  type(T_MixPhase), intent(in)   :: vMixFas(:)
  type(T_Species),  intent(inout):: vSpc(:)
  type(T_Component),intent(inout):: vCpn(:)
  !
  integer :: iCp,iSp
  !
  if(iDebug>0) write(fTrc,'(/,A)') "<== Basis_MixPhase_CpnUpdate"
  !
  !------------------update activities of externally buffered components
  do iCp=1,size(vCpn)
    iSp= vCpn(iCp)%iSpc
    !
    if(vCpn(iCp)%iMix/=0) then
    !--- if prim'species is in non'aqu'sol', then component is mobile --
      vSpc(iSp)%Dat%LAct= vMixFas(vCpn(iCp)%iMix)%vLnAct(vCpn(iCp)%iPol)
      vCpn(iCp)%LnAct= vSpc(iSp)%Dat%LAct
    end if
    !
    if(iDebug>0) write(fTrc,'(A,I3,A1,A,A1,G15.6)') &
    & "CPN=",iCp,T_,vSpc(iSp)%NamSp,T_,vSpc(iSp)%Dat%LAct
  end do
  !-----------------/
  !
  if(iDebug>0) write(fTrc,'(A,/)') "<==/ Basis_MixPhase_CpnUpdate"
  !
  return
end subroutine Basis_MixPhase_CpnUpdate

subroutine Basis_SpcTypes_Build( & !
& vSpc,vCpn, & ! IN
& vLCi,vLAx,vLMx,vLAs,vLMs) ! OUT
!--
!-- build the logical arrays vLCi,vLAx,vLMx,vLAs,vLMs,...
!--
  use M_T_Element,  only: T_Element
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component
  !
  type(T_Species),     intent(in) :: vSpc(:)
  type(T_Component),   intent(in) :: vCpn(:)
  logical,dimension(:),intent(out):: vLCi,vLAx,vLMx,vLAs,vLMs
  !
  integer:: iCp,iSp
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_SpcTypes_Build"
  !
  !-------------------------------------------update logical tables vLxx
  vLCi=.false.
  vLAx=.false.
  vLMx=.false.
  vLAs(:)= vSpc(:)%Typ(1:3)=="AQU"
  vLMs(:)= vSpc(:)%Typ(1:3)=="MIN" .or. vSpc(:)%Typ(1:3)=="GAS"
  !vLGs(:)= vSpc(:)%Typ=="GAS"
  !
  do iCp=1,size(vCpn)
    !
    iSp= vCpn(iCp)%iSpc !-> iSp is primary species
    !
    vLAs(iSp)= .false.
    vLMs(iSp)= .false.
    !
    select case(trim(vCpn(iCp)%Statut))
    !
    case("INERT","BALANCE") !Intern'Species,
      vLCi(iSp)=.true.
      !includes "balance" species for charge adjustment (in case pH fixed)
    !
    case("MOBILE","BUFFER")
      select case(vSpc(iSp)%Typ(1:3))
        case("AQU")       ; vLAx(iSp)=.true. !eXtern'Aqu.Species
        case("MIN","GAS") ; vLMx(iSp)=.true. !eXtern'Min.Species
      end select
    !
    end select
  end do
  !-----------------------------------------------/update logical tables
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_SpcTypes_Build"
  return
end subroutine Basis_SpcTypes_Build

subroutine Basis_Cpn_Order( & !
& vSpc,nCi,nAx, & !in
& vCpn)           !out
!--
!-- sort the components, inert first, mobile next --
!--
  use M_Numeric_Const,only: Ln10
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component
  !~ use M_Basis_Vars, only: iEleCpn,iCpnEle
  !~ use M_Basis_Vars, only: Basis_CpnPrm_Alloc
  !
  type(T_Species),  intent(inout):: vSpc(:) !vSpc(iSp)%Dat%LAct=vCpn(iCp)%LnAct
  integer,          intent(in)   :: nCi,nAx
  type(T_Component),intent(inout):: vCpn(:)
  !
  type(T_Component),dimension(:),allocatable::vC !-> re-ordered vCpn
  integer:: iCp,iCi,iAx,iMx,iSp !,nSp,nCp
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Cpn_Order"
  !
  if(iDebug>0) then
    write(fTrc,'(A)') "Components->Species"
    do iCp=1,size(vCpn)
      write(fTrc,'(A,I3,A1,A)') &
      & "Cpn",iCp,T_,trim(vSpc(vCpn(iCp)%iSpc)%NamSp)
    end do
  end if
  !
  iCi=0 !index for Internal Prim'Aqu'Species
  iAx=0 !index for mobile Prim'Aqu'species, begins after nCi in vCpn
  iMx=0 !index of mobile Min.species, begins after Ci+nAx in vEle
  !
  !--------------------------reorder components, Inert / MobAqu / MobMin
  allocate(vC(1:size(vCpn)))
  !
  do iCp=1,size(vCpn)
    !
    iSp=vCpn(iCp)%iSpc !-> iSp is primary species
    !
    select case(trim(vCpn(iCp)%Statut))

    case("INERT","BALANCE") !Intern'Species
      !--------------------------------------- inert components first --
      iCi=iCi+1
      vC(iCi)=vCpn(iCp)
      !
    case("MOBILE","BUFFER")
      vSpc(iSp)%Dat%LAct= vCpn(iCp)%LnAct
      select case(vSpc(iSp)%Typ(1:3))

      case("AQU") !eXtern'Aqu.Species
        !---------------------------------------- then mobile aqueous --
        iAx=iAx+1
        vC(nCi+iAx)=vCpn(iCp)

      case("MIN","GAS") !eXtern'Min.Species
        !----------------------------------------- mobile non'aqueous --
        iMx=iMx+1
        vC(nCi+nAx+iMx)=vCpn(iCp)

      end select

    end select
    !
  end do
  !
  vCpn(:)= vC(:)
  !
  deallocate(vC)
  !-------------------------/reorder components, Inert / MobAqu / MobMin
  !
  !-------------------------------------------- build iEleCpn,iCpnEle --
  !~ call Basis_CpnPrm_Alloc(size(vCpn))
  !~ call Basis_CpnPrm_Calc(vCpn,iEleCpn,iCpnEle)
  !-------------------------------------------/ build iEleCpn,iCpnEle --
  !
  !----------------------------------------------------------------trace
  if(iDebug>0) then
    write(fTrc,'(/,A,/)') "Components, Re-Ordered"
    write(fTrc,'(2(A,/))') &
    & "iCp,Cpn%Statut,vSpc(Cpn%iSpc)%NamSp,vCpn(iCp)%Mole,Cpn%LnAct/Ln10"
    !
    do iCp=1,size(vCpn)
      write(fTrc,'(I3,A1,2(A,A1),2(G15.6,A1))') &
      & iCp,                        T_,&
      & vCpn(iCp)%Statut,           T_,&
      & vSpc(vCpn(iCp)%iSpc)%NamSp, T_,&
      & vCpn(iCp)%Mole,             T_,&
      & vCpn(iCp)%LnAct/Ln10,       T_
    end do
    !call Pause_
  end if
  !---------------------------------------------------------------/trace
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Cpn_Order"
  !
end subroutine Basis_Cpn_Order

! subroutine Basis_CpnPrm_Calc( & !
! & vCpn,             & ! IN
! & iEleCpn,iCpnEle)    ! OUT
!   use M_T_Element,  only: T_Element
!   use M_T_Component,only: T_Component
!   !
!   type(T_Component),intent(in) :: vCpn(:)
!   integer,          intent(out):: iEleCpn(:),iCpnEle(:)
!   !
!   integer:: I,J
!   !
!   iCpnEle(:)=vCpn(:)%iEle
!   !-> iEl=iCpnEle(iCp)= which Ele is associated with vCpn(iCp)
!   !
!   do I=1,size(vCpn)
!     do J=1,size(vCpn)
!       if(J==iCpnEle(I)) iEleCpn(I)=J
!     end do
!     !-> vCpn(iEleCpn(iEl)) is associated with element iEl
!   end do
!   !
! end subroutine Basis_CpnPrm_Calc

subroutine Basis_SpcOrdPrm_Build( &
& vSpc,vCpn, &
& vLCi,vLAx,vLMx,vLAs,vLMs, &
& vOrdPr,vOrdAq,vOrdAs,vOrdMs, &
& vPrmFw,vPrmBk)
!--
!-- build vOrdPr, vOrdAq, vOrdAs, vOrdMs, vPrmFw, vPrmBk,
!-- variables used- vLAs,vLMs,.., nCi,nAx,nMx,nAs,nMs, ...
!-- vPrmFw(J) is index, in vSpc, of the J-th species in the current array
!-- vPrmBk(J) is index in the current array of the J-th species of vSpc
!--
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component
  !
  type(T_Species),  dimension(:),intent(in) :: vSpc
  type(T_Component),dimension(:),intent(in) :: vCpn
  logical,          dimension(:),intent(in) :: vLCi,vLAx,vLMx
  logical,          dimension(:),intent(in) :: vLAs,vLMs
  integer,          dimension(:),intent(out):: vOrdPr,vOrdAq,vOrdAs,vOrdMs
  integer,          dimension(:),intent(out):: vPrmFw,vPrmBk
  !
  integer:: I,J,nCp,nSp
  integer:: nAq,nCi,nAs,nAx,nMx,nMs
  ! integer:: iCi,iAx,iMx
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_SpcOrdPrm_Build"
  !
  nSp= size(vSpc)
  nCp= size(vCpn)
  nAq= count(vSpc%Typ(1:3)=="AQU")
  !
  nCi= count(vLCi)
  nAs= count(vLAs)
  nAx= count(vLAx)
  nMx= count(vLMx)
  nMs= count(vLMs)
  !
  !--- Aqueous Inert Prim'Species --
  vOrdPr(1:nCi)= vCpn(1:nCi)%iSpc
  !--- Aqueous Mobile Species --
  if(nAx>0) vOrdPr(nCi+1:nCi+nAx)= vCpn(nCi+1:nCi+nAx)%iSpc
  !--- non Aqu'Mobile Species --
  if(nMx>0) vOrdPr(nCi+nAx+1:nCi+nAx+nMx)= vCpn(nCi+nAx+1:nCi+nAx+nMx)%iSpc
  !
  !--- Aqueous Inert Prim'Species --
  vOrdAq(1:nCi)= vCpn(1:nCi)%iSpc
  !--- Aqueous Mobile Species --
  if(nAx>0) &
  & vOrdAq(nCi+nAs+1:nCi+nAs+nAx)= vCpn(nCi+1:nCi+nAx)%iSpc
  !
  !------------------ scan vSpc for aqu. second. species -> build vOrdAs 
  if(nAs>0) then !
    J=0
    do I=1,nSp
      if(vLAs(I)) then !vSpc(I)%Typ=="AQU"
        J=J+1
        vOrdAs(J)=I
        !-> vOrdAs(J, J=1:nAs) points to index of Jth sec'aqu'species in vSpc
        vOrdAq(nCi+J)=I
      end if
    end do
  end if
  !------------------/
  !
  !-------------------scan vSpc for min. second. species -> build vOrdMs
  if(nMs>0) then
    J=0
    do I=1,nSp
      if(vLMs(I)) then !vSpc(I)%Typ=="MIN" .or. vSpc(:)%Typ=="GAS"
        J=J+1
        vOrdMs(J)=I
        !-> vOrdMs(J, J=1:nMs) points to index of Jth sec'min'species in vSpc
      end if
    end do
  end if
  !-------------------/
  !
  !nCi          //nAs           //nAx           //nMx       //nMs
  !Aqu.Inert.Spc//Aqu.Second.Spc//Aqu.Mobile.Spc//Min.Mobile//Min.notMobile
  !1            //nCi+1         //nCi+nAs+1     //nAq+1     //nAq+nMx+1
  vPrmFw(1:nCi)=                            vOrdPr(1:nCi)
  if(nAx>0) vPrmFw(nCi+nAs+1:nCi+nAs+nAx)=  vOrdPr(nCi+1:nCi+nAx)
  if(nMx>0) vPrmFw(nAq+1:nAq+nMx)=          vOrdPr(nCi+nAx+1:nCi+nAx+nMx)
  if(nAs>0) vPrmFw(nCi+1:nCi+nAs)=          vOrdAs(1:nAs)
  if(nMs>0) vPrmFw(nAq+nMx+1:nAq+nMx+nMs)=  vOrdMs(1:nMs)
  !
  do I=1,nSp
    vPrmBk(vPrmFw(I))=I
  end do
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_SpcOrdPrm_Build"
  !
end subroutine Basis_SpcOrdPrm_Build

subroutine AlfaTable_Build( &
& nCp,nSp,tStoikio, &
& vIndexBuffer, &
& vTot, &
& tAlfSp)
!--
!-- build tAlfSp from tStoikio(1:nCp,1:nSp)
!--
  use M_Basis_Ortho
  !
  integer,  intent(in) :: nCp,nSp
  real(dp), intent(in) :: tStoikio(:,:)
  integer,  intent(in) :: vIndexBuffer(:)
  !logical,  intent(in) :: HasBuffer
  real(dp), intent(inout):: vTot(:)
  real(dp), intent(out)  :: tAlfSp(:,:)
  !
  real(dp),allocatable:: tXOrthoBasis(:,:)
  != indexes of buffer components in tStoikio
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< AlfaTable_Build"
  !
  !tStoikio:column iSp = stoikio of species iSp in terms of elements (1:nCp)
  !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !
  tAlfSp(1:nCp,1:nSp)= tStoikio(1:nCp,1:nSp)
  !
  if(ANY(vIndexBuffer(:)>0)) then
    !
    allocate(tXOrthoBasis(nCp,nCp))
    !
    if(iDebug>2) then
      call Basis_Ortho_Compute_Debug(vIndexBuffer,tAlfSp,tXOrthoBasis)
    else
      call Basis_Ortho_Compute(vIndexBuffer,tAlfSp,tXOrthoBasis)
    end if
    !
    tAlfSp= matmul(tXOrthoBasis,tAlfSp)
    vTot=   matmul(tXOrthoBasis,vTot)
    !
    deallocate(tXOrthoBasis)
    !
  end if
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ AlfaTable_Build"
  !
end subroutine AlfaTable_Build

subroutine AlfaTable_Extract( &
& vCpn,tAlfSp,vLCi,vLAx,vLMx,vLAs,vLMs, &
& tAlfPr,tAlfAs,tAlfMs)
!-----------------------------------------------------------------------
!.extract blocks from tAlfSp(1:nCp,1:nSp) to form tAlfPr,tAlfAs,tAlfMs
!.according to vLCi,vLAx,vLMx,vLAs,vLMs
!-----------------------------------------------------------------------
  use M_T_Component,only: T_Component
  !
  type(T_Component),intent(in):: vCpn(:)
  real(dp), intent(in):: tAlfSp(:,:)
  logical,  intent(in):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !
  real(dp), intent(out):: tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !
  integer :: I,J
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< AlfaTable_Extract"
  !
  !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !
  !--extract blocks from tStoikio(1:nCp,1:nSp) to tAlfPr,tAlfAs,tAlfMs--
  J= 0
  !-------------------------------------subblock of inert prim'species--
  do I=1,count(vLCi)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  end do
  !-------------------------------------subblock of mobile aqu'species--
  do I=1,count(vLAx)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  end do
  !-------------------------------------subblock of mobile min'species--
  do I=1,count(vLMx)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  end do
  !
  tAlfAs(:,:)= Zero
  !---------------------------------extract columns of sec'aqu.species--
  if(count(vLAs)>0) call Extract(2,tAlfSp,vLAs,tAlfAs)
  !
  tAlfMs(:,:)= Zero
  !-----------------------------extract columns of sec'min.gas.species--
  if(count(vLMs)>0) call Extract(2,tAlfSp,vLMs,tAlfMs)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ AlfaTable_Extract"
  !
end subroutine AlfaTable_Extract

subroutine AlfaTable_CalcBuffered( & !
& vCpn,vSpc,tStoikio,vLCi,vLAs,vLMs, & !
& vTot, &
& tAlfPr,tAlfAs,tAlfMs) !!! NOT useD ???
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  use M_Numeric_Mat,only: LU_BakSub
  !
  type(T_Component),intent(in):: vCpn(:)
  type(T_Species),  intent(in):: vSpc(:)
  real(dp),         intent(in):: tStoikio(:,:)
  logical,          intent(in):: vLCi(:),vLAs(:),vLMs(:)
  real(dp),         intent(inout):: vTot(:)
  real(dp),         intent(inout):: tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !
  real(dp),allocatable:: tTransform(:,:)
  integer, allocatable:: vIndx(:)
  real(dp),allocatable:: vY(:)
  logical :: bSingul
  integer :: I,nCp,nSp,nCi
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< AlfaTable_CalcBuffered"
  !
  nCp= size(vCpn)
  nSp= size(vSpc)
  nCi= count(vLCi)
  !
  !for mole balance operations, replace the elements by "standard" primary species
  != "INERT" component: single ion with default valency, e.g. Na->Na+,... -> no modif'
  != "mobile" component: use species as component
  !  -> elements replaced by corresponding species, e.g. C by CO2
  !
  allocate(tTransform(nCp,nCp))
  allocate(vIndx(nCp))
  allocate(vY(nCp))
  !
  tTransform(1:nCp,1:nCp)=Zero
  do I=1,nCp; tTransform(I,I)=One; end do !identity matrix
  !
  if(nCp>nCi) then
  ! for buffered species, replace element by species
  ! using corresponding column of formula matrix
    do I=nCi+1,nCp
      tTransform(:,I)= tStoikio(:,vCpn(I)%iSpc)
    end do
  end if
  !
  call Compute_Transform(tTransform,vIndx,bSingul)
  if(bSingul) call Stop_("SINGUL IN CalcAlfaTableBuffered")
  !-> tTransform gives Elements as functions of Int'Elements/Ext'Species
  !
  !-----------------------compute new total amounts of component species
  call LU_BakSub(tTransform,vIndx,vTot(1:nCp))
  !----------------------/
  !do I=1,nSp
  !  vY=tAlfSp(1:nCp,I)
  !  call LU_BakSub(tTransform,vIndx,vY)
  !  tAlfSp(1:nCp,I)=vY(1:nCp)
  !end do
  !
  do I=1,nCp
    vY=tAlfPr(1:nCp,I)
    call LU_BakSub(tTransform,vIndx,vY)
    tAlfPr(1:nCp,I)=vY(1:nCp)
  end do
  !
  do I=1,count(vLAs)
    vY=tAlfAs(1:nCp,I)
    call LU_BakSub(tTransform,vIndx,vY)
    tAlfAs(1:nCp,I)=vY(1:nCp)
  end do
  !
  do I=1,count(vLMs)
    vY=tAlfMs(1:nCp,I)
    call LU_BakSub(tTransform,vIndx,vY)
    tAlfMs(1:nCp,I)=vY(1:nCp)
  end do
  !
  deallocate(tTransform)
  deallocate(vIndx)
  deallocate(vY)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ AlfaTable_CalcBuffered"
end subroutine AlfaTable_CalcBuffered

subroutine NuTable_Build( & !
& vIndxPrimSpc,tStoikio,  & !IN
& tNuSp)                    !OUT
  use M_Numeric_Mat,only: LU_BakSub
  !
  integer, intent(in) :: vIndxPrimSpc(:)
  real(dp),intent(in) :: tStoikio(:,:)
  real(dp),intent(out):: tNuSp(:,:)
  !
  real(dp),allocatable:: tTransform(:,:)
  integer, allocatable:: vIndx(:)
  real(dp),allocatable:: vY(:)
  logical :: bSingul
  integer :: I,J,nCp
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< NuTable_Build"
  !
  nCp= size(tStoikio,1)
  !
  allocate(tTransform(nCp,nCp))
  !
  !-------------------assemble tTransform from the columns of tStoikio--
  !-----------------------------------------that point to Prim'Species--
  J=0
  do I=1,nCp
    J=J+1
    tTransform(:,J)= tStoikio(:,vIndxPrimSpc(I))
  end do
  !
  !--------------------------> tTransform is tStoikio(nCi | nAx | nMx)--
  !! iCi= 0
  !! iAx= count(vLCi)
  !! iMx= count(vLCi) +count(vLAx)
  !! do I=1,size(tNuSp,1)
  !!   if(vLCi(I)) then
  !!     iCi= iCi+1
  !!     tTransform(:,iCi)= tStoikio(:,vCpn(iCi)%iSpc)
  !!   elseif(vLAx(I)) then
  !!     iAx= iAx+1
  !!     tTransform(:,iAx)= tStoikio(:,vCpn(iAx)%iSpc)
  !!   elseif(vLMx(I)) then
  !!     iMx= iMx+1
  !!     tTransform(:,iMx)= tStoikio(:,vCpn(iMx)%iSpc)
  !!   end if
  !! end do
  !
  allocate(vIndx(nCp))
  call Compute_Transform(tTransform,vIndx,bSingul)
  !
  if(bSingul) call Stop_("SINGUL IN NuTable_Build")
  !
  allocate(vY(nCp))
  !
  !-- in tNuSp, species are in same order as in vSpc --
  do I=1,size(tNuSp,1)
    vY= tStoikio(1:nCp,I)
    call LU_BakSub(tTransform,vIndx,vY)
    !Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    tNuSp(I,1:nCp)= vY(1:nCp)
    !-> this includes transposition between tAlfa and tNu !!
  end do
  !
  deallocate(tTransform)
  deallocate(vIndx)
  deallocate(vY)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ NuTable_Build"
  !
end subroutine NuTable_Build
  !
subroutine NuTable_Extract(  & !
& tNuSp,vLAs,vLMs, & !
& tNuAs,tNuMs)
  !
  real(dp),intent(in) :: tNuSp(:,:)
  logical, intent(in) :: vLAs(:),vLMs(:)
  real(dp),intent(out):: tNuAs(:,:),tNuMs(:,:)
  !
  !----extract blocks from tNuSp, for the different types of sec'species
  !----------------------------------------------- second'aqu'species --
  if(count(vLAs)>0) call Extract( & !
  & 1,tNuSp,vLAs, & !IN
  & tNuAs)          !OUT
  !----------------------------------------------- second'min'species --
  if(count(vLMs)>0) call Extract( & !
  & 1,tNuSp,vLMs, & !IN
  & tNuMs)          !OUT
  !------------------------------------------------------/extract blocks
  !
end subroutine NuTable_Extract

! subroutine Basis_NuTable_Change(  & !
! & vIndxPrimSpc, & !IN
! & tStoikio,     & !IN
! & tNuSp)          !OUT
! !--
! !-- NuTable_Build, simplified version for general use
! !-- - no extraction of tNuAs, tNuMs
! !--
!   use M_Numeric_Mat,only: LU_BakSub
!   !
!   integer, intent(in) :: vIndxPrimSpc(:)
!   real(dp),intent(in) :: tStoikio(:,:)
!   real(dp),intent(out):: tNuSp(:,:)
!   !
!   real(dp),allocatable:: tTransform(:,:)
!   integer, allocatable:: vIndx(:)
!   real(dp),allocatable:: vY(:)
!   logical :: bSingul
!   integer :: I,J,nCp
!   !
!   if(iDebug>0) write(fTrc,'(/,A)') "< Basis_NuTable_Change"
!   !
!   nCp= size(tStoikio,1)
!   !
!   allocate(tTransform(nCp,nCp))
!   !--------------------- assemble tTransform from columns of tStoikio --
!   J=0
!   do I=1,nCp
!     J=J+1
!     tTransform(:,J)= tStoikio(:,vIndxPrimSpc(I))
!   end do
!   !-------------------------> tTransform is tStoikio(nCi | nAx | nMx) --
!   !
!   allocate(vIndx(nCp))
!   call Compute_Transform(tTransform,vIndx,bSingul)
!   !
!   if(bSingul) call Stop_("SINGUL IN Basis_NuTable_Change")
!   !
!   allocate(vY(nCp))
!   !
!   !-- in tNuSp, species will be in same order as in tStoikio, same as in vSpc --
!   do I=1,size(tNuSp,1)
!     vY= tStoikio(1:nCp,I)
!     call LU_BakSub(tTransform,vIndx,vY)
!     !Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
!     tNuSp(I,1:nCp)= vY(1:nCp)
!     !-> this includes transposition between tAlfa and tNu !!
!   end do
!   !
!   deallocate(tTransform)
!   deallocate(vIndx)
!   deallocate(vY)
!   !
!   if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_NuTable_Change"
!   !
! end subroutine Basis_NuTable_Change

subroutine Compute_Transform(tTransform,vIndx,Error)
  use M_Numeric_Mat,only: LU_Decomp
  real(dp),intent(inout):: tTransform(:,:)
  integer, intent(out)  :: vIndx(:)
  logical, intent(out)  :: Error
  !
  real(dp):: D
  !
  !--the transformation matrix, tTransform, must be invertible !!
  call LU_Decomp(tTransform,vIndx,D,Error)
  !
end subroutine Compute_Transform

subroutine Fas_NuTable_Build( & !
& vFas,vMixFas,vMixModel,tNuSp, & !in
& tNuFas)   !out
  use M_T_Phase,    only: T_Phase
  use M_T_MixPhase, only: T_MixPhase
  use M_T_MixModel, only: T_MixModel
  !
  type(T_Phase),     intent(in) :: vFas(:)
  type(T_MixPhase),  intent(in) :: vMixFas(:)
  type(T_MixModel),  intent(in) :: vMixModel(:)
  real(dp),          intent(in) :: tNuSp(:,:)
  !
  real(dp),          intent(out):: tNuFas(:,:)
  !
  type(T_MixModel):: SM
  type(T_MixPhase):: SF
  integer:: I,J,K,iCp,NPole !,nFs
  integer:: nCp
  !
  nCp=size(tNuSp,2)
  !
  do I=1,size(vFas)
    !
    if(vFas(I)%iSpc /= 0) then  ! "PURE","DISCRET"
      !
      K=vFas(I)%iSpc
      tNuFas(I,:)= tNuSp(K,:)
      !
    else if(vFas(I)%iMix /= 0) then ! "MIXT"
      !
      SF=vMixFas(vFas(I)%iMix) !-> the solution phase
      SM=vMixModel(SF%iModel)  !-> the corresponding model
      NPole=SM%NPole
      !
      ! where(.not.SM%vHasPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
      ! where(.not.SF%vLPole(1:NPole))   SF%vXPole(1:SM%NPole)=Zero
      !
      do iCp=1,nCp
        tNuFas(I,iCp)= Zero
        !
        do J=1,NPole
          if(SF%vLPole(J)) &
          & tNuFas(I,iCp)= tNuFas(I,iCp) &
          &              + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
        end do
        !
      end do
      !
    else if (vFas(I)%iSol /= 0) then
      !
      K=vFas(I)%iSpc
      tNuFas(I,:)= tNuSp(K,:)
    end if
  
    !! select case(vFas(I)%Typ)
    !! !
    !! case("PURE") !,"DISCRET")
    !!   K=vFas(I)%iSpc
    !!   tNuFas(I,:)= tNuSp(K,:)
    !! !--
    !! case("MIXT")
    !!   !
    !!   SF=vMixFas(vFas(I)%iMix) !-> the solution phase
    !!   SM=vMixModel(SF%iModel)  !-> the corresponding model
    !!   NPole=SM%NPole
    !!   !
    !!   where(.not.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
    !!   !
    !!   do iCp=1,nCp
    !!     tNuFas(I,iCp)= Zero
    !!     !
    !!     do J=1,NPole
    !!       if(SF%vLPole(J)) &
    !!       & tNuFas(I,iCp)= tNuFas(I,iCp) &
    !!       &              + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
    !!     end do
    !!     !
    !!   end do
    !! !--
    !! case("SOLU")
    !!   K=vFas(I)%iSpc
    !!   tNuFas(I,:)= tNuSp(K,:)
    !!   
    !! !--
    !! end select
    !
  end do
  !
end subroutine Fas_NuTable_Build

subroutine Fas_AlfaTable_Build( & !
& vFas,vMixFas,vMixModel, tAlfSp, & !in
& tAlfFs)   !out
  use M_T_Phase,    only: T_Phase
  use M_T_MixPhase, only: T_MixPhase
  use M_T_MixModel, only: T_MixModel
  !
  type(T_Phase),     intent(in) :: vFas(:)
  type(T_MixPhase),  intent(in) :: vMixFas(:)
  type(T_MixModel),  intent(in) :: vMixModel(:)
  real(dp),          intent(in) :: tAlfSp(:,:)
  !
  real(dp),          intent(out):: tAlfFs(:,:)
  !
  type(T_MixModel):: SM
  type(T_MixPhase):: SF
  integer:: I,J
  integer:: iPur,iMix,iSol,iCp
  integer:: nCp,NPole
  !
  nCp=size(tAlfSp,1)
  !
  do I=1,size(vFas)
    !
    iPur= vFas(I)%iSpc
    iMix= vFas(I)%iMix
    !
    if(iPur>0) then
      tAlfFs(:,I)= tAlfSp(:,iPur)
    end if
    !
    if(iMix>0) then
      !
      SF=vMixFas(iMix) !-> the solution phase
      SM=vMixModel(SF%iModel)  !-> the corresponding model
      NPole=SM%NPole
      !
      where(.not.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
      !
      do iCp=1,nCp
        tAlfFs(iCp,I)= Zero
        !
        do J=1,NPole
          if(SF%vLPole(J)) &
          & tAlfFs(iCp,I)= tAlfFs(iCp,I) &
          &              + tAlfSp(iCp,SM%vIPole(J))*SF%vXPole(J)
        end do
     end do
      !
    end if
  end do
  !
end subroutine Fas_AlfaTable_Build

subroutine Extract(cod,tIn,vL,tOut)
!--
!-- from matrix tIn,
!-- build the matrix tOut
!-- by extracting the rows (cod--1) or the columns (cod--2)
!-- that have the logical vL true
!--
  integer,                intent(in) :: cod
  real(dp),dimension(:,:),intent(in) :: tIn
  logical, dimension(:),  intent(in) :: vL
  real(dp),dimension(:,:),intent(out):: tOut
  integer::I,J
  J=0
  select case(cod)
    case(1)
      do I=1,size(tIn,1)
        if(vL(I)) then  ; J=J+1; tOut(J,:)=tIn(I,:)  ; end if
      end do
    case(2)
      do I=1,size(tIn,2)
        if(vL(I)) then  ; J=J+1; tOut(:,J)=tIn(:,I)  ; end if
      end do
  end select
end subroutine Extract

subroutine Basis_Cpn_ISpc_Select( &
& vSpc,isW,vCpn, & !in
& vIndex) !out
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  use M_T_Species,  only: T_Species
  use M_T_Component,only: T_Component
  use M_Numeric_Tools !-> sort, iMinLoc_R, iFirstLoc
  !
  type(t_Species),  intent(in) :: vSpc(:)
  integer,          intent(in) :: isW
  type(t_Component),intent(in) :: vCpn(:)
  integer,          intent(out):: vIndex(:)
  !
  integer :: nCp,nAq
  integer :: iAq,iPr,I,J,K,JJ
  real(dp):: Y,Pivot
  real(dp),allocatable:: tStoik(:,:)
  integer, allocatable:: vIPivot(:),vSpPrm(:)
  logical, allocatable:: vSpDone(:)
  logical, allocatable:: vCpDone(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Cpn_ISpc_Select"
  !
  nAq= count(vSpc%Typ=="AQU")
  nCp= size(vCpn)
  !
  vIndex(:)= vCpn(:)%iSpc
  !
  allocate(vSpPrm(1:nAq))
  allocate(vSpDone(1:nAq)); vSpDone(1:nAq)=.false.
  allocate(vCpDone(1:nCp)); vCpDone(1:nCp)=.false.
  !
  call CalcPermut(vSpc(1:nAq)%Dat%Mole,vSpPrm)
  ! -> sort aqu.species in order of increasing abudance
  ! -> permutation array vPrm
  ! to apply permutation vSpPrm to array Arr: Arr=Arr(vSpPrm)
  !
  allocate(tStoik(nCp,nCp))   ;  tStoik(:,:)= Zero
  allocate(vIPivot(nCp))      ;  vIPivot(:)= 0
  !
  K=   0
  iPr= 0
  !
  !--------------------- species 1 is solvent, component 1 is solvent --
  K=   K+1
  iPr= iPr+1
  !vCpn(iPr)%iSpc= isW
  vIndex(iPr)= isW
  vSpDone(isW)= .true.
  vCpDone(iPr)= .true.
  tStoik(:,K)= vSpc(isW)%vStoikio(1:nCp) !tFormula(:,isW)
  vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
  !! vIPivot(K)= 1 !!???!!!
  Pivot= tStoik(vIPivot(K),K)
  tStoik(:,K)= tStoik(:,K) /Pivot
  !
  !--------------------------------scan the BUFFER and MOBILE components
  do i=1,nCp
    !
    !-- read BUFFER and MOBILE stoikiometries in order to select,
    !-- in do00 block,
    !-- primary aqueous species independent from these species
    !
    if(vCpn(i)%Statut=="BUFFER" .or. vCpn(i)%Statut=="MOBILE") then
      !
      K= K+1
      vIndex(i)= vCpn(i)%iSpc
      tStoik(:,K)= vSpc(vCpn(i)%iSpc)%vStoikio(1:nCp)
      do J=1,K-1
        Y= tStoik(vIPivot(J),K)
        tStoik(:,K)= tStoik(:,K)- Y*tStoik(:,J)
        !-> substract column vIPivot(J) from column K
        !-> tStoik(vIPivot(J),K) will be 0
      end do
      vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
      !-> index of first non zero element
      if(iDebug>2) then
        do JJ=1,nCp
          write(fTrc,'(F7.2,1X)',advance="NO") tStoik(JJ,K)
        end do
        write(fTrc,'(2(A,I3),1X,A)') &
        & "K=",K, " /  iPr=",i,trim(vSpc(vIndex(i))%NamSp)
        !pause_
      end if
      Pivot= tStoik(vIPivot(K),K)
      tStoik(:,K)= tStoik(:,K) /Pivot
      !-> first non-zero element, tStoik(vIPivot(K),K), will be 1.
      !
      vCpDone(i)= .true.
      if(vSpc(vCpn(i)%iSpc)%Typ=="AQU") vSpDone(vCpn(i)%iSpc)= .true.
      !
    end if
    !
  end do
  !-------------------------------/scan the BUFFER and MOBILE components
  !
  Do00: do iPr=1,nCp
    !
    if(vCpDone(iPr)) cycle
    !
    K= K+1
    if(K>nCp) exit do00
    !
    I= nAq+1
    !-----------------------------------------------scan aqueous species
    DoAqu: do
      I= I-1
      if(I==0) exit DoAqu
      !
      !--scan species beginning with most abundant----------------------
      iAq= vSpPrm(i)
      !print *,"Spc=",trim(vSpc(iAq)%NamSp),vSpc(iAq)%Dat%Mole  ;  pause_
      if(vSpDone(iAq)) cycle DoAqu
      !
      tStoik(:,K)=  vSpc(iAq)%vStoikio(1:nCp) !tFormula(:,iAq)
      !
      if(tStoik(vCpn(iPr)%iEle,K)==0) cycle DoAqu
      !
      do J=1,K-1
        Y= tStoik(vIPivot(J),K)
        tStoik(:,K)= tStoik(:,K)- Y*tStoik(:,J)
        !-> tStoik(vIPivot(J),K) will be 0
      end do
      ! vIPivot(K)= ifirstLoc(abs(tStoik(:,K)) >1.D-9)
      vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
      !-> index of first non zero element
      !
      if(vIPivot(K)>nCp) then
      !--all coeff's are 0 -> this species is not independent-----------
      !---> cycle doAqu --
        !
        if(iDebug==4) then
          do JJ=1,nCp; write(fTrc,'(F7.2,1X)',advance="NO") tStoik(JJ,K); end do
          write(fTrc,'(A,I3,1X,2A)') &
          & "iFirst=",vIPivot(K),trim(vSpc(iAq)%NamSp),"-> removed"
        end if
        !
        vSpDone(iAq)= .true.
        cycle doAqu
      end if
      !
      if(iDebug==4) then
        do JJ=1,nCp; write(fTrc,'(F7.2,1X)',advance="NO") tStoik(JJ,K); end do
        write(fTrc,'(2(A,I3),1X,A)') &
        & "K=",K," /  iPr=",iPr,trim(vSpc(iAq)%NamSp)
      end if
      !
      Pivot= tStoik(vIPivot(K),K)
      tStoik(:,K)= tStoik(:,K) /Pivot
      !
      !iPr= iPr + 1
      !vCpn(iPr)%iSpc= iAq
      vIndex(iPr)=  iAq
      vCpDone(iPr)= .true.
      vSpDone(iAq)= .true.
      !
      exit DoAqu
      !
    end do DoAqu
    !----------------------------------------------/scan aqueous species
  end do Do00
  !
  if(iDebug==4) then
    write(71,'(A,1X)',advance="NO") "BASIS="
    do I=1,nCp
      write(71,'(A15,1X)',advance="NO") vSpc(vIndex(I))%NamSp
    end do
    write(71,*)
  end if
  !
  deallocate(vSpDone)
  deallocate(vCpDone)
  deallocate(vSpPrm)
  deallocate(tStoik)
  deallocate(vIPivot)
  !
  if(iDebug==4) then
    print *,"< Basis_Cpn_ISpc_Select"
    do i=1,nCp
      write(6,'(G15.6,1X,A)') &
      & vSpc(vIndex(I))%Dat%Mole, &
      & trim(vSpc(vIndex(I))%NamSp)
    end do
    print *,"</ Basis_Cpn_ISpc_Select"
    call Pause_
  end if
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Cpn_ISpc_Select"
  !
  return
end subroutine Basis_Cpn_ISpc_Select

subroutine Basis_Species_Select_Test( & !
& vSpc,isW,vCpn,tFormul, & !in
& vIndex) !out
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  use M_Numeric_Tools,only: CalcPermut
  use M_T_Species,    only: T_Species
  use M_T_Component,  only: T_Component
  !
  type(t_Species),  intent(in) :: vSpc(:)
  integer,          intent(in) :: isW
  type(t_Component),intent(in) :: vCpn(:)
  real(dp),         intent(in) :: tFormul(:,:)
  integer,          intent(out):: vIndex(:)
  !
  integer :: nCp,nAq,nCi
  integer :: i
  integer, allocatable:: vOrder(:)
  !
  nAq= count(vSpc%Typ=="AQU")
  nCp= size(vCpn)
  nCi= count(vCpn(:)%Statut=="INERT")
  !
  allocate(vOrder(nAq))
  !--------------------sort aqu.species in order of decreasing abudance,
  !------------------------------------- -> permutation array vSpvPrm --
  !-------- to apply permutation vSpPrm to array Arr: Arr-Arr(vSpPrm) --
  call CalcPermut(-vSpc(1:nAq)%Dat%Mole,vOrder)
  !--------------------/
  !
  vIndex(:)= 0
  vIndex(1)= isW
  vIndex(nCi+1:nCp)= vCpn(nCi+1:nCp)%iSpc
  !
  call Basis_Species_Select( &
  !& vSpc,isW,vCpn, & !in
  & tFormul,  & !in
  & vOrder,   & !in
  & vIndex)       !inout
  !! do i=2,nCp
  !!   if(vCpn(i)%Statut/="INERT") vIndex(i)= vCpn(i)%iSpc
  !! end do
  !
  !! if(iDebug==4) then
    print *,"< Basis_Cpn_ISpc_Select"
    do i=1,nAq
      write(6,'(G15.6,1X,A)') &
      & vSpc(vOrder(i))%Dat%Mole,trim(vSpc(vOrder(i))%NamSp)
    end do
    print *,"== NEW BASIS =="
    do i=1,nCp
      write(6,'(G15.6,1X,A)') &
      & vSpc(vIndex(I))%Dat%Mole,trim(vSpc(vIndex(I))%NamSp)
    end do
    print *,"</Basis_Cpn_ISpc_Select"
    call Pause_
  !! end if
  !
  deallocate(vOrder)
end subroutine Basis_Species_Select_Test

subroutine Basis_Species_Select( &
!& vSpc,isW,vCpn, & !in
& tFormulSpc, & !in
& vOrder,    &  !in
& vIndex)       !inout
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  !use M_T_Species,    only: T_Species
  !use M_T_Component,  only: T_Component
  !
  !type(t_Species),  intent(in) :: vSpc(:)
  !integer,          intent(in) :: isW
  !type(t_Component),intent(in) :: vCpn(:)
  real(dp),intent(in)   :: tFormulSpc(:,:)
  integer, intent(in)   :: vOrder(:)
  integer, intent(inout):: vIndex(:)
  !
  ! vIndex(1)= isW
  ! vIndex(nCi+1:nCp)= vCpn(nCi+1:nCp)%iSpc
  !
  integer :: nCp,nAq,nCx
  integer :: I,K
  logical :: Error
  real(dp),allocatable:: tStoikSpc(:,:)
  real(dp),allocatable:: tStoikCpn(:,:)
  !
  if(iDebug==4) write(fTrc,'(/,A)') "<-------------Basis_Species_Select"
  !
  nAq= size(vOrder)
  nCp= size(vIndex)
  nCx= count(vIndex(:)>0)
  !
  allocate(tStoikCpn(nCp,nCx))
  allocate(tStoikSpc(nCp,nAq))
  !
  !------------------------build the stoikio table of imposed components
  K= 0
  do i=1,nCp
    if(vIndex(i)>0) then
      K= K+1
      tStoikCpn(1:nCp,K)= tFormulSpc(1:nCp,vIndex(i))
    end if
  end do
  !-----------------------/build the stoikio table of imposed components
  do i=1,nAq
    tStoikSpc(1:nCp,i)= tFormulSpc(1:nCp,i)
  end do
  !
  call Basis_FreeSet_Select( & !
  & tStoikSpc,    & !IN
  & vIndex,       & !OUT
  & Error,        & !OUT
  & vOrder=      vOrder,    & !IN
  & tStoikioCpn= tStoikCpn)   !IN
  !
  deallocate(tStoikSpc)
  deallocate(tStoikCpn)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</-------------Basis_Species_Select"
  !
  return
end subroutine Basis_Species_Select

subroutine Basis_FreeSet_Select( & !
& tStoikio,    & !IN
& vIndex,      & !INOUT
& Error,       & !OUT
& vOrder,      & !IN
& tStoikioCpn)   !IN
!--
!-- using Gaussian elimination (without permutations),
!-- find a set of independent species (of stoikios tStoikio)
!-- independent also of those in tStoikioCpn
!--
  use M_Numeric_Tools,only: iMaxLoc_R
  !
  real(dp),intent(in)   :: tStoikio(:,:)
  integer, intent(inout):: vIndex(:)
  logical, intent(out)  :: Error
  integer, intent(in),optional:: vOrder(:)
  real(dp),intent(in),optional:: tStoikioCpn(:,:)
  !
  integer :: nCp,nSp,nX
  integer :: iSp,I,J,K,N
  real(dp):: Y,Pivot
  real(dp),allocatable:: tStoikTmp(:,:)
  integer, allocatable:: vIPivot(:)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "<--------------Basis_FreeSet_Select"
  !
  Error= .false.
  !
  if(present(tStoikioCpn)) then
    nX= size(tStoikioCpn,2)
  else
    nX= 0
  end if
  nCp= size(tStoikio,1)
  nSp= size(tStoikio,2)
  !
  allocate(tStoikTmp(nCp,nCp))   ;  tStoikTmp(:,:)= Zero
  allocate(vIPivot(nCp))         ;  vIPivot(:)= 0
  !
  K= 0
  !----------------------------------------scan imposed basis components
  !--e.g. SOLVENT / BUFFER / MOBILE components--------------------------
  !--(they have vIndex(:)>0)--------------------------------------------
  do i=1,nX
    !
    !-- read stoikiometries of species already imposed, in order to select,
    !-- in do_0 block,
    !-- primary aqueous species independent from these species
    !
    K= K+1
    tStoikTmp(:,K)= tStoikioCpn(:,i)
    do J=1,K-1
      Y= tStoikTmp(vIPivot(J),K)
      tStoikTmp(:,K)= tStoikTmp(:,K)- Y*tStoikTmp(:,J)
      !-> substract column vIPivot(J) from column K
      !-> tStoikTmp(vIPivot(J),K) will be 0
    end do
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(:,K)))
    Pivot= tStoikTmp(vIPivot(K),K)
    if(ABS(Pivot)<1.D-9) then
      != all coeff's are 0 -> this species is not independent -> exit ==
      Error= .true.
      return
    end if
    !
    tStoikTmp(:,K)= tStoikTmp(:,K) /Pivot
    !-> first non-zero element, tStoikTmp(vIPivot(K),K), will be 1.
    !
  end do
  !---------------------------------------/scan imposed basis components
  !
  N=1
  !-----------------------------------------------scan the other species
  do_0: do iSp=1,nSp
    !
    K= K+1
    if(K>nCp) exit do_0
    !
    if(present(vOrder)) then ;  tStoikTmp(:,K)= tStoikio(:,vOrder(iSp))
    else                     ;  tStoikTmp(:,K)= tStoikio(:,iSp)
    end if
    !
    do J=1,K-1
      Y= tStoikTmp(vIPivot(J),K)
      tStoikTmp(:,K)= tStoikTmp(:,K)- Y*tStoikTmp(:,J)
      !-> tStoikTmp(vIPivot(J),K) will be 0
    end do
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(:,K)))
    Pivot= tStoikTmp(vIPivot(K),K)
    if(ABS(Pivot)<1.D-9) then
      !-- all coeff's are 0 -> this species is not independent -> cycle do_0 --
      ! if(iDebug>2) then
      !   do JJ=1,nCp; write(6,'(F7.2,1X)',advance="NO") tStoikTmp(JJ,K); end do
      !   write(6,'(2(A,I3),1X,A)') "K=",K," /  iPr=",K,"Reject"
      ! end if
      K= K-1
      cycle do_0
    end if
    !
    ! if(iDebug>2) then
    !   do JJ=1,nCp; write(6,'(F7.2,1X)',advance="NO") tStoikTmp(JJ,K); end do
    !   write(6,'(2(A,I3),1X,A)') "K=",K," /  iPr=",K,"Accept"
    ! end if
    !
    tStoikTmp(:,K)= tStoikTmp(:,K) /Pivot
    !
    !--TRICK: inert aqueous species are placed--------------------------
    !--AFTER solvent and BEFORE mobile species !!!----------------------
    ! vIndex(K-nX+1)= vOrder(iSp)
    do while(vIndex(N)>0)
      N= N+1
    end do
    if(present(vOrder)) then ;  vIndex(N)= vOrder(iSp)
    else                     ;  vIndex(N)= iSp
    end if
    !
  end do do_0
  !----------------------------------------------/scan the other species
  !
  deallocate(tStoikTmp)
  deallocate(vIPivot)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</-------------Basis_FreeSet_Select"
  !
  return
end subroutine Basis_FreeSet_Select

subroutine Basis_AlfaNu_Build( & !
& vEle,vCpn,vSpc,tStoikio,     & !IN
& vFas,vMixFas,vMixModel,      & !
& vLCi,vLAx,vLMx,vLAs,vLMs,    & !
& LSho,LBuffer,                & !
& vOrdPr,vOrdAs,vOrdMs,        & !
& vTot,                        & !
& tAlfSp,                      & !OUT
& tAlfPr,tAlfAs,tAlfMs,        & !
& tNuSp,tNuAs,tNuMs,           & !
& tNuFas,tAlfFs)                 !

  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  use M_T_Phase,    only: T_Phase
  use M_T_MixPhase, only: T_MixPhase
  use M_T_MixModel, only: T_MixModel
  !
  type(T_Element),  dimension(:),  intent(in):: vEle
  type(T_Component),dimension(:),  intent(in):: vCpn
  type(T_Species),  dimension(:),  intent(in):: vSpc
  real(dp),         dimension(:,:),intent(in):: tStoikio
  type(T_Phase),    dimension(:),  intent(in):: vFas
  type(T_MixPhase), dimension(:),  intent(in):: vMixFas
  type(T_MixModel), dimension(:),  intent(in):: vMixModel
  logical,          dimension(:),  intent(in):: vLCi,vLAx,vLMx,vLAs,vLMs
  logical,                         intent(in):: LSho,LBuffer
  integer,          dimension(:),  intent(in):: vOrdPr,vOrdAs,vOrdMs
  !
  real(dp),         dimension(:),  intent(inout):: vTot
  real(dp),         dimension(:,:),intent(out)::   tAlfSp !,tAlfEle
  real(dp),         dimension(:,:),intent(out)::   tAlfPr,tAlfAs,tAlfMs
  real(dp),         dimension(:,:),intent(out)::   tNuSp,tNuAs,tNuMs
  real(dp),         dimension(:,:),intent(out)::   tNuFas,tAlfFs
  !
  integer:: nCi,nAs,nMs,nAx,nMx
  integer:: I,J
  integer:: vIndxPrimSpc(size(vCpn))
  integer:: vIndxBuffers(size(vCpn))
  !
  if(iDebug>0) write(fTrc,'(/,A)') "<----------------Basis_AlfaNu_Build"
  !
  nCi= count(vLCi)
  nAx= count(vLAx)
  nMx= count(vLMx)
  nAs= count(vLAs)
  nMs= count(vLMs)
  !
  vIndxBuffers(:)= 0  !
  if(LBuffer) then
    do I=1,size(vCpn)
      if(vCpn(I)%Statut=="BUFFER") vIndxBuffers(I)= vCpn(I)%iSpc
      !if(vCpn(I)%Statut=="BUFFER" .or. vCpn(I)%Statut=="MOBILE") &
      !& vIndxBuffers(I)= vCpn(I)%iSpc
    end do
  end if
  !
  call AlfaTable_Build( &
  & size(vCpn),size(vSpc),tStoikio, & !
  & vIndxBuffers, & !
  & vTot(:), & !
  & tAlfSp) !,tAlfEle)
  !
  call Fas_AlfaTable_Build( & !
  & vFas,vMixFas,vMixModel,tAlfSp, & !IN
  & tAlfFs)   !OUT
  !
  call AlfaTable_Extract( &
  & vCpn,tAlfSp,vLCi,vLAx,vLMx,vLAs,vLMs, &
  & tAlfPr,tAlfAs,tAlfMs)
  !
  !--------- component order in vCpn is Inert, Aqu'Mobile, Min'Mobile --
  J= 0
  !-------------------------------------- index of inert prim'species --
  do I=1,nCi  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(I)%iSpc  ;  end do
  !-------------------------------------- index of mobile aqu'species --
  do I=1,nAx  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(J)%iSpc  ;  end do
  !-------------------------------------- index of mobile min'species --
  do I=1,nMx  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(J)%iSpc  ;  end do
  !
  call NuTable_Build( & !
  & vIndxPrimSpc,tStoikio, & !IN
  & tNuSp) !OUT
  !
  call NuTable_Extract( & !
  & tNuSp,vLAs,vLMs, & !IN
  & tNuAs,tNuMs) !OUT
  !
  call Fas_NuTable_Build( & !
  & vFas,vMixFas,vMixModel,tNuSp, & !IN
  & tNuFas)   !OUT
  !
  if(LSho) call AlfaTable_Sho( &
  & vEle,vCpn,vSpc,vFas,vOrdPr,vOrdAs,vOrdMs, &
  & tAlfPr,tAlfAs,tAlfMs,tAlfFs)
  !
  if(LSho) call NuTable_Sho( &
  & vCpn,vSpc,vLCi,vLAx,vLMx,vLAs,vLMs, &
  & vOrdPr,vOrdAs,vOrdMs,tNuSp,tNuAs,tNuMs, &
  & .true.)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</---------------Basis_AlfaNu_Build"
  !
end subroutine Basis_AlfaNu_Build

subroutine Basis_StoikTable_Sho( &
& vEle,vCpn,vSpc, &         !in
& tStoikio, &               !in
& vLCi,vLAx,vLMx,vLAs,vLMs) !in

  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut,Files_Index_Write
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  !
  type(T_Element),  intent(in):: vEle(:)
  type(T_Component),intent(in):: vCpn(:)
  type(T_Species),  intent(in):: vSpc(:)
  real(dp),         intent(in):: tStoikio(:,:)
  logical,          intent(in):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !
  integer:: F,I,iPr,nCp,nSp
  !
  nSp=size(vSpc)
  nCp=size(vCpn)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_stoik_for.log")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_stoik_for.log",&
  & "STOIKIO: stoichiometry matrix tStoikio, transposed")
  !
  !--- header
  do iPr=1,nCp
    write(F,'(A3,A1)',advance="no") vEle(vCpn(iPr)%iEle)%NamEl, T_
  end do
  write(F,'(3(A,A1))') "Chg",T_,"_____species",T_,"STATUS",T_
  !---/
  !
  do I=1,nSp
    do iPr=1,nCp
      write(F,'(F7.2,A1)',advance="no") tStoikio(iPr,I), T_
    end do
    write(F,'(I3,A1,A12,A1)',advance="no") vSpc(I)%Z, T_, trim(vSpc(I)%NamSp),T_
    if(vLCi(I)) write(F,'(A)') "INT_PRIMARY"
    if(vLAx(I)) write(F,'(A)') "EXT_AQUEOUS"
    if(vLMx(I)) write(F,'(A)') "EXT_MINERAL"
    if(vLAs(I)) write(F,'(A)') "INT_AQUEOUS"
    if(vLMs(I)) write(F,'(A)') "INT_MINERAL"
  end do
  !
  close(F)
  !
end subroutine Basis_StoikTable_Sho

subroutine NuTable_Sho( &
& vCpn,vSpc, &
& vLCi,vLAx,vLMx,vLAs,vLMs, &
& vOrdPr,vOrdAs,vOrdMs, &
& tNuSp,tNuAs,tNuMs, &
& LogK_Sho)
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut,Files_Index_Write
  use M_Numeric_Const,  only: Ln10
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  !
  type(T_Component),dimension(:),  intent(in):: vCpn
  type(T_Species),  dimension(:),  intent(in):: vSpc
  logical,          dimension(:),  intent(in):: vLCi,vLAx,vLMx,vLAs,vLMs
  integer,          dimension(:),  intent(in):: vOrdPr,vOrdAs,vOrdMs
  real(dp),         dimension(:,:),intent(in):: tNuSp,tNuAs,tNuMs
  logical,                         intent(in):: LogK_Sho
  !
  integer :: F,iPr,I
  integer :: nSp,nCp,nAs,nMs
  real(dp):: X,rScale
  !
  nSp=size(vSpc)
  nCp=size(vCpn)
  nAs=count(vLAs)
  nMs=count(vLMs)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_stoik_nu.log") !"OUT_STOIK2.TXT")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_stoik_nu.log",&
  & "STOIKIO: Nu Table: stoikio of second'species in terms of prim'species")
  !
  write(F,'(2A)') &
  & "!.tNuSp, ALL SPECIES", &
  & "!.stoichiometry, logK_dtb, logK_reaction, logK_reaction/scale"
  !
  do iPr=1,nCp
    write(F,'(A7,A1)',advance="no") trim(vSpc(vOrdPr(iPr))%NamSp),T_
  end do
  !
  write(F,'(A)') "Status and Name(iSpc)"
  !
  do I=1,nSp
    !
    do iPr=1,nCp
      if(tNuSp(I,iPr)/=0) then; write(F,'(F7.2,A1)',advance="no") tNuSp(I,iPr), T_
      else;                     write(F,'(A7,  A1)',advance="no")          ".", T_
      end if
    end do
    !
    if(vLCi(I)) write(F,'(A,A1)',advance="no") "INT_PRIMARY",T_
    if(vLAx(I)) write(F,'(A,A1)',advance="no") "EXT_AQUEOUS",T_
    if(vLMx(I)) write(F,'(A,A1)',advance="no") "EXT_MINERAL",T_
    if(vLAs(I)) write(F,'(A,A1)',advance="no") "INT_AQUEOUS",T_
    if(vLMs(I)) write(F,'(A,A1)',advance="no") "INT_MINERAL",T_
    !
    if(LogK_Sho) then
      X= vSpc(I)%G0rt - dot_product(tNuSp(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
      rScale= SUM(ABS(tNuSp(I,1:nCp))) !-> scaling factor !!!
      write(F,'(G15.6,A1)',advance="no") rScale,T_
      if(rScale==Zero) rScale= One
      write(F,'(3(F15.6,A1))',advance="no") &
      & -vSpc(I)%G0rt /Ln10,T_,X/Ln10,T_,X/Ln10/rScale,T_
      != +log_10(K of Formation Reaction)
    end if
    !
    write(F,'(A)') trim(vSpc(I)%NamSp)
    !
  end do
  !
  if(nAs>0) then

    write(F,'(/,A,/)') "!.tNuAs, SECOND.AQU, logK_dtb, logK_reaction"

    do iPr=1,nCp
      write(F,'(A7,A1)',advance="no") trim(vSpc(vOrdPr(iPr))%NamSp),T_
    end do

    write(F,'(A)') "Name(iAs)"

    do I=1,nAs
      do iPr=1,nCp !---------------------------------------write stoikio
        write(F,'(F7.2,A1)',advance="no") tNuAs(I,iPr), T_
      end do
      !
      if(LogK_Sho) then
        X=  vSpc(vOrdAs(I))%G0rt &
        & - dot_product(tNuAs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
        !X=vG0As(I) - dot_product(tNuAs(I,1:nCp),vG0Pr(1:nCp))
        write(F,'(2(F15.6,A1))',advance="no") &
        & -vSpc(vOrdAs(I))%G0rt /Ln10,T_,X/Ln10,T_ !=+log10K of Formation Reaction
      end if
      !
      write(F,'(A)') trim(vSpc(vOrdAs(I))%NamSp)
      !
    end do

  end if
  !
  if(nMs>0) then

    write(F,'(/,A,/)') "!.tNuMs, SECOND.MIN, logK_dtb, logK_reaction"

    do iPr=1,nCp
      write(F,'(A7,A1)',advance="no") trim(vSpc(vOrdPr(iPr))%NamSp),T_
    end do

    write(F,'(A)') "Name(iMs)"

    do I=1,nMs
      !
      do iPr=1,nCp !===================================< write stoikio==
        write(F,'(F7.2,A1)',advance="no") tNuMs(I,iPr), T_
      end do
      !
      if(LogK_Sho) then
        X=  vSpc(vOrdMs(I))%G0rt &
        & - dot_product(tNuMs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
        !X=vG0As(iAs) - dot_product(tNuAs(iAs,1:nCp),vG0Pr(1:nCp))
        write(F,'(2(F15.6,A1))',advance="no") &
        & -vSpc(vOrdMs(I))%G0rt /Ln10,T_,X/Ln10,T_ !=+log10K of Formation Reaction
      end if
      !
      write(F,'(A)') trim(vSpc(vOrdMs(I))%NamSp)
      !
    end do

  end if
  !
  close(F)
end subroutine NuTable_Sho

subroutine AlfaTable_Sho( &
& vEle,vCpn,vSpc,vFas,    &
& vOrdPr,vOrdAs,vOrdMs,   &
& tAlfPr,tAlfAs,tAlfMs,tAlfFs)
  !---------------------------------------------------------------------
  use M_IoTools,    only: GetUnit
  use M_Files,      only: DirOut,Files_Index_Write
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  use M_T_Phase,    only: T_Phase
  !
  type(T_Element),  dimension(:),  intent(in):: vEle
  type(T_Component),dimension(:),  intent(in):: vCpn
  type(T_Species),  dimension(:),  intent(in):: vSpc
  type(T_Phase),    dimension(:),  intent(in):: vFas
  integer,          dimension(:),  intent(in):: vOrdPr,vOrdAs,vOrdMs
  real(dp),         dimension(:,:),intent(in):: tAlfPr,tAlfAs,tAlfMs,tAlfFs
  !
  integer:: nCp,nCi,nAs,nMs,nFs
  integer:: F
  integer:: I,J
  !---------------------------------------------------------------------
  nCp= size(vCpn)
  nFs= size(vFas)
  nCi= count(vCpn(:)%Statut=="INERT") + count(vCpn(:)%Statut=="BUFFER")
  nAs= size(tAlfAs,2)
  nMs= size(tAlfMs,2)
  !nAs=count(vLAs)
  !nMs=count(vLMs)
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_stoik_alfa.log") !"out_stoik1.txt")
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_stoik_alfa.log",&
  & "STOIKIO: Alfa Tables: stoikio of all species in terms of elements / prim'species")
  !
  write(F,'(/,A,/)') "tAlfPr, Prim'Species"
  call NamCp_Sho
  do I=1,nCp
    do J=1,nCp; write(F,'(F7.3,A1)',advance="no") tAlfPr(J,I), T_; end do
    write(F,'(I3,A1,A12)') vSpc(vOrdPr(I))%Z, T_, trim(vSpc(vOrdPr(I))%NamSp)
  end do
  !
  if(nAs>0) then
    write(F,'(/,A,/)') "tAlfAs, Aqu'Sec'Species"
    call NamCp_Sho
    do I=1,nAs
      do J=1,nCp; write(F,'(F7.3,A1)',advance="no") tAlfAs(J,I), T_; end do
      write(F,'(I3,A1,A12)') vSpc(vOrdAs(I))%Z, T_, trim(vSpc(vOrdAs(I))%NamSp)
    end do
  end if
  if(nMs>0) then
    write(F,'(/,A,/)') "tAlfMs, Min'Sec'Species"
    call NamCp_Sho
    do I=1,nMs
      do J=1,nCp; write(F,'(F7.3,A1)',advance="no") tAlfMs(J,I), T_; end do
      write(F,'(I3,A1,A12)') vSpc(vOrdMs(I))%Z, T_, trim(vSpc(vOrdMs(I))%NamSp)
    end do
  end if
  ! if(nFs>0) then
  !   write(F,'(/,A,/)') "tAlfFs, Phases"  
  !   do J=1,nCi;     write(F,'(A3,A1)',advance="no") vEle(vCpn(J)%iEle)%NamEl, T_; end do
  !   do J=nCi+1,nCp; write(F,'(A7,A1)',advance="no") vEle(vCpn(J)%iEle)%NamEl, T_; end do
  !   write(F,*)
  !   do I=1,nFs
  !     do J=1,nCp; write(F,'(F7.3,A1)',advance="no") tAlfFs(J,I), T_; end do
  !     write(F,'(A)') trim(vFas(I)%NamFs)
  !   end do
  ! end if
  close(F)

contains

subroutine NamCp_Sho
  integer:: K
  do K=1,nCi;     write(F,'(A3,A1)',advance="no") vCpn(K)%NamCp, T_; end do
  do K=nCi+1,nCp; write(F,'(A7,A1)',advance="no") vCpn(K)%NamCp, T_; end do
  write(F,*)
end subroutine NamCp_Sho

end subroutine AlfaTable_Sho

subroutine Basis_FindIBal( & !
& vCpn,vSpc,               & ! IN
& iH_,iBal)                  ! OUT
!--
!-- find index of the BALANCE element -> iBal
!--
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species
  !
  type(T_Component),intent(in):: vCpn(:)
  type(T_Species),  intent(in):: vSpc(:)
  integer,          intent(in) :: iH_
  integer,          intent(out):: iBal
  !---------------------------------------------------------------------
  integer::I
  !
  !-------------- finds for which element the mole balance constraint --
  !------------------------------ is replaced by a charge balance one --
  iBal= 0
  do I=1,size(vCpn)
    if(trim(vCpn(I)%Statut)=="BALANCE") then
      iBal=I
      exit
    end if
  end do
  !
  !<new, deleted 200911>
  !~ if(iBal>0) then
    !~ if (vSpc(vCpn(iBal)%iSpc)%Z==0) &
    !~ & call Stop_("Balance Species MUST BE Charged Species")
  !~ end if
  !</new>
  !
  if(iBal==0 .and. vCpn(iH_)%Statut=="INERT") iBal= iH_
  !
  !---------------------------- option ELECTRONEUTRALITY NOT ENFORCED --
  !=< used mainly for construction of pH diagram using a PATH CHANGE, ==
  !----------------------------------- not for speciation of a system --
  if(iBal==0 &
  & .and. vCpn(iH_)%Statut/="INERT" &
  & .and. iDebug>0) &
  & print '(A)', &
  & "CAVEAT: NO Element For Electron Balance -> ELECTRONEUTRALITY NOT ENFORCED"
  !! & call Stop_("NO Element For Electron Balance ????")
  !
end subroutine Basis_FindIBal

end module M_Basis_Tools

!~ subroutine AlfaTable_Calc( &
!~ & vCpn,vSpc,tStoikio,vLCi,vLAx,vLMx,vLAs,vLMs,bOH2, &
!~ & tAlfSp,tAlfPr,tAlfAs,tAlfMs)
!~ !.called in Basis_AlfaNu_Build
  !~ use M_T_Component,only: T_Component
  !~ use M_T_Species,  only: T_Species
  !~ use M_Numeric_Mat,      only: LU_Decomp, LU_BakSub
  !~ !
  !~ type(T_Component),intent(in):: vCpn(:)
  !~ type(T_Species),  intent(in):: vSpc(:)
  !~ real(dp), intent(in):: tStoikio(:,:)
  !~ logical,  intent(in):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !~ logical,  intent(in):: bOH2
  !~ !
  !~ real(dp), intent(out):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !~ !
  !~ real(dp),allocatable:: tTransform(:,:)
  !~ integer, allocatable:: Indx(:)
  !~ real(dp),allocatable:: Y(:)
  !~ !
  !~ logical :: bSingul
  !~ !
  !~ real(dp):: D
  !~ integer :: I,nCp,nSp,nCi,nAx,nMx,nAs,nMs
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A)') "< AlfaTable_Calc"
  !~ !
  !~ nCp= size(vCpn)
  !~ nSp= size(vSpc)
  !~ nCi= count(vLCi)
  !~ nAs= count(vLAs)
  !~ nAx= count(vLAx)
  !~ nMx= count(vLMx)
  !~ nMs= count(vLMs)
  !~ !
  !~ allocate(tTransform(nCp,nCp))
  !~ allocate(Indx(nCp))
  !~ allocate(Y(nCp))
  !~ !
  !~ !tStoikio:column iSp = stoikio of species iSp in terms of elements (1:nCp)
  !~ !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !~ !
  !~ !extract blocks from tStoikio(1:nCp,1:nSp) to tAlfPr,tAlfAs,tAlfMs
  !~ !
  !~ tAlfSp(1:nCp,1:nSp)= tStoikio(1:nCp,1:nSp)
  !~ !
  !~ !do I=1,nCp
  !~ !  tTransform(:,I)= vCpn(I)%vStoikio(1:nCp)
  !~ !end do
  !~ !
  !~ if(bOH2) then
    !~ tTransform(1:nCp,1:nCp)=Zero
    !~ do I=1,nCp; tTransform(I,I)=One; end do !-> identity matrix
    !~ tTransform(1:nCp,1)= tStoikio(1:nCp,1) !-> column 1, O replaced by OH2
    !~ call LU_Decomp(tTransform,Indx,D,bSingul)
    !~ if(bSingul) call Stop_("SINGUL IN AlfaTableCalc")
    !~ !-> tTransform gives Elements as Functions of Int'Elements/Ext'Species
    !~ !
    !~ do I=1,nSp
      !~ Y= tStoikio(1:nCp,I)
      !~ call LU_BakSub(tTransform,Indx,Y)
      !~ tAlfSp(1:nCp,I)=Y(1:nCp)
    !~ end do
    !~ !
  !~ end if
  !~ !
  !~ deallocate(tTransform)
  !~ deallocate(Indx)
  !~ deallocate(Y)
  !~ !
  !~ !__________________________________________subblock of inert prim'species
  !~ do I=1,nCi
    !~ tAlfPr(:,I)=tAlfSp(:,vCpn(I)%iSpc)
  !~ end do
  !~ !__________________________________________subblock of mobile aqu'species
  !~ do I=1,nAx
    !~ tAlfPr(:,nCi+I)=tAlfSp(:,vCpn(nCi+I)%iSpc)
  !~ end do
  !~ !__________________________________________subblock of mobile min'species
  !~ do I=1,nMx
    !~ tAlfPr(:,nCi+nAx+I)=tAlfSp(:,vCpn(nCi+nAx+I)%iSpc)
  !~ end do
  !~ !
  !~ if(nAs>0) call Extract(2,tAlfSp,vLAs,tAlfAs) !extract columns of sec'aqu.species
  !~ if(nMs>0) call Extract(2,tAlfSp,vLMs,tAlfMs) !extract columns of sec'min.gas.species
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A,/)') "</ AlfaTable_Calc"
  !~ !
!~ end subroutine AlfaTable_Calc

!~ subroutine FasTable_Calc( & !
!~ & vCpn,vFas,vMixFas,vMixModel, & !in
!~ & tNuSp, tAlfSp, & !in
!~ & tNuFas,tAlfFs)   !out

  !~ use M_T_Component,only: T_Component
  !~ use M_T_Phase,    only: T_Phase
  !~ use M_T_MixPhase, only: T_MixPhase
  !~ use M_T_MixModel, only: T_MixModel
  !~ !
  !~ type(T_Component), intent(in) :: vCpn(:)
  !~ type(T_Phase),     intent(in) :: vFas(:)
  !~ type(T_MixPhase),  intent(in) :: vMixFas(:)
  !~ type(T_MixModel),  intent(in) :: vMixModel(:)
  !~ real(dp),          intent(in) :: tNuSp(:,:), tAlfSp(:,:)
  !~ !
  !~ real(dp),          intent(out):: tNuFas(:,:),tAlfFs(:,:)
  !~ !
  !~ type(T_MixModel):: SM
  !~ type(T_MixPhase):: SF
  !~ integer:: I,J,K,iCp,NPole !,nFs
  !~ integer:: nCp
  !~ !
  !~ nCp=size(vCpn)
  !~ !nFs=size(vFas)
  !~ !
  !~ do I=1,size(vFas)
    !~ !
    !~ select case(vFas(I)%Typ)
      !~ !
      !~ case("PURE","DISCRET")
        !~ K=vFas(I)%iSpc
        !~ tNuFas(I,1:nCp)= tNuSp(K,1:nCp)
        !~ tAlfFs(1:nCp,I)= tAlfSp(1:nCp,K)
      !~ !
      !~ case("MIXT")
        !~ !
        !~ SF=vMixFas(vFas(I)%iMix) !-> the solution phase
        !~ SM=vMixModel(SF%iModel)  !-> the corresponding model
        !~ NPole=SM%NPole
        !~ !
        !~ where(.not.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
        !~ !
        !~ do iCp=1,nCp
          !~ tNuFas(I,iCp)= Zero
          !~ tAlfFs(iCp,I)= Zero
          !~ !
          !~ do J=1,NPole
            !~ if(SF%vLPole(J)) then
              !~ tNuFas(I,iCp)= tNuFas(I,iCp) &
              !~ &            + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
              !~ tAlfFs(iCp,I)= tAlfFs(iCp,I) &
              !~ &            + tAlfSp(iCp,SM%vIPole(J))*SF%vXPole(J)
            !~ end if
          !~ end do
          !~ !
        !~ end do
      !~ !
    !~ end select
  !~ end do
  !~ !
!~ end subroutine FasTable_Calc

!~ subroutine Basis_Calc_iW_iH_iOH_iO2( & !
!~ & vSpc, &                  !in
!~ & isW,isH_,isOH,isO2,MWsv) !out
!~ !--
!~ !-- find indexes of species H+, OH-, O2 in vSpc (unchanged for given vSpc)
!~ !--
  !~ use M_T_Species,only: T_Species,Species_Index
  !~ !
  !~ type(T_Species),intent(in) :: vSpc(:)
  !~ integer,        intent(out):: isW,isH_,isOH,isO2
  !~ real(dp),       intent(out):: MWsv
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A)') "< Basis_Calc_iW_iH_iOH_iO2"
  !~ !
  !~ !-------------- ranks of species H+ and OH- in current species list --
  !~ isW=  Species_Index("H2O",vSpc)
  !~ isOH= Species_Index("OH-",  vSpc); if(isOH==0) isOH=Species_Index("OH[-]",vSpc)
  !~ isH_= Species_Index("H+",   vSpc); if(isH_==0) isH_=Species_Index("H[+]",vSpc)
  !~ isO2= Species_Index("O2,AQ",vSpc); if(isO2==0) isO2=Species_Index("O2(AQ)",vSpc)
  !~ isO2= Species_Index("O2_AQ",vSpc)
  !~ !
  !~ if(isW--0)  call Stop_("species H2O not found !!!") !-------------stop
  !~ if(isH_--0) call Stop_("species H+ not found  !!!") !-------------stop
  !~ if(isOH--0) call Stop_("species OH- not found !!!") !-------------stop
  !~ !
  !~ MWsv=vSpc(isW)%WeitKg
  !~ !
  !~ if(iDebug>0) write(fTrc,'(2(A,1X,I3))') "  isH_",isH_,"  isOH",isOH
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A,/)') "</ Basis_Calc_iW_iH_iOH_iO2"
!~ end subroutine Basis_Calc_iW_iH_iOH_iO2

!~ subroutine Basis_RedoxAssign(vCpn,vEle,vSpc)
!~ !------------------------------------ assign redox states of elements --
  !~ use M_T_Element,  only: T_Element
  !~ use M_T_Component,only: T_Component
  !~ use M_T_Species,  only: T_Species
  !~ !
  !~ type(T_Component),intent(in) :: vCpn(:)
  !~ type(T_Element),  intent(in) :: vEle(:)
  !~ type(T_Species),  intent(in) :: vSpc(:)
  !~ !
  !~ type(T_Species):: S_
  !~ integer:: nCp,I, ZSp, Z_
  !~ logical:: fOk
  !~ !
  !~ nCp=size(vCpn)
  !~ !
  !~ do I=1,nCp
    !~ S_=vSpc(vCpn(I)%iSpc)
    !~ !
    !~ !calculate redox state of element in its Prim'Species -> basis valency
    !~ Z_= dot_product( S_%vStoikio(1:nCp),vEle(1:nCp)%Z ) &
    !~ & - S_%vStoikio(I)*vEle(I)%Z
    !~ !
    !~ !restrict this procedure to ele's with no valency assigned in dtb (Fe,S,...)
    !~ if(vEle(I)%Redox=="VAR") then
      !~ !vEle(I)%Z=(ZSp-Z_)/vStoik(I)
      !~ if(S_%vStoikio(I)>0) &
      !~ write(fTrc,'(2A,2(A,I3),/)') &
      !~ & "valency of ",vEle(I)%NamEl,&
      !~ & ", Old=",vEle(I)%Z,&
      !~ & ", New=",(ZSp-Z_) /S_%vStoikio(I)
    !~ end if
    !~ !
  !~ end do
  !~ !
!~ end subroutine Basis_RedoxAssign

!~ subroutine Basis_Calc_ieO_ieH_ieOx(vEle,ieO_,ieH_,ieOx)
!~ != find indexes of elements H and OX in vEle
!~ != -> "static" (stable for a given database)
!~ != used for an aqueous system
  !~ use M_T_Element,only: T_Element,Element_Index
  !~ !
  !~ type(T_Element),dimension(:),intent(in) :: vEle
  !~ integer,                     intent(out):: ieO_,ieH_,ieOx
  !~ !
  !~ ieO_= Element_Index("O__",vEle)
  !~ ieH_= Element_Index("H__",vEle)
  !~ ieOx= Element_Index("OX_",vEle)
  !~ !
!~ end subroutine Basis_Calc_ieO_ieH_ieOx

