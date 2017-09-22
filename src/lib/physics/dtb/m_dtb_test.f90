module M_Dtb_Test
!--
!-- routines or testing databases, solution models, etc
!--
  use M_Kinds
  use M_T_Tpcond,only: T_TPCond
  use M_Trace,   only: iDebug,fTrc,fHtm,Stop_,T_,Pause_,Warning_

  implicit none
  !
  private
  !
  public:: Dtb_Tabulate
  public:: Dtb_Tabulate_ForSystem
  public:: Dtb_Test_Fluid
  public:: Dtb_Test_H2OHkf
  public:: Dtb_Test_AquHkf
  public:: Dtb_Test_Species
  public:: Dtb_Tabulate_PSatH2O
  public:: Dtb_Test_EQ36
  public:: Dtb_Test_Mixture_Gmix
  !
contains

subroutine Dtb_Test_Mixture_Gmix
  use M_Dtb_Const,   only: T_CK
  use M_IoTools,     only: GetUnit
  use M_TPcond_Read, only: TPpath_Read
  use M_T_MixModel,  only: T_MixModel, MixModel_GibbsMixRT 
  use M_Global_Tools,only: Species_TP_Update
  !
  use M_Global_Vars,only: vSpcDtb,vSpc,vMixModel
  use M_Global_Vars,only: vDiscretModel,vDiscretParam
  use M_Path_Vars,  only: vTPpath
  !
  type(T_MixModel):: MM
  real(dp),allocatable:: vX(:)
  logical, allocatable:: vL(:)
  real(dp):: TdgK, Pbar
  real(dp):: Gmix
  integer :: i, F, NP
  !
  call GetUnit(F)
  open(F,file="test_mixmodel.tab")
  !
  MM= vMixModel(1)
  NP= MM%NPole
  allocate(vX(NP))
  allocate(vL(NP))
  !
  vX(:)= 0.D0
  vL(:)= .false.
  vL(1:3)= .true.
  !
  vX(1)= 0.D0
  vX(2)= 0.01D0
  do
    vX(1)= vX(1)+1.D-2
    vX(3)= 1.D0 - vX(1)-vX(2)
    write(F,'(G15.6,A1,$)') vX(1),T_
    if(vX(1)>0.999D0) exit
  end do
  write(F,*)
  !
  !--read the TP path
  call TPpath_Read
  !
  do i=1,size(vTPpath)
    print *,TdgK
    !
    TdgK= vTPpath(i)%TdgC +T_CK
    Pbar= vTPpath(i)%Pbar
    !
    call Species_TP_Update( &
    & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
    & vSpc,vMixModel) !inout
    !
    vX(1)= 0.D0
    vX(3)= 0.01D0
    do
      !
      vX(1)= vX(1)+1.D-2
      vX(2)= 1.D0 - vX(1)-vX(3)
      if(vX(1)>0.999D0) exit
      !
      GMix= MixModel_GibbsMixRT( & !
      & TdgK,Pbar,  & !IN
      & MM,         & !IN, mixing model
      & vL,         & !IN
      & vX)           !IN
      !! GMix= GMix *R_jk*TdgK
      write(F,'(G15.6,A1,$)') GMix,T_
      !
    end do
    write(F,*)
    !
  end do
  !
  deallocate(vX)
  deallocate(vL)
  !
  close(F)
  !
  return
end subroutine Dtb_Test_Mixture_Gmix

subroutine Dtb_Tabulate(sCode)
!-----------------------------------------------------------------------
!--tabulate the species logK's along the vTPpath
!--retrieved from the TP.PATH block
!-----------------------------------------------------------------------
  use M_Dtb_Const,   only: T_CK
  use M_T_SolModel,  only: T_SolModel,T_SolModelDat
  use M_SolModel_Tools,only: SolModel_TP_Update
  use M_Dtb_Calc,    only: DtbSpc_Table_Build
  use M_TPcond_Read, only: TPpath_Read
  use M_Path_Vars,   only: vTPpath
  !
  use M_Global_Vars,only: vSpcDtb,vSpc
  !
  character(len=*),intent(in):: sCode
  !
  real(dp),allocatable:: tGrt(:,:),tVol(:,:)
  type(T_SolModelDat),allocatable:: vSolvDat(:)
  type(T_SolModel):: Slv
  integer :: i,N,M
  real(dp):: TdgK,Pbar

  !--read the TP path
  call TPpath_Read
  !
  N= size(vTPpath)
  M= size(vSpc)
  !
  allocate(vSolvDat(1:N))
  allocate(tGrt(1:size(vSpc),1:N))
  tGrt= Zero
  allocate(tVol(1:size(vSpc),1:N))
  tVol= Zero
  !
  !----------compute the solvent properties for the series of T,P points
  do i=1,N

    TdgK= vTPpath(i)%TdgC +T_CK
    Pbar= vTPpath(i)%Pbar

    call SolModel_TP_Update(TdgK,Pbar,Slv)

    vSolvDat(i)= Slv%Dat

  end do
  !-----------/
  !
  !---------compute the species' thermodyn' properties- build table tGrt
  call DtbSpc_Table_Build(vTPpath,vSpcDtb,vSpc,            tGrt,tVol)
  !--------/
  !
  !--- write tVol or tGrt to a file --
  if(trim(sCode)=="VOLUME") then
    call DtbSpc_Table_Write( &
    & trim(sCode),vTPpath,vSolvDat,vSpcDtb,vSpc,tVol)
  else
    call DtbSpc_Table_Write( &
    & trim(sCode),vTPpath,vSolvDat,vSpcDtb,vSpc,tGrt)
  end if
  !
  deallocate(tGrt)
  deallocate(tVol)
  deallocate(vSolvDat)
  deallocate(vTPpath)
  !
end subroutine Dtb_Tabulate

!-----------------------------------------------------------------------
!-- setting logK=0 for primary species
!-- write a logK database for current run,
!-----------------------------------------------------------------------
subroutine Dtb_Tabulate_ForSystem( &
& vTPpath,   &
& vCpn,vSpc,vSpcDtb, &
& vLCi,vLAx, &
& vOrdPr,    &
& tNuSp)
  use M_IoTools,    only: GetUnit, OutStrVec
  use M_Files,      only: DirOut,Files_Index_Write
  use M_Dtb_Const,  only: T_CK
  use M_T_Tpcond,   only: T_TPCond
  use M_T_SolModel, only: T_SolModel,T_SolModelDat
  use M_SolModel_Tools,only: SolModel_TP_Update
  use M_Numeric_Const,only: Ln10
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_Species,T_SpeciesDtb
  !
  use M_Dtb_Calc,   only: DtbSpc_Table_Build
  !
  !---------------------------------NEW 2011-03, for discretized species
  use M_Global_Vars,only: vMixModel,vDiscretModel,vDiscretParam
  use M_T_MixModel, only: MixModel_Param_Update
  use M_DiscretModel_Tools,only: DiscretSpecies_TP_Update
  !--/
  !
  type(T_TPCond),    intent(in):: vTPpath(:)
  type(T_Component), intent(in):: vCpn(:)
  type(T_Species),   intent(in):: vSpc(:)
  type(T_SpeciesDtb),intent(in):: vSpcDtb(:)
  logical,           intent(in):: vLCi(:),vLAx(:)
  integer,           intent(in):: vOrdPr(:)
  real(dp),          intent(in):: tNuSp(:,:)
  !
  integer :: iTP,iSp,I,J,N,nCp,F,iMM
  real(dp):: X,Tk,Pb !, T, P
  character(len=15):: Cc
  type(T_SolModel):: Slv
  type(T_Species):: S
  type(T_SolModelDat),allocatable:: vSolvDat(:)
  type(T_Species),    allocatable:: vSpcTmp(:)
  real(dp),           allocatable:: tGrt(:,:),tVol(:,:)
  !
  N= size(vTPpath)
  !
  allocate(vSolvDat(1:N))
  allocate(tGrt(1:size(vSpc),1:N))
  tGrt= Zero
  allocate(tVol(1:size(vSpc),1:N))
  tVol= Zero
  !
  do iTP=1,N
    !! call SolventDat_Default(vSolvDat(iTP))
    Tk= vTPpath(iTP)%TdgC +T_CK
    Pb= vTPpath(iTP)%Pbar
    call SolModel_TP_Update(Tk,Pb,Slv)
    vSolvDat(iTP)= Slv%Dat
  end do
  !
  call DtbSpc_Table_Build(vTPpath,vSpcDtb,vSpc,               tGrt,tVol)

  !---NEW 2011-03
  !-------------------------------update species built by discretization
  if(size(vDiscretModel)>0) then

    allocate(vSpcTmp(size(vSpc)))  ;  vSpcTmp(:)= vSpc(:)

    do iTP=1,N

      Tk= vTPpath(iTP)%TdgC +T_CK
      Pb= vTPpath(iTP)%Pbar

      !-----------------update T,P dependent parameters in mixing models
      do iMM=1,size(vMixModel)
        call MixModel_Param_Update( &
        & Tk,Pb,       &  !in
        & vMixModel(iMM)) !out
      end do
      !--------------/

      do iSp= 1,size(vSpcTmp)
        if(vSpc(iSp)%iDtb>0) &
        & vSpcTmp(iSp)%G0rt= tGrt(iSp,iTP)
      end do

      call DiscretSpecies_TP_Update( &
      & vMixModel,     & !IN
      & Tk,Pb,         & !IN
      & vDiscretModel, & !IN
      & vDiscretParam, & !IN
      & vSpcTmp)         !INOUT

      do iSp= 1,size(vSpcTmp)
        if(vSpc(iSp)%iDiscret/=0) &
        & tGrt(iSp,iTP)= vSpcTmp(iSp)%G0rt
      end do

    end do

    deallocate(vSpcTmp)

  end if
  !------------------------------/update species built by discretization
  !
  nCp=size(vCpn)
  !
  !--------------------------------logK from current prim'species set--/
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_logk_prim.dtb")
  !
  !--------------------------------------------------------index in html
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_logk_prim.dtb",&
  & "table of logK from prim'species, for current run")
  !----------------------------------------------------/
  !
  !-------------------------------------------------------TP.TABLE block
  write(F,'(/,A)') "TP.TABLE"
  call OutStrVec(F,vTPpath(1:N)%TdgC,Opt_S="TdgC")
  call OutStrVec(F,vTPpath(1:N)%Pbar,Opt_S="Pbar")
  write(F,'(A,/)') "END TP.TABLE"
  !---------------------------------------------------/
  !
  !--------------------------------------------------------SOLVENT block
  write(F,'(/,A)') "SOLVENT"
  write(F,'(A)') "NAME H2O"
  call OutStrVec(F,vSolvDat(1:N)%Rho, Opt_S="Rho")
  call OutStrVec(F,vSolvDat(1:N)%Eps, Opt_S="Eps")
  call OutStrVec(F,vSolvDat(1:N)%dhA, Opt_S="DHA")
  call OutStrVec(F,vSolvDat(1:N)%dhB, Opt_S="DHB")
  call OutStrVec(F,vSolvDat(1:N)%Bdot,Opt_S="BDOT")
  write(F,'(A,/)') "END SOLVENT"
  !----------------------------------------------------/
  !
  !------------------------------------------------------- SPECIES block
  write(F,'(/,A)') "SPECIES" !
  !-------------------------------------first, write the primary species
  DoPrim: do I=1,size(vSpc)
    !
    if(.not. (vLCi(I) .or. vLAx(I))) cycle DoPrim
    !
    S=vSpc(I)
    if(S%iDtb>0) then
      Cc= trim(vSpcDtb(S%iDtb)%DtbTrace)
    else
      Cc= "_"
    end if
    !
    !MODif!2016/02!
    !write(F, '(A3,A1, A15,A1, A23,A1, A63,A1,$)') &
    write(F, '(4(A,A1),$)') &
    & trim(S%Typ),    T_, &
    & trim(Cc),       T_, &
    & trim(S%NamSp),  T_, &
    & trim(S%Formula),T_
    
    !write(F,'(4(A,A1),$)') &
    !& trim(S%Typ),    T_, &
    !& trim(Cc),       T_, &
    !& trim(S%NamSp),   T_, &
    !& trim(S%Formula),T_
    
    if(S%Typ=="AQU") then
      !for aqu'species, write size parameter
      !print *,"S%NamSp=",trim(S%NamSp)
      write(F, '(F15.1,A1,$)') S%AquSize, T_
    else
      !for non aqueous species, write density
      write(F, '(G15.8,A1,$)') S%WeitKg / S%V0, T_
    end if
    
    !MODif!2016/02!
    !write(F,'(*(A,A1))') ("0.000", T_,J=1,N)
    write(F,'(*(A,A1))') ("0.000", T_,J=1,N)
    
  end do DoPrim
  !------------------------------------/first, write the primary species
  !
  !------------------------------------then, write the secondary species
  DoSecond: do I=1,size(vSpc)
    !
    if(vLCi(I) .or. vLAx(I)) cycle DoSecond
    !
    S=vSpc(I)
    if(S%iDtb>0) then  ;  Cc= trim(vSpcDtb(S%iDtb)%DtbTrace)
    else               ;  Cc= "_"
    end if
    !
    !MODif!2016/02!
    !write(F, '(A3,A1, A15,A1, A23,A1, A63,A1,$)') &
    write(F, '(4(A,A1),$)') &
    & trim(S%Typ),     T_, &
    & trim(Cc),        T_, &
    & trim(S%NamSp),   T_, &
    & trim(S%Formula), T_
    !
    if(S%Typ=="AQU") then
      !-- for aqu'species, write size parameter
      write(F, '(F15.1,A1,$)') S%AquSize, T_
    else
      !-- for non aqueous species, write density
      if(S%V0<1.D-9) then
        write(F,'(A,A1,$)') "_?_",T_
      else
        write(F, '(G15.8,A1,$)') S%WeitKg / S%V0, T_
      end if
    end if
    !
    do J=1,N
      ! X=  vSpc(vOrdAs(I))%G0rt &
      ! & - dot_product(tNuAs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
      X=  tGrt(I,J) &
      & - dot_product(tNuSp(I,1:nCp),tGrt(vOrdPr(1:nCp),J))
      write(F,'(G15.8,A1,$)') -X /Ln10, T_
    end do
    !
    write(F,*)
    !
  end do DoSecond
  !-----------------------------------/then, write the secondary species
  !
  write(F,'(A,/)') "END SPECIES"
  !-------------------------------------------------------/SPECIES block
  !
  close(F)
  !
  if(idebug>2) print '(A)',"Table of LogK in "//trim(DirOut)//"_logk.dtb"
  !------------------------------/ logK from current prim'species set --
  !
  !---------------------------------------------------------------------
  !-------------------------------------------------- logK from database
  !-------------------------without refering to current prim'species set
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_logk_elem.dtb")
  !
  !-------------------------------------------------------index in html
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_logk_elem.dtb",&
  & "table of logK from elements, for current run")
  !----------------------------------------------------/
  !
  !-------------------------------------------------------TP.TABLE block
  write(F,'(/,A)') "TP.TABLE"
  call OutStrVec(F,vTPpath(1:N)%TdgC,Opt_S="TdgC")
  call OutStrVec(F,vTPpath(1:N)%Pbar,Opt_S="Pbar")
  write(F,'(A,/)') "END TP.TABLE"
  !------------------------------------------------------/
  !
  !--------------------------------------------------------SOLVENT block
  write(F,'(/,A)') "SOLVENT"
  write(F,'(A)') "NAME H2O"
  call OutStrVec(F,vSolvDat(1:N)%Rho, Opt_S="Rho")
  call OutStrVec(F,vSolvDat(1:N)%Eps, Opt_S="Eps")
  call OutStrVec(F,vSolvDat(1:N)%dhA, Opt_S="DHA")
  call OutStrVec(F,vSolvDat(1:N)%dhB, Opt_S="DHB")
  call OutStrVec(F,vSolvDat(1:N)%Bdot,Opt_S="BDOT")
  write(F,'(A,/)') "END SOLVENT"
  !-------------------------------------------------------/SOLVENT block
  !
  !--------------------------------------------------------SPECIES block
  write(F,'(/,A)') "SPECIES" !
  !-------------------------------------first, write the primary species
  DoPrim2: do I=1,size(vSpc)
    !
    if(.not. (vLCi(I) .or. vLAx(I))) cycle DoPrim2
    !
    S=vSpc(I)
    if(S%iDtb>0) then
      Cc= trim(vSpcDtb(S%iDtb)%DtbTrace)
    else
      Cc= "_"
    end if
    !
    write(F, '(4(A,A1),$)') &
    & trim(S%Typ),      T_, &
    & trim(Cc),         T_, &
    & trim(S%NamSp),    T_, &
    & trim(S%Formula),  T_
    !write(F,'(4(A,A1),$)') &
    !& trim(S%Typ),     T_, &
    !& trim(Cc),        T_, &
    !& trim(S%NamSp),   T_, &
    !& trim(S%Formula), T_
    if(S%Typ=="AQU") then
      !for aqu'species, write size parameter
      write(F, '(F15.1,A1,$)') S%AquSize, T_
    else
      !for non aqueous species, write density, kg/m3
      write(F, '(G15.8,A1,$)') S%WeitKg / S%V0, T_
    end if
    do J=1,N
      X=  tGrt(I,J)
      write(F,'(G15.8,A1,$)') -X /Ln10, T_
    end do
    !
    write(F,*)
  end do DoPrim2
  !------------------------------------/first, write the primary species
  !
  !------------------------------------then, write the secondary species
  DoSecond2: do I=1,size(vSpc)
    !
    if(vLCi(I) .or. vLAx(I)) cycle DoSecond2
    !
    S=vSpc(I)
    if(S%iDtb>0) then
      Cc= trim(vSpcDtb(S%iDtb)%DtbTrace)
    else
      Cc= "_"
    end if
    !
    write(F, '(4(A,A1),$)') &
    & trim(S%Typ),    T_, &
    & trim(Cc),       T_, &
    & trim(S%NamSp),  T_, &
    & trim(S%Formula),T_
    !
    if(S%Typ=="AQU") then
      !-- for aqu'species, write size parameter
      write(F, '(F15.1,A1,$)') S%AquSize, T_
    else
      !-- for non aqueous species, write density
      if(S%V0<1.D-9) then
        write(F,'(A,A1,$)') "_?_",T_
      else
        write(F, '(G15.8,A1,$)') S%WeitKg / S%V0, T_
      end if
    end if
    !
    do J=1,N
      X=  tGrt(I,J)
      write(F,'(G15.8,A1,$)') -X /Ln10, T_
    end do
    !
    write(F,*)
    !
  end do DoSecond2
  !-----------------------------------/then, write the secondary species
  !
  write(F,'(A,/)') "END SPECIES"
  !-------------------------------------------------------/SPECIES block
  close(F)
  !
  if(idebug>2) print '(A)',"Table of LogK in "//trim(DirOut)//"_logk_elem.dtb"
  !--------------------------------------------------/logK from database
  !
  deallocate(tGrt)
  deallocate(vSolvDat)
  !
end subroutine Dtb_Tabulate_ForSystem

subroutine Dtb_Tabulate_EntropyZero(vS0Ele)
  use M_T_Species,  only: Species_EntropyZero
  use M_Global_Vars,only: vEle,vSpc
  !
  real(dp),intent(out):: vS0Ele(:)
  !
  integer :: I
  !
  do I=1,size(vSpc)
    vS0Ele(I)= Species_EntropyZero(vEle,vSpc(I))
    if(iDebug>2) write(71,'(A,A1,G15.6)') vSpc(I)%NamSp,T_,vS0Ele(I)
  end do
  !
  return
end subroutine Dtb_Tabulate_EntropyZero

subroutine DtbSpc_Table_Write(sCode,vTPCond,vSolModlDat,vSpcDtb,vSpc,tGrt)
  use M_Numeric_Const,only: Ln10
  use M_Dtb_Const, only: T_CK,R_jK,Tref
  use M_IoTools,   only: GetUnit, OutStrVec
  use M_Files,     only: DirOut,Files_Index_Write
  !
  use M_T_Species, only: T_Species, T_SpeciesDtb
  use M_T_Tpcond,  only: T_TPCond
  use M_T_SolModel,only: T_SolModelDat
  !---------------------------------------------------------------------
  character(len=*),   intent(in):: sCode
  type(T_SolModelDat),intent(in):: vSolModlDat(:)
  type(T_TPCond),     intent(in):: vTPCond(:)
  type(T_SpeciesDtb), intent(in):: vSpcDtb(:)
  type(T_Species),    intent(in):: vSpc(:)
  real(dp),           intent(in):: tGrt(:,:)
  !---------------------------------------------------------------------
  type(T_Species):: S
  integer :: f,iSp,jTP,N,M
  real(dp):: X
  character(len=15):: Cc
  real(dp),allocatable:: vS0Ele(:)
  !---------------------------------------------------------------------
  N= size(vTPCond)
  M= size(vSpc)

  call GetUnit(F)
  open(F,file=trim(DirOut)//"_"//trim(sCode)//".dtb")
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_"//trim(sCode)//".dtb",&
  & trim(sCode)//" table for current run")

  if(trim(sCode)=="GIBBS2") then
    allocate(vS0Ele(1:M))
    call Dtb_Tabulate_EntropyZero(vS0Ele)
  end if
  !
  !---------------------------------------------------------write header
  select case(trim(sCode))
  !
  case("LOGK","LOGK_EQ36")
    !
    write(F,'(/,A)') "TP.TABLE"
    !
    call OutStrVec(F,vTPCond(:)%TdgC,Opt_S="TdgC")
    call OutStrVec(F,vTPCond(:)%Pbar,Opt_S="Pbar")
    !
    write(F,'(A,/)') "END TP.TABLE"
    !
    write(F,'(/,A)') "SOLVENT"
    !
    write(F,'(A)') "NAME H2O"
    call OutStrVec(F,vSolModlDat(:)%Rho, Opt_S="Rho")
    call OutStrVec(F,vSolModlDat(:)%Eps, Opt_S="Eps")
    call OutStrVec(F,vSolModlDat(:)%dhA, Opt_S="DHA")
    call OutStrVec(F,vSolModlDat(:)%dhB, Opt_S="DHB")
    call OutStrVec(F,vSolModlDat(:)%Bdot,Opt_S="BDOT")
    !
    write(F,'(A,/)') "END SOLVENT"
    !
    write(F,'(/,A)') "SPECIES" !
    !
  !
  case("GIBBS","GIBBS2","VOLUME")
    !
    select case(trim(sCode))
    case("GIBBS")   ;  write(F,'(A)') "!Gibbs_FreeEnergy_alaBenson"
    case("GIBBS2")  ;  write(F,'(A)') "!Gibbs_FreeEnergy_alaBerman"
    case("VOLUME")  ;  write(F,'(A)') "!Molar Volume, cm3"
    end select

    !-- write array of Temperatures --
    write(F,'(4(A,A1),$)') &
    & ".Type",T_,".Trace",T_,".Species",T_,".ECFORM",T_
    ! call OutStrVec(F,vTPCond(:)%TdgC,Opt_S="TdgC")
    call OutStrVec(F,vTPCond(:)%TdgC)
    
    !-- write array of Pressures --
    write(F,'(4(A,A1),$)') &
    & ".Type",T_,".Trace",T_,".Species",T_,".ECFORM",T_
    ! call OutStrVec(F,vTPCond(:)%Pbar,Opt_S="Pbar")
    call OutStrVec(F,vTPCond(:)%Pbar)
  !
  end select
  !--------------------------------------------------------/write header
  !
  !--------------------------------------------------------write results
  do iSp=1,M
    !
    S=vSpc(iSp)
    !
    if(S%iDtb>0) then  ;  Cc= trim(vSpcDtb(S%iDtb)%DtbTrace)
    else               ;  Cc= "_"
    end if
    !
    write(F,'(4(A,A1),$)') &
    & trim(S%Typ),     T_, &
    & trim(Cc),        T_, &
    & trim(S%NamSp),   T_, &
    & trim(S%Formula), T_
    !
    select case(sCode)
    !
    case("LOGK")
      if(S%Typ=="AQU") then
        !-- for aqu'species, write size parameter --
        write(F, '(F15.1,A1,$)') S%AquSize, T_
      else
        !-- for non aqueous species, write density --
        write(F, '(G15.8,A1,$)') S%WeitKg /S%V0, T_
      end if
      !-- write array of log10(-G/RT), i.e. "logK" --
      do jTP=1,N
        write(F,'(G15.8,A1,$)') -tGrt(iSp,jTP) /Ln10, T_
      end do
    !
    case("GIBBS")
      !--- write Gibbs free energy --
      do jTP= 1,N
        X= tGrt(iSp,jTP) *(vTPCond(jTP)%TdgC +T_CK) *R_jK
        ! X= X - Tref *vS0Ele(iSp) !-> to Berman-Brown
        write(F,'(G15.8,A1,$)') X, T_
      end do
    !
    case("GIBBS2")
      !--- write Gibbs free energy, in Berman-Brown convention --
      do jTP= 1,N
        X= tGrt(iSp,jTP) *(vTPCond(jTP)%TdgC +T_CK) *R_jK
        X= X - Tref *vS0Ele(iSp) !-> to Berman-Brown
        write(F,'(G15.8,A1,$)') X, T_
      end do
    !
    case("VOLUME")
      !--- write molar volume --
      do jTP= 1,N
        X= tGrt(iSp,jTP)*1.0D6
        write(F,'(G15.8,A1,$)') X, T_
      end do
    !
    end select
    !
    write(F,*)
    !
  end do
  !-------------------------------------------------------/write results
  !
  if (trim(sCode)=="LOGK") write(F,'(A,/)') "END SPECIES"
  !
  if(trim(sCode)=="GIBBS2") deallocate(vS0Ele)
  !
  close(F)
  !
  if(idebug>1) then
    print '(A)',"Table of "//trim(sCode)//" in "//trim(DirOut)//"_"//trim(sCode)//".dtb"
    if(iDebug>2) call Pause_
  end if
  !
end subroutine DtbSpc_Table_Write

subroutine Dtb_Test_EQ36
!--
!-- tabulate logK of formation of all species of vSpc
!-- along the EQ36 series of T,P points
!--
  use M_Dtb_Const, only: T_CK
  use M_Dtb_Calc,  only: DtbSpc_Table_Build
  use M_T_SolModel,only: T_SolModel,T_SolModelDat !,SolventDat_Default
  use M_SolModel_Tools,only: SolModel_TP_Update
  !
  use M_Global_Vars,only: vSpcDtb,vSpc
  !
  type(T_TPCond)      :: vTPcond(8)
  type(T_SolModel)     :: Slv
  type(T_SolModelDat)  :: vSolvDat(8)
  real(dp),allocatable:: tGrt(:,:),tVol(:,:)
  integer:: i
  !
  allocate(tGrt(1:size(vSpc),1:8))
  tGrt= Zero
  allocate(tVol(1:size(vSpc),1:8))
  tVol= Zero
  !
  vTPcond(1)%TdgC= 0.01D0  ;  vTPcond(1)%Pbar= 1.013D0
  vTPcond(2)%TdgC= 25.0D0  ;  vTPcond(2)%Pbar= 1.013D0
  vTPcond(3)%TdgC= 60.0D0  ;  vTPcond(3)%Pbar= 1.013D0
  vTPcond(4)%TdgC= 100.D0  ;  vTPcond(4)%Pbar= 1.013D0
  vTPcond(5)%TdgC= 150.D0  ;  vTPcond(5)%Pbar= 4.757D0
  vTPcond(6)%TdgC= 200.D0  ;  vTPcond(6)%Pbar= 15.34D0
  vTPcond(7)%TdgC= 250.D0  ;  vTPcond(7)%Pbar= 39.74D0
  vTPcond(8)%TdgC= 300.D0  ;  vTPcond(8)%Pbar= 85.84D0
  !
  do i=1,8
    !! call SolventDat_Default(vSolvDat(i))
    !! provisional !! should call calculation of vSolvDat !!
    call SolModel_TP_Update(vTPcond(i)%TdgC +T_CK,vTPcond(i)%Pbar,Slv)
    vSolvDat(i)= Slv%Dat
  end do
  !
  call DtbSpc_Table_Build(vTPCond,vSpcDtb,vSpc,tGrt,tVol)
  !
  call DtbSpc_Table_Write("LOGK",vTPCond,vSolvDat,vSpcDtb,vSpc,tGrt)
  !
  deallocate(tGrt)
  !
end subroutine Dtb_Test_EQ36

subroutine Dtb_Tabulate_PSatH2O
!--
!-- tabulate (T,P) series a la EQ36, up to critical temp' of water
!--
  use M_Dtb_Const,only: T_CK
  use M_IOTools,  only: GetUnit
  use M_Fluid_Calc
  !
  real(dp):: TdgC,TdgK,Pbar,Psat
  integer :: f
  !
  call GetUnit(f)
  open(f,file="psath2o.restab")
  !call Files_Index_Write(fHtm,&
  !& "psath2o.restab",&
  !& "T(deg)C-P(bar) along gas saturation")

  TdgC=Zero
  do

    TdgK= TdgC +T_CK

    if (TdgK<=647.25D0) call Eos_H2O_Psat(TdgK,Psat)

    if(TdgC<=100._dp) then  ;  Pbar= 1.013D0
    else                    ;  Pbar= Psat
    end if

    write(f,'(F7.2,A1,F7.2)') TdgC,t_,Pbar

    TdgC= TdgC +5.D0
    if(TdgK>647.25D0) exit

  end do

  close(f)
end subroutine Dtb_Tabulate_PSatH2O

subroutine Dtb_Test_Fluid
  use M_Numeric_Const,only: ln10
  use M_Dtb_Const,  only: T_CK, R_JK,Tref,Pref
  use M_Files,      only: DirOut
  use M_IOTools,    only: GetUnit
  use M_TPcond_Read,only: TPpath_Read,TPgrid_Build
  use M_Path_Vars,  only: vTPpath
  !
  use M_Fluid_Calc
  ! use M_FluidMix_KerJac
  !
  logical :: Ok
  real(dp):: TdgC,TdgK,Pbar,Grt,H,S,V_m3
  real(dp):: LnFug
  integer :: i,f1,f2,f3
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! real! (dp):: XCO2
  ! real(dp):: FugC_H2O,FugC_CO2,ActH2O,ActCO2
  ! 
  ! TdgC= 400.0D0
  ! Pbar= 10.0D3
  ! TdgK= TdgC + T_CK
  ! 
  ! XCO2= 0.01D0
  ! do
  ! 
  !   call EoS_KerJac_CalcMix( &
  !   & TdgK,Pbar,XCO2, & !in
  !   & FugC_H2O,FugC_CO2,ActH2O,ActCO2) !out
  ! 
  !   write(11,'(5(G15.6,1X))') XCO2,FugC_H2O,FugC_CO2,ActH2O,ActCO2
  !   
  !   XCO2= XCO2 + 0.01D0 
  ! 
  !  if(XCO2>0.995D0) exit
  ! 
  ! end do
  !
  ! return
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(idebug>1) &
  & print *,"Compute the thermodynamic properties of H2O"

  TdgK= Tref
  Pbar= Pref
  !
  call TPgrid_Build(Ok)
  if(.not. Ok) call TPpath_Read(TdgK,Pbar)
  !
  call GetUnit(f1)
  open(f1,file=trim(DirOut)//"_dtbfluid_ghio.restab")
  write(f1,'(10(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_, &
  & "G",    t_,"H",   t_,"S",        t_,"H-TS",t_,"Rho", t_
  !
  call GetUnit(f2)
  open(f2,file=trim(DirOut)//"_dtbfluid_haar.restab")
  write(f2,'(7(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_,"G",t_,"Rho",t_
  !
  call GetUnit(f3)
  open(f3,file=trim(DirOut)//"_dtbfluid_supcrt.restab")
  write(f3,'(10(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_, &
  & "G",    t_,"H",   t_,"S",        t_,"H-TS",t_,"Rho", t_
  !
  do i=1,size(vTPpath)

    TdgC= vTPpath(i)%TdgC
    Pbar= vTPpath(i)%Pbar
    TdgK= TdgC +T_CK

    call EosFluid_Calc( &
    & "H2O","HAAR.GHIORSO",TdgK,Pbar, &
    & Grt,H,S,V_m3,LnFug,Ok)
    if(Ok) then
      write(f1,'(A,A1,3(G15.6,A1),$)') &
      & "HAARGHIO",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      write(f1,'(6(G15.6,A1))') &
      & -Grt/ln10,t_,  Grt*R_JK*TdgK, t_,  &
      &  H,       t_,  S,             t_,  &
      &  H-TdgK*S,t_,  0.01852D0/V_m3,t_
    end if

    call EosFluid_Calc( &
    & "H2O","HAAR",TdgK,Pbar, &
    & Grt,H,S,V_m3,LnFug,Ok)
    if(Ok) then
      write(f2,'(A,A1,3(G15.6,A1),$)') &
      & "HAAR",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      write(f2,'(3(G15.6,A1))') &
      &  -Grt/ln10,t_,  Grt*R_JK*TdgK,t_, 0.01852D0/V_m3,t_
    end if

    call EosFluid_Calc( &
    & "H2O","SUPCRT",TdgK,Pbar, &
    & Grt,H,S,V_m3,LnFug,Ok)
    if(Ok) then
      write(f3,'(A,A1,3(G15.6,A1),$)') &
      & "SUPCRT",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      write(f3,'(6(G15.6,A1))') &
      & -Grt/ln10,t_,  Grt*R_JK*TdgK,  t_,  &
      &  H,       t_,  S,              t_,  &
      &  H-TdgK*S,t_,  0.01852D0/V_m3, t_
    end if

  end do

  close(f1)
  close(f2)
  close(f3)

  if(idebug>1) then
    print *,">> results in files dtbfluid_xxx.restab"
    if(iDebug>2) call Pause_
  end if

  deallocate(vTPpath)
  !
end subroutine Dtb_Test_Fluid

subroutine Dtb_Test_Species
!--
!-- tabulate the properties of all species
!--
  use M_TPcond_Read,only: TPpath_Read
  use M_Dtb_Write
  !
  use M_Dtb_Vars
  use M_Path_Vars, only: vTPpath
  !
  real(dp):: TdgK,Pbar
  !---------------------------------------------------------------------
  TdgK= 298.15D0
  Pbar= 1.0D0
  !
  call TPpath_Read(TdgK,Pbar)
  !
  if(size(vDtbAquHkf)>0) call DtbAquHkf_Tabulate(vTPpath)
  ! if(size(vDtbAquThr)>0) call DtbAquThr_Tabulate(vTPpath)
  ! if(size(vDtbMinHkf)>0) call DtbMinHkf_Tabulate(vTPpath)
  if(size(vDtbMinThr)>0) call DtbMin_Tabulate("THR",vTPpath)
  if(size(vDtbMinHkf)>0) call DtbMin_Tabulate("HKF",vTPpath)
  !
  if(size(vDtbAquHkf)>0) then
    if(size(vDtbMinThr)>0) call DtbSys_Tabulate("THR",vTPpath)
    if(size(vDtbMinHkf)>0) call DtbSys_Tabulate("HKF",vTPpath)
  end if
  !
  return
end subroutine Dtb_Test_Species

subroutine Dtb_Test_H2OHkf
!--
!-- tabulate the properties of H2O used by the HKF model
!--
  use M_Dtb_Const,  only: T_CK
  use M_IOTools,    only: GetUnit
  use M_TPcond_Read,only: TPpath_Read,TPgrid_Build
  use M_Path_Vars,  only: vTPpath
  use M_Fluid_Calc
  use M_T_DtbH2OHkf
  !---------------------------------------------------------------------
  logical :: Ok
  real(dp):: TdgC,TdgK,Pbar,Psat
  integer :: f1,f2
  integer :: i
  type(T_DtbH2OHkf):: pW
  !---------------------------------------------------------------------
  TdgK= 298.15D0
  Pbar= 1.01325D0
  !
  call TPgrid_Build(Ok)
  if(.not. Ok) call TPpath_Read(TdgK,Pbar)
  !
  call GetUnit(f1)
  open(f1,file="dtbh2o.restab")
  write(f1,'(11(A,A1))') &
  & "TdgC",t_,"log(Pbar)",t_,"Pbar",t_,&
  & "Rho",t_,"Eps",t_,"Q",t_,"X",t_,"Y",t_,"Z",t_,"dhA",t_,"dhB",t_
  !
  call GetUnit(f2)
  open(f2,file="dtbh2o_alfa.restab")
  write(f2,'(9(A,A1))') &
  & "TdgC",t_,"log(Pbar)",t_,"Pbar",t_,&
  & "Rho",t_,"Alfa",t_,"dAlfdT",t_,"Beta",t_,"dBetdT",t_,"dBetdP",t_
  !
  do i=1,size(vTPpath)
    !
    TdgC= vTPpath(i)%TdgC
    TdgK= TdgC +T_CK
    Pbar= vTPpath(i)%Pbar
    !
    Psat= Zero
    if (TdgK<=647.25D0) call Eos_H2O_Psat(TdgK,Psat)
    !
    if(TdgK>647.25D0 .or. Pbar>Psat) then
    !-> restrict to liquid or supercritic'domain
      !
      !!call CalcRhoH2O( &
      !!& TdgK,Pbar, & !IN
      !!& RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP) !OUT
      !!write(f1,'(8(G15.6,A1))') &
      !!& TdgC,t_,log10(Pbar),t_, &
      !!& RH2O,t_,Alfa,t_,Beta,t_,dAlfdT,t_,dBetdT,t_,dBetdP,t_
      !
      call DtbH2OHkf_Calc(TdgK,Pbar,pW)
      !
      write(f1,'(11(G15.6,A1))') &
      & TdgC,t_,        log10(Pbar),t_, Pbar,t_, &
      & pW%Rho*1.D3,t_, pW%Eps,t_,      pW%Q,t_, &
      & pW%X,t_,        pW%Y,t_,        pW%Z,t_, &
      & pW%dhA,t_,      pW%dhB,t_
      !
      write(f2,'(9(G15.6,A1))') &
      & TdgC,t_,        log10(Pbar),t_, Pbar,t_, &
      & pW%Rho*1.D3,t_, &
      & pW%Alfa,t_,     pW%dAlfdT,t_, &             !thermal expansion.
      & pW%Beta,t_,     pW%dBetdT,t_, pW%dBetdP,t_  !compressibility
      !
    end if
    !
  end do
  close(f1)
  close(f2)
  !
  deallocate(vTPpath)
  !
  return
end subroutine Dtb_Test_H2OHkf

subroutine Dtb_Test_AquHkf
!--
!-- test computations on aqueous species (HKF model, HSV data)
!-- tabulate values of logK and partial molar volume
!-- if vSpc contains also minerals,
!-- their logK and density (kg/m3) are tabulated
!-- along points of a TP.GRID
!--
  use M_Dtb_Const,    only: T_CK,Tref,Pref
  use M_Files,        only: DirOut
  use M_IOTools,      only: GetUnit
  use M_Numeric_Const,only: Ln10
  use M_T_Species,    only: T_Species
  use M_Dtb_Calc,     only: Species_TP_Update_fromDtb
  use M_T_DtbH2OHkf,  only: T_DtbH2OHkf,DtbH2OHkf_Calc
  use M_Fluid_Calc,   only: Eos_H2O_Psat
  use M_TPcond_Read,  only: TPgrid_Build,TPpath_Read
  !
  use M_Global_Vars,  only: vSpcDtb,vSpc
  use M_Dtb_Vars,     only: vDtbAquHkf
  use M_Path_Vars,    only: vTPpath
  !
  real(dp):: TdgC,TdgK,Pbar,X
  !
  logical:: Ok
  integer:: N, I, J
  integer:: f1,f2
  type(T_DtbH2OHkf):: PropsH2O
  !
  call TPgrid_Build(Ok)
  if(.not. Ok) call TPpath_Read(Tref,Pref)
  !
  N= size(vSpc)
  !
  if(N<=0) return
  !
  if(size(vDtbAquHkf)<1) then
    print '(A)',"no HSV data on species"
    return
  end if
  !
  call GetUnit(f1)
  open(f1,file=trim(DirOut)//"_tpgrid_logk.restab")
  write(f1,'(3(A,A1),$)') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !
  call GetUnit(f2)
  open(f2,file=trim(DirOut)//"_tpgrid_dens.restab")
  write(f2,'(3(A,A1),$)') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !
  do I=1,size(vSpc) !N
    !!if(vSpc(I)%Typ=="AQU") write(f1,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
    write(f1,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
  end do
  write(f1,*)
  do I=1,size(vSpc) !N
    !!if(vSpc(I)%Typ=="AQU") write(f1,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
    write(f2,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
  end do
  write(f2,*)
  !
  do N=1,size(vTPpath)
    !
    TdgC= vTPpath(N)%TdgC
    TdgK= TdgC +T_CK
    Pbar= vTPpath(N)%Pbar
    !
    !! print '(A,F7.2,1X,F7.2)',"Dtb_Test_AquHkf ",TdgC,Pbar
    write(f1,'(3(G15.6,A1),$)') TdgC,t_,LOG10(Pbar),t_,Pbar,t_
    write(f2,'(3(G15.6,A1),$)') TdgC,t_,LOG10(Pbar),t_,Pbar,t_
    !
    call DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
    !-> solvent properties, for aqu'species
    !
    do J=1,size(vSpc)
      !
      call Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,vSpc(J))
      !
      write(f1,'(G15.6,A1,$)') -vSpc(J)%G0rt/Ln10,t_
      !
      if(vSpc(J)%Typ=="AQU") then
        write(f2,'(G15.6,A1,$)')  vSpc(J)%V0,t_
      else
        X= vSpc(J)%WeitKg
        write(f2,'(G15.6,A1,$)')  X/vSpc(J)%V0,t_
      end if
      !
    end do
    !
    write(f1,*)
    write(f2,*)
    !
  end do
  !
  close(f1)
  close(f2)
  !
  deallocate(vTPpath)
  !
  return
end subroutine Dtb_Test_AquHkf

end module M_Dtb_Test

! subroutine Dtb_Test_AquHkf_
!   use M_Dtb_Const,  only: T_CK
!   use M_IOTools,    only: GetUnit
!   use M_Files,      only: DirOut
!   use M_Numeric_Const,only: Ln10
!   use M_T_Species,  only: T_Species
!   use M_Dtb_Vars,   only: vDtbAquHkf
!   use M_Global_Vars,only: vSpc
!   use M_Dtb_Calc,   only: DtbSpc_TP_Update
!   use M_Dtb_Read,   only: Dtb_Read_TPgrid
!   use M_T_DtbH2OHkf,only: T_DtbH2OHkf,DtbH2OHkf_Calc
!   !use M_T_DtbAquHkf,only: DtbAquHkf_Calc,DtbAquHkf_CalcThr
!   use M_Fluid_Calc ,only: Eos_H2O_Psat
!   !
!   !real(dp),intent(in):: T_Min, T_Max, T_delta
!   !real(dp),intent(in):: P_Min, P_Max, P_delta, P_ratio
!   !
!   logical :: Ok
!   real(dp):: T_Min, T_Max, T_delta, T_ratio
!   real(dp):: P_Min, P_Max, P_delta, P_ratio
!   real(dp):: TdgC,TdgK,Pbar,Psat,X
!   !
!   integer:: N, I, J
!   integer:: f1,f2
! 
!   type(T_DtbH2OHkf):: PropsH2O
!   !
!   T_Min=   0._dp     !5._dp     !
!   T_Max=   200._dp   !500._dp   !
!   T_delta= 5._dp     !5._dp     !
!   T_ratio= One
!   P_Min=   1._dp     !1._dp     !
!   P_Max=   1000._dp  !1000._dp  !
!   P_ratio= 1.2_dp    !1.1_dp    !
!   P_delta= 5._dp
!   !
!   call Dtb_Read_TPgrid( &
!   & Ok,T_Min,T_Max,T_ratio,T_delta,P_Min,P_Max,P_ratio,P_delta)
!   !
!   if(iDebug>3) then
!     print '(A,4G15.6)',"T_Min,T_Max,T_ratio,T_delta",T_Min,T_Max,T_ratio,T_delta
!     print '(A,4G15.6)',"P_Min,P_Max,P_ratio,P_delta",P_Min,P_Max,P_ratio,P_delta
!   call Pause_
!   end if
!   !
!   N= size(vSpc)
!   if(N<=0) return
!   if(size(vDtbAquHkf)<1) then
!     print '(A)',"no HSV data on species"
!     return
!   end if
!   !
!   call GetUnit(f1) ; open(f1,file=trim(DirOut)//"_tpgrid_logk.restab")
!   call GetUnit(f2) ; open(f2,file=trim(DirOut)//"_tpgrid_dens.restab")
!   write(f1,'(3(A,A1),$)') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
!   write(f2,'(3(A,A1),$)') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
!   do I=1,size(vSpc) !N
!     !!if(vSpc(I)%Typ=="AQU") write(f1,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
!     write(f1,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
!     write(f2,'(A,A1,$)') trim(vSpc(I)%NamSp),t_
!   end do
!   write(f1,*)
!   write(f2,*)
!   !
!   TdgC= T_Min
!   Do_T: do while(TdgC<T_Max)
!     Pbar= P_Min
!     Do_P: do while(Pbar<=P_Max)
!       TdgK= TdgC +T_CK
!       !print '(F7.2,1X,F7.2)',TdgC,Pbar
!       !write(f1,'(2(G15.6,A1),$)') TdgC,t_,Pbar,t_
!       !
!       Psat= Zero
!       if (TdgK<=647.25D0) call Eos_H2O_Psat(TdgK,Psat)
!       if(TdgK>647.25D0 .or. Pbar>Psat) then
!         write(f1,'(3(G15.6,A1),$)') TdgC,t_,log10(Pbar),t_,Pbar,t_
!         write(f2,'(3(G15.6,A1),$)') TdgC,t_,log10(Pbar),t_,Pbar,t_
!         !
!         call DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O) !-> solvent properties, for aqu'species
!         !
!         do J=1,size(vSpc)
!           call DtbSpc_TP_Update(TdgK,Pbar,PropsH2O,vSpc(J))
!           write(f1,'(G15.6,A1,$)') -vSpc(J)%G0rt/Ln10,t_
!           if(vSpc(J)%Typ=="GAS") then
!             X= exp(vSpc(J)%LnPhi -log(Pbar))
!             write(f1,'(G15.6,A1,$)') X,t_
!           end if
!           if(vSpc(J)%Typ=="AQU") then
!             write(f2,'(G15.6,A1,$)')  vSpc(J)%V0,t_
!           else
!             X= vSpc(J)%WeitKg
!             write(f2,'(G15.6,A1,$)')  X/vSpc(J)%V0,t_
!           end if
!         end do
!         write(f1,*)
!         write(f2,*)
!       end if
!       !write(f1,'(2(G15.6,A1))') Grt,t_,H,t_ !,S,t_,18.052D0/V_m3,t_
!       !
!       if(P_ratio/=One) then; Pbar= Pbar *P_ratio
!       else                 ; Pbar= Pbar + P_delta
!       end if
!       !
!     end do Do_P
!     TdgC= TdgC + T_delta
!   end do Do_T
!   close(f1)
!   close(f2)
!   !
!   return
! end subroutine Dtb_Test_AquHkf_

