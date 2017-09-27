module M_Dynam_Files
  use M_Kinds
  use M_Trace,only: fTrc,fHtm,T_,iDebug
  implicit none
  !
  private
  !
  public:: Dynam_Files_Init
  public:: Dynam_Files_Close
  public:: Dynam_Files_OpenLogs
  public:: Dynam_Files_WriteFasAff
  !
  integer,public:: &
  & fDynAct= 0, & !results: log10(activities) 
  & fDynMol= 0, & !mole numbers
  & fDynEle= 0, & !component mole numbers (different from element, in orthogonal option)
  & fDynElement= 0, & !element mole numbers
  & fDynQsK= 0, & !minerals
  & fDynMnK= 0, & !minerals
  & fDynGam= 0    !gammas

contains

subroutine Dynam_Files_Init
  use M_Files,      only: DirOut,Files_Index_Write
  use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
  use M_IoTools,    only: GetUnit
  use M_Equil_Write,only: Equil_Write_EnTete
  !
  use M_Global_Vars,only: vEle,vSpc
  use M_System_Vars,only: vCpn
  use M_Basis_Vars, only: vOrdAq,vPrmFw
  !
  use M_Dynam_Vars, only: TUnit,dTSav,fSavTime,fSavRate,DebNewt,DebJacob
  !
  integer:: I,nAq
  !
  nAq= count(vSpc%Typ(1:3)=="AQU")
  !
  if(fDynGam==0) then
    !
    call GetUnit(fDynGam)
    open(fDynGam,file=trim(DirOut)//"_gamma.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_gamma.restab", &
    & "DYNAMIC: activity coeff's of aqueous species")
    !
    write(fDynGam,'(2(A,A1))',advance="no") "STEP",T_,"TIME",T_
    call Equil_Write_EnTete(fDynGam,vSpc,vPrmFw,vSpc%Typ=="AQU") 
    !
  end if
  !
  if(fDynAct==0) then
    !
    call GetUnit(fDynAct)
    open(fDynAct,file=trim(DirOut)//"_activ.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_activ.restab", &
    & "DYNAMIC: at each time step, species activities")
    !
    write(fDynAct,'(2(A,A1))',advance="no") "STEP",T_,"TIME/"//TUnit,T_
    do I=1,nAq; write(fDynAct,'(A12,A1)',advance="no") vSpc(vOrdAq(I))%NamSp,T_; end do
    write(fDynAct,*)
    !
  end if
  !
  if(fDynMol==0) then
    !
    call GetUnit(fDynMol)
    open(fDynMol,file=trim(DirOut)//"_moleinbox.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_moleinbox.restab", &
    & "DYNAMIC: mole numbers of aqu.species in box")
    !
    write(fDynMol,'(2(A,A1))',advance="no") "STEP",T_,"TIME/"//TUnit,T_
    do I=1,nAq; write(fDynMol,'(A,A1)',advance="no") trim(vSpc(vOrdAq(I))%NamSp),T_; end do
    write(fDynMol,*)
    !
  end if
  !
  if(iDebug>2 .and. fDynElement==0) then
    !
    call GetUnit(fDynElement)
    open(fDynElement,file=trim(DirOut)//"_elements.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_elements.restab", &
    & "DYNAMIC: at each time step,"//  &
    & " mole numbers of elements within Box, in Fluid and in Minerals")
    !
    !-- element order will be the order in tStoikioAqu / tStoikioKin,
    !-- which are retrieved from vSpc(:)%vStoikio(:)
    !-- it is thus the element order in original vEle(:)
    !
    write(fDynElement,'(2(A,A1))',advance="no") "STEP",T_,"TIME/"//TUnit,T_
    do I=1,size(vEle); write(fDynElement,'(A,A1)',advance="no") &
    & vEle(I)%NamEl//"inFluid", T_; end do
    do I=1,size(vEle); write(fDynElement,'(A,A1)',advance="no") &
    & vEle(I)%NamEl//"inMin", T_; end do
    do I=1,size(vEle); write(fDynElement,'(A,A1)',advance="no") &
    & vEle(I)%NamEl//"Total", T_; end do
    write(fDynElement,*)
    !
  end if
  !
  if(fDynEle==0) then
    !
    call GetUnit(fDynEle)
    open(fDynEle,file=trim(DirOut)//"_elem.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_elem.restab", &
    & "DYNAMIC: at each time step,"//  &
    & " mole numbers of elements within Box, in Fluid and in Minerals")
    !
    !-- this header will be valid only when component--element !
    !-- -> to be MODIFIED in other cases,
    !-- e.g. when using orthogonal components !!
    !
    write(fDynEle,'(4(A,A1))',advance="no") &
    & "STEP",T_,"TIME/"//TUnit,T_,"DeltaDarcy",T_
    !! do I=1,size(vEle); write(fDynEle,'(A,A1)',advance="no") vEle(iCpnEle(I))%NamEl//"totF", T_; end do
    !! do I=1,size(vEle); write(fDynEle,'(A,A1)',advance="no") vEle(iCpnEle(I))%NamEl//"molal",T_; end do
    !! do I=1,size(vEle); write(fDynEle,'(A,A1)',advance="no") vEle(iCpnEle(I))%NamEl//"inMin",T_; end do
    write(fDynEle,'(100(A,A1))',advance="no") &
    & (trim(vCpn(I)%NamCp)//"totF", T_,I=1,size(vCpn))
    write(fDynEle,'(100(A,A1))',advance="no") &
    & (trim(vCpn(I)%NamCp)//"molal",T_,I=1,size(vCpn))
    write(fDynEle,'(100(A,A1))',advance="no") &
    & (trim(vCpn(I)%NamCp)//"inMin",T_,I=1,size(vCpn))
    write(fDynEle,*)
    !
  end if
  !
  if(fDynMnK==0) then
    !
    call GetUnit(fDynMnK)
    open(fDynMnK,file=trim(DirOut)//"_minmol.restab")
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_minmol.restab", &
    & "DYNAMIC: at each time step, pH, vol'fractions fluid and minerals, logQsK, etc")
    !
    call EnTeteFMnk(fDynMnK) !,bCell)
    !
  end if
  !
  if(dTSav>Zero) then
    call GetUnit(fSavTime)
    open(fSavTime,file=trim(DirOut)//"_time.restab")
    !
    call EnTeteFMnk(fSavTime) !,bCell)
    !
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_time.restab", &
    & "DYNAMIC: at regular time laps, pH, vol'fractions fluid and minerals, logQsK, etc")
  end if
  !
  if(iDebug>2) then
    if(fSavRate==0) then
      !
      call GetUnit(fSavRate)
      open(fSavRate,file=trim(DirOut)//"_rate.restab")
      !
      call EnTeteFRate(fSavRate)
      !
      call Files_Index_Write(fHtm, &
      & trim(DirOut)//"_rate.restab", &
      & "DYNAMIC: at each time step, details on mineral dissol/precip rates (surface,...)")
     end if
  end if
  !
  if(DebNewt) then
    !
    !fNewtF:  trace for Function values on aqu'species
    !fNewtR:  trace for Residual values on aqu'species
    !fNewtFm: trace for Function values on minerals
    !fNewtRm: trace for Residual values on minerals
    !
    call GetUnit(fNewtF);  open(fNewtF, file=trim(DirOut)//"_newtaqu1.log")
    call GetUnit(fNewtR);  open(fNewtR, file=trim(DirOut)//"_newtaqu2.log") !trace for Residual values
    !call GetUnit(fNewtFm); open(fNewtFm,file=trim(DirOut)//"_newtmin1.log")
    !call GetUnit(fNewtRm); open(fNewtRm,file=trim(DirOut)//"_newtmin2.log")
    !
    fNewt_I= 0
    !
    write(fNewtF,'(A,A1)',advance="no") "Index",T_
    call Equil_Write_EnTete(fNewtF,vSpc,vPrmFw,vSpc%Typ=="AQU",Str1="indx") 
    !
    write(fNewtR,'(A,A1)',advance="no") "Index",T_
    call Equil_Write_EnTete(fNewtR,vSpc,vPrmFw,vSpc%Typ=="AQU",Str1="indx") 
    !
    DebJacob= .true.
    !
  end if
  !
end subroutine Dynam_Files_Init

subroutine EnTeteFMnk(f) 
  use M_Global_Vars,only: vKinFas
  use M_Dynam_Vars, only: TUnit
  !
  integer,intent(in):: f
  !!logical,intent(in):: bCell
  !
  integer::i,N
  !
  N=size(vKinFas)
  !
  !write(f,'(4(A,A1))',advance="no") "!iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !do J=1,3
  !  do i=1,N; write(f,'(A,A1)', advance="no") trim(vKinFas(i)%Name),T_; end do
  !end do
  !write(f,*)
  !
  !write(f,'(4(A,A1))',advance="no") "!iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !do i=1,N; write(f,'(A,A1)',advance="no") "PhiM",T_;  end do
  !do i=1,N; write(f,'(A,A1)',advance="no") "LogQsK",T_;  end do
  !do i=1,N; write(f,'(A,A1)',advance="no") "MolNr",T_;  end do
  !write(f,*)
  !
  !!if(bCell) write(f,'(A,A1)',advance="no") "iCell",T_
  write(f,'(4(A,A1))',advance="no") "iStep", T_,"Time/"//TUNit,T_,"pH",T_,"PhiFluid",T_
  do i=1,N; write(f,'(A,A1)',advance="no") "PhiM_"//trim(vKinFas(i)%NamKF),T_;   end do
  do i=1,N; write(f,'(A,A1)',advance="no") "LogQsK_"//trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "MolNr"//trim(vKinFas(i)%NamKF),T_;   end do
  write(f,*)
end subroutine EnTeteFMnk

subroutine EnTeteFRate(f)
  use M_Global_Vars,only: vKinFas
  !
  integer,intent(in)::f
  !
  integer::i,N
  !
  N=size(vKinFas)
  !
  write(f,'(4(A,A1))',advance="no") "iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !
  do i=1,N; write(f,'(A,A1)',advance="no") "Phi_"    //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "SurfR_"  //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "SurfKg_" //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "lgQsK_"  //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "ratQsK_" //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "ratAct_" //trim(vKinFas(i)%NamKF),T_; end do
  do i=1,N; write(f,'(A,A1)',advance="no") "Radius_" //trim(vKinFas(i)%NamKF),T_; end do
  !do i=1,N; write(f,'(A,A1)',advance="no") "SferNum_"//trim(vKinFas(i)%NamKF),T_; end do
  write(f,*)
  !
end subroutine EnTeteFRate

subroutine Dynam_Files_Close
  use M_Numeric_Tools,only: fNewtF,fNewtR
  use M_Dynam_Vars, only: fSavTime,fSavRate
  !! use M_Stockvar_Kinxim,only: LSTOCK,DEL_STOCKVAR,print_STOCKVAR !,SET_STOCKVAR,INIT_STOCKVAR
  !
  if(fSavTime>0) then; close(fSavTime); fSavTime= 0; end if
  if(fSavRate>0) then; close(fSavRate); fSavRate= 0; end if
  !
  if(fDynMol>0) then; write(fDynMol,*); close(fDynMol);  fDynMol= 0; end if 
  if(fDynAct>0) then; write(fDynAct,*); close(fDynAct);  fDynAct= 0; end if
  if(fDynGam>0) then; write(fDynAct,*); close(fDynGam);  fDynGam= 0; end if
  !
  if(fNewtF>0)  then ; close(fNewtF);  fNewtF=  0; end if
  if(fNewtR>0)  then ; close(fNewtR);  fNewtR=  0; end if
  !
  !if(fNewtFm>0)  then; close(fNewtFm); fNewtFm= 0; end if
  !if(fNewtRm>0)  then; close(fNewtRm); fNewtRm= 0; end if
  !
  if(fDynEle>0)  then ; close(fDynEle) ; fDynEle= 0 ; end if 
  if(fDynMnK>0)  then ; close(fDynMnK) ; fDynMnK= 0 ; end if
  !
  if(fDynElement>0)  then; close(fDynElement); fDynElement=   0; end if 
  
  !---del table to store the time evolution of species 
  !--- stockage solution
  !! if (LSTOCK) then
  !!    if ( DebugCoores ) call Print_Stockvar('output.var');
  !!    call Del_Stockvar
  !! end if
  
end subroutine Dynam_Files_Close

subroutine Dynam_Files_OpenLogs(f1,f2,f3,f4)
  use M_IOTools,only: GetUnit
  use M_Files,  only: DirOut,DirLog,Files_Index_Write
  use M_Global_Vars,only: vEle
  !! use M_Basis_Vars, only: iCpnEle
  use M_Dynam_Vars,only: vCpnBox
  use M_Dynam_Vars, only: TUnit
  !
  integer,intent(out):: f1,f2,f3,f4
  !
  integer:: I
  !
  !--------------------------------------------------------------------!
  call GetUnit(f1)
  !!open(f1,file=trim(DirLog)//"calcdyn_totout.log")
  open(f1,file=trim(DirOut)//"_bilans.restab")
  !
  call Files_Index_Write(fHtm,&
  !! & trim(DirLog)//"calcdyn_totout.log",& !
  & trim(DirOut)//"_bilans.restab", &
  & "DYNAMIC/LOG: balances on elements")
  !
  write(f1,'(2(A,A1))',advance="no") "STEP",T_,"TIME/"//TUnit,T_
  do I=1,size(vCpnBox)
    write(f1,'(A,A1)',advance="no") &
    & trim(vEle(vCpnBox(I)%iEle)%NamEl)//'inj',T_; end do 
  do I=1,size(vCpnBox)
    write(f1,'(A,A1)',advance="no") &
    & trim(vEle(vCpnBox(I)%iEle)%NamEl)//'bilan',T_; end do 
  do I=1,size(vCpnBox)
    write(f1,'(A,A1)',advance="no") &
    & trim(vEle(vCpnBox(I)%iEle)%NamEl)//'diff',T_; end do 
  write(f1,*)
  !
  !--------------------------------------------------------------------!
  call GetUnit(f2)
  !open(f2,file=trim(DirLog)//"calcdyn_steps.log")
  open(f2,file=trim(DirOut)//"_calcdyn_steps.log")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_calcdyn_steps.log",&
  & "DYNAMIC/LOG: iMaxDelta,VarMax,iDo1,iDo2,etc.")
  !
  write(f2,'(14(A,A1))') &
  & ".MaxVal",T_,"iStep",T_,"Time",T_,"dTime",T_,"PhiF",T_,&
  & "iMaxDelta",T_,"VarMAx",T_,&
  & "iDo1",T_,"iDo2",T_,"Newt_iDo",T_,&
  & "Newt_iErr",T_,"NewtErrF",T_,"NewtErrX",T_,"NewtErrG",T_
  !
  !--------------------------------------------------------------------!
  call GetUnit(f3)
  !open(f3,file=trim(DirLog)//"calcdyn_maxvars.log")
  open(f3,file=trim(DirOut)//"_calcdyn_maxvars.log")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_calcdyn_maxvars.log",&
  & "DYNAMIC/LOG: MaxAqu,MaxMin,MaxAquRel,MaxMinRel,...")
  !
  write(f3,'(9(A,A1))') &
  & ".MaxAqu",T_,".MaxMin",T_,&
  & "iStep",T_,"Time",T_,"dTime",T_,&
  & "MaxAqu",T_,"MaxMin",T_,&
  & "MaxAquRel",T_,"MaxMinRel",T_
  !
  !--------------------------------------------------------------------!
  call GetUnit(f4)
  open(f4,file=trim(DirLog)//"calcdyn_satur.log")
  !
  return
end subroutine Dynam_Files_OpenLogs

subroutine Dynam_Files_WriteFasAff(iStep,Time,vLnAct,vLnBuf)
!.write affinity of pure phases (i.e. non-aqu'species)
  use M_IOTools
  use M_Numeric_Const,only: Ln10
  use M_Files,        only: DirOut,Files_Index_Write
  use M_Basis_Vars,   only: vOrdPr,tNuFas,nCi
  use M_Global_Vars,  only: vFas,vSpc
  !
  integer, intent(in):: iStep
  real(dp),intent(in):: Time
  real(dp),intent(in):: vLnAct(:)
  real(dp),intent(in):: vLnBuf(:)
  !
  integer :: iFs,nFs,nCp,nCx
  real(dp):: X
  !
  nFs= size(vFas)
  nCp= size(vOrdPr)
  nCx= nCp -nCi
  !
  if(iStep>9999) return
  !
  if (nFs>0) then
    !
    !--------------------------------------- open fDynQsK, write en tete
    if(fDynQsK==0) then
      call GetUnit(fDynQsK)
      open(fDynQsK,file=trim(DirOut)//"_minqsk.restab")
      call Files_Index_Write(fHtm, &
      & trim(DirOut)//"_minqsk.restab", &
      & "DYNAMIC: Q/K of all phases from database")
      write(fDynQsK,'(2(A,A1))',advance= 'NO') "Step",T_,"Time",T_ !,"pH",T_
      do iFs=1,nFs
        write(fDynQsK,'(A,A1)',advance= 'NO') trim(vFas(iFs)%NamFs),T_
      end do
      write(fDynQsK,*)
    end if
    !------------------------------------------------------------------!
    !
    write(fDynQsK,'(I5,A1,G15.6,A1)',advance= 'NO') iStep,T_,Time,T_
    do iFs=1,nFs
      X= &
      & vFas(iFs)%Grt &
      - dot_product( tNuFas(iFs,1:nCi), &
      &              vLnAct(1:nCi) +vSpc(vOrdPr(1:nCi))%G0rt )
      if(nCx>0) &
      & X= X &
      &  - dot_product( tNuFas(iFs,nCi+1:nCp), &
      &              vLnBuf(1:nCx) +vSpc(vOrdPr(nCi+1:nCp))%G0rt )
      !
      write(fDynQsK,'(F12.3,A1)',advance="no") -X /Ln10,T_
      !                                       log10(QsK)= - Affinity /Ln10
    end do
    write(fDynQsK,*)
  end if
end subroutine Dynam_Files_WriteFasAff    

end module M_Dynam_Files
