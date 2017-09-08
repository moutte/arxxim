module M_GEM_Write
  !
  use M_Kinds
  use M_Trace, only: fTrc,iDebug,T_,Stop_
  !
  implicit none
  !
  private
  !
  public:: GEM_Path_Write_Phases
  public:: GEM_Path_Write_Mixtures
  public:: GEM_Write_Log_Entete
  public:: GEM_Single_Write_Mixtures_Moles
  public:: GEM_Single_Write_Mixtures_Title
  !
contains

subroutine GEM_Path_Write_Phases( &
& DimPath, &
& vFasIsPresent, &
& vSimplex_Ok, &
& tResult, &
& tVolume, &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  use M_T_MixModel,   only: T_MixModel
  use M_GEM_Vars,     only: T_SavModel
  !
  use M_Global_Vars,only: vSpc,vFas,vMixModel
  use M_GEM_Vars,   only: vCpnGEM
  !---------------------------------------------------------------------
  integer,       intent(in):: DimPath
  logical,       intent(in):: vFasIsPresent(:)
  logical,       intent(in):: vSimplex_Ok(DimPath)
  real(dp),      intent(in):: tResult(:,:)
  real(dp),      intent(in):: tVolume(:,:)
  type(T_SavModel),intent(in):: tResultMix(:,:)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K
  integer :: nC,nF,nP,nFpur,nFmix
  integer :: fMol,fVol,fGrm,fMov
  !real(dp):: vX(MaxPole)
  real(dp):: Tot,X,Y
  type(T_SavModel):: Fas0
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  !---------------------------------------------------------------------
  nC= size(vCpnGEM)
  nF= size(vFasIsPresent)
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  call GetUnit(fMol)
  open(fMol,file=trim(DirOut)//"_phase_mole.restab")
  call  Write_Title(fMol)
  !
  call GetUnit(fMov)
  open(fMov,file=trim(DirOut)//"_phase_molarvol.restab")
  call  Write_Title(fMov)
  !
  call GetUnit(fVol)
  open(fVol,file=trim(DirOut)//"_phase_cm3.restab")
  call  Write_Title(fVol)
  !
  call GetUnit(fGrm)
  open(fGrm,file=trim(DirOut)//"_phase_grams.restab")
  call  Write_Title(fGrm)
  !
  K=0
  do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !
      call Write_Left(fMol)
      call Write_Left(fMov)
      call Write_Left(fVol)
      call Write_Left(fGrm)
      !
      !! Tot=sum(tResult(1:nC,iPath))
      !! do iFs=1,nC
      !!   write(fMol,'(G15.6,A1)',advance="no") tResult(iFs,iPath)/Tot,T_
      !! end do
      
      !-----------------------------------------------------mole numbers
      do iFs=1,nF != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) &
        & write(fMol,'(G15.6,A1)',advance="no") &
        ! nr'moles phase          *nr'oxygen in formula
        ! & tResult(nC+iFs,iPath)*tFormula(iFs,1),T_
        & tResult(2+nC+iFs,iPath),T_
      end do
      write(fMol,*)
      !----------------------------------------------------/mole numbers

      !---------------------------------------------volume/cm3, weight/g
      do iFs=1,nFpur != (nC+3:nC+nF+2)
        if(vFasIsPresent(iFs)) then
          X= tResult(2+nC+iFs,iPath)
          Y= tVolume(iFs,iPath)
          write(fVol,'(G15.6,A1)',advance="no") X*Y*1.D6,T_
          write(fGrm,'(G15.6,A1)',advance="no") X*vFas(iFs)%WeitKg*1.D3,T_
          write(fMov,'(G15.6,A1)',advance="no") Y*1.D6,T_
        end if
      end do
      do iFs=1,nFmix
        if(vFasIsPresent(nFpur +iFs)) then
        
          Fas0= tResultMix(iFs,iPath) != mole nr of phase iFs
          MM= vMixModel(Fas0%iModel)
          X= 0.D0
          Y= 0.D0
          do J=1,MM%nPole
            do I=1,MM%nPole
              X= X + vSpc(MM%vIPole(I))%V0     * Fas0%tXPole(J,I) !volume/cm3
              Y= Y + vSpc(MM%vIPole(I))%WeitKg * Fas0%tXPole(J,I) !weight/g
            end do
          end do
          
          !X= dot_product( &
          !& tResultMix(iFs,iPath)%tXPole(1,:), &
          !& tResultMix(iFs,iPath)%vVol0(:))
          
          write(fVol,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*X*1.D6,T_
          write(fGrm,'(G15.3,A1)',advance="no") &
          & tResult(2+nC+nFpur+iFs,iPath)*Y*1.D3,T_
          
        end if
      end do
      write(fVol,*)
      write(fGrm,*)
      write(fMov,*)
      !--------------------------------------------/volume cm3, weight/g
    end if !!if(vSimplex_Ok(iPath))

  end do !!do iPath
  !
  write(fMol,'(A1)') "_"
  !
  close(fMol)
  close(fVol)
  close(fGrm)
  !
  if(iDebug>2) &
  & print '(2A)', "Phases: Results in ",trim(DirOut)//"_phase_xx.restab"
  !
contains

  subroutine Write_Left(FF)
    integer,intent(in):: FF
    integer:: i
    !
    write(FF,'(I3,A1)', advance="no") K,T_    != count
    write(FF,'(2(F12.3,A1))',advance="no") &  !
    & tResult(1,iPath), T_, &                 != TdgC
    & tResult(2,iPath), T_                    != Pbar
    do i=1,nC != (3:nC+2)                     != components
      write(FF,'(G15.6,A1)',advance="no") tResult(2+i,iPath),T_
    end do
    !
  end subroutine Write_Left

  subroutine Write_Title(FF)
    integer,intent(in):: FF
    !--------------------------------------------------------title lines
    write(FF,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
    !
    !--component names
    do i=1,nC
      write(FF,'(A,A1)',advance="no") trim(vCpnGEM(i)%NamCp)//"_cpn", T_
    end do
    !
    !--pure phase names
    do iFs=1,nFpur
      if(vFasIsPresent(iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vFas(iFs)%NamFs), T_
    end do
    !--mixture phase names
    do iFs=1,nFmix
      if(vFasIsPresent(nFpur +iFs)) &
      & write(FF,'(A,A1)',advance="no") trim(vMixModel(iFs)%Name), T_
    end do
    write(FF,*)
    !-------------------------------------------------------/title lines
  end subroutine Write_Title

end subroutine GEM_Path_Write_Phases

subroutine GEM_Single_Write_Mixtures_Title
  use M_IoTools,      only: GetUnit
  use M_T_MixModel,   only: T_MixModel,MaxPole,Mix_Molecular
  !
  use M_Global_Vars,  only: vSpc,vFas,vMixModel
  !
  integer :: iFs,I,P,nP,nPP
  integer :: F
  type(T_MixModel):: MM
  !
  call GetUnit(F)
  open(F,file="tmp_mixmodels_names.tab")
  !
  do iFs=1,size(vMixModel)
    !
    !if(vFasIsPresent(nFpur +iFs)) then
    !
    MM= vMixModel(iFs)
    nP= MM%nPole
    !if(Fas0%iModel /= 0) then
    if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
      nPP= 1
    else
      nPP= MM%nPole
    end if
    do I=1,nPP
      write(F,'(A,A1)',advance="NO") trim(MM%Name),T_
      do P=1,nP
        if(MM%vIPole(P)>0) &
        & write(F,'(A,A1)',advance="NO") &
        & trim(vSpc(MM%vIPole(P))%NamSp),T_
      end do
    end do
    !
    !end if
    !
  end do
  !
  close(F)
end subroutine GEM_Single_Write_Mixtures_Title

subroutine GEM_Single_Write_Mixtures_Moles(vSavModel)
  use M_IoTools,      only: GetUnit
  use M_Global_Vars,  only: vSpc,vFas,vMixModel
  use M_T_MixModel,   only: T_MixModel,MaxPole,Mix_Molecular
  use M_GEM_Vars,     only: T_SavModel
  !
  type(T_SavModel),intent(in):: vSavModel(:)
  !
  integer :: iFs,I,P,nP,nPP
  integer :: F
  type(T_MixModel):: MM
  type(T_SavModel):: Fas0
  !
  call GetUnit(F)
  open(F,file="tmp_mixmodels_moles.tab")
  !
  do iFs=1,size(vMixModel)
    !!if(vFasIsPresent(nFpur +J)) then
    Fas0= vSavModel(iFs)
    MM= vMixModel(Fas0%iModel)
    nP= MM%nPole
    !if(Fas0%iModel /= 0) then
    if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
      nPP= 1
    else
      nPP= MM%nPole
    end if
    ! print *,"Fas0%iModel=",Fas0%iModel
    ! print *,"nP=",nP ;  pause
    do I=1,nPP
      write(F,'(G15.6,A1)',advance="NO") Fas0%vMole(I),T_
      do P=1,nP
        if(MM%vIPole(P)>0) &
        & write(F,'(F7.3,A1)',advance="NO") &
        & Fas0%tXPole(I,P),T_
      end do
    end do
    !end if
    !!end if
  end do
  !
  close(F)
end subroutine GEM_Single_Write_Mixtures_Moles

subroutine GEM_Path_Write_Mixtures( &
& DimPath,       &
& vFasIsPresent, &
& vSimplex_Ok,   &
& TolX,          &
& tResult,       &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  use M_Files,        only: DirOut
  use M_Dtb_Const,    only: T_CK
  use M_IoTools,      only: GetUnit
  use M_Numeric_Tools,only: iMaxLoc_R
  !
  use M_T_Phase,      only: T_Phase
  use M_T_MixModel,   only: T_MixModel,MaxPole,Mix_Molecular
  !
  use M_Global_Vars,  only: vSpc,vFas,vMixModel
  use M_GEM_Vars,     only: T_SavModel
  use M_GEM_Vars,     only: vCpnGEM
  !---------------------------------------------------------------------
  integer, intent(in):: DimPath
  logical, intent(in):: vFasIsPresent(:)
  logical, intent(in):: vSimplex_Ok(DimPath)
  real(dp),intent(in):: TolX ! <-MixMinim_TolX
  real(dp),intent(in) :: tResult(:,:) ! for T,P values
  type(T_SavModel),intent(inout):: tResultMix(:,:)
  !---------------------------------------------------------------------
  integer :: iPath,iFs,I,J,K,P,Q
  integer :: nP,nFmix,nFPur,nPP
  integer :: FMIX
  real(dp):: vX(MaxPole)
  real(dp):: Tot
  type(T_SavModel):: Fas0,Fas1,SavModelZero
  type(T_MixModel):: MM
  !
  logical,parameter:: ComputeXMean= .true. !.false.
  !---------------------------------------------------------------------
  call GetUnit(FMIX)
  open(FMIX,file=trim(DirOut)//"_mixtures.restab")
  !
  nFpur= size(vFas)
  nFmix= size(vMixModel)
  !
  SavModelZero%iModel=      0
  SavModelZero%nFas=        0
  SavModelZero%tXPole(:,:)= Zero
  SavModelZero%vMole(:)=    Zero
  SavModelZero%vGrt0(:)=    Zero
  SavModelZero%vVol0(:)=    Zero
  !
  !-----------------------------------------------------------title line
  write(FMIX,'(3(A,A1))',advance="no") "count",T_, "TdgC",T_, "Pbar",T_
  !
  do iFs=1,nFmix
    !
    if(vFasIsPresent(nFpur +iFs)) then
      !
      !Fas= tResultMix(iFs,iPath)
      MM= vMixModel(iFs)
      nP= MM%nPole
      !if(Fas0%iModel /= 0) then
      if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
        nPP= 1
      else
        nPP= MM%nPole
      end if
      do I=1,nPP
        write(FMIX,'(A,A1)',advance="NO") trim(MM%Name),T_
        do P=1,nP
          write(FMIX,'(A,A1)',advance="NO") &
          & trim(vSpc(MM%vIPole(P))%NamSp),T_
        end do
      end do
      !
    end if
    !
  end do
  write(FMIX,*)
  !----------------------------------------------------------/title line
  
  K=0
  doPath: do iPath=1,DimPath

    if(vSimplex_Ok(iPath)) then
      !
      K=K+1
      !------------------------------compute mean comp'n of each mixture
      !----------------------------------------------and its free energy
      if(ComputeXMean) then
        !
        do J=1,nFmix
          !
          Fas0=        tResultMix(J,iPath)
          Fas1=        SavModelZero
          Fas1%iModel= Fas0%iModel
          nP=          vMixModel(Fas0%iModel)%nPole
          Fas1%nFas=   nP
          !
          do P=1,Fas0%nFas
            do Q=1,Fas1%nFas
              !
              vX(1:nP)= abs( Fas0%tXPole(P,1:nP) - Fas0%tXPole(Q,1:nP) )
              !
              if(maxval(vX(1:nP)) < TolX *1.0D2) then
                !
                Tot= Fas0%vMole(P) +Fas1%vMole(Q)
                vX(1:nP)= Fas0%vMole(P) *Fas0%tXPole(P,1:nP) &
                &       + Fas1%vMole(Q) *Fas0%tXPole(Q,1:nP)
                Fas1%tXPole(Q,1:nP)= vX(1:nP) /Tot
                Fas1%vMole(Q)= Tot
                !
                exit
                !
              end if
              !
            end do
          end do
          !
          tResultMix(J,iPath)= Fas1
          !
        end do
        !
      end if
      !----------------------------------/compute mean comp'n of mixture

      !---------------------------------------write mixture compositions
      !-----------------------------------------------------in file FMIX
      write(FMIX,'(I3,A1)',     advance="no") K,T_   != count
      !
      write(FMIX,'(2(F7.2,A1))',advance="no") &      !
      & tResult(1,iPath), T_, &                      != TdgC
      & tResult(2,iPath), T_                         != Pbar

      do J=1,nFmix
        if(vFasIsPresent(nFpur +J)) then
          Fas0= tResultMix(J,iPath)
          MM= vMixModel(Fas0%iModel)
          nP= MM%nPole
          !if(Fas0%iModel /= 0) then
          if( MM%Model==Mix_Molecular .and. MM%NMarg==0 ) then
            nPP= 1
          else
            nPP= MM%nPole
          end if
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          do I=1,nPP
            write(FMIX,'(G15.6,A1)',advance="NO") Fas0%vMole(I),T_
            do P=1,nP
              write(FMIX,'(F7.3,A1)',advance="NO") &
              & Fas0%tXPole(I,P),T_
            end do
          end do
          !end if
        end if
      end do
      !
      write(FMIX,*)
      !--------------------------------------/write mixture compositions
      !
    end if !!if(vSimplex_Ok(iPath))

  end do doPath
  !
  if(iDebug>2) &
  & print '(2A)', "Mixtures: Results in ",trim(DirOut)//"_mixtures.restab"
  !
end subroutine GEM_Path_Write_Mixtures

subroutine GEM_Write_Log_Entete(vFas)
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirOut
  use M_T_Phase,only: T_Phase
  !
  type(T_Phase),intent(in):: vFas(:)
  !
  integer:: fVol,fGrm
  !
  call GetUnit(fVol)
  open(fVol,file=trim(DirOut)//"_gibbs.log")
  write(fVol,'(A)') "g of all phases, relative to current assemblage"
  !
  call GetUnit(fGrm)
  open(fGrm,file=trim(DirOut)//"_mixcomp.log")
  write(fGrm,'(A)')   "properties of mixtures with min(g)<0 at each iteration"
  write(fGrm,'(A,/)') "mix.model.name, Xi, Gmin"
  !
end subroutine GEM_Write_Log_Entete

end module M_GEM_Write

!! !!unused!!
!! subroutine Simplex_WriteTable(iFil,M,N,tSimplex)
!!   use M_IoTools,     only: OutStrVec 
!!   ! use M_Simplex_Vars,only: tSimplex 
!!   integer,  intent(in):: iFil
!!   integer,  intent(in):: M,N
!!   real(dp), intent(in):: tSimplex(0:M+1,0:N)
!!   !
!!   real(dp),dimension(1:N+1)::V
!!   integer:: iLin !, iCol
!!   !
!!   do iLin=0,M+1
!!     V(1:N+1)=tSimplex(iLin,0:N)
!!     call OutStrVec(iFil,V,Opt_C="7")
!!   end do
!!   !
!!   write(iFil,'(A1)') "_"
!!   !
!! end subroutine Simplex_WriteTable


!Simplex, notes from NumRecp:
!
!linear programming= 
!maximize the linear function 
!  Z=  a(0,1).x(1) + ... + a(0,n).x(n)
!subject to 
!  primary constraints: a(0,1:n)>=0
!and simultaneously subject to 
!  M=  m1 + m2 + m3 additional constraints,
!  m1 of them of the form: a(i,1).x1 + ... + a(i,N).x(N) >= b(i)   //i=1..m1
!  m2 of them of the form: a(j,1).x1 + ... + a(j,N).x(N) <= b(j)   //j=m1+1..m1+m2
!  m3 of them of the form: a(k,1).x1 + ... + a(k,N).x(N) <= b(k)   //k=m1+m2+1..m1+m2+m3

!On output, the tableau a is indexed by two returned arrays of integers
!iposv(j) contains, for j= 1...M, 
!the number i whose original variable x(i) is now represented by row j+1 of a
!These are thus the left-hand variables in the solution
!(The first row of a is of course the z-row.)
!A value i > N indicates that the variable is a y(i) rather than an x(i), x(N+j)==y(j)
!Likewise, izrov(j) contains, for j= 1..N,
!the number i whose original variable x(i) is now a right-hand variable, represented by column j+1 of a. 
!These variables are all zero in the solution. 
!The meaning ofi > N is the same as above, 
!except that i >N +m1 +m2 denotes an artificial or slack variable 
!which was used only internally and should now be entirely ignored.


