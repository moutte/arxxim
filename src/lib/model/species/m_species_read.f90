module M_SpeciesDtb_Read
!--
!-- routines for reading data for species from the databases
!--
  use M_Kinds
  use M_Trace,    only: iDebug,fTrc,fHtm,Stop_,T_,Warning_,Pause_
  use M_T_Species,only: T_SpeciesDtb
  implicit none
  !
  private
  !
  public:: SpeciesDtb_BuildLnk
  public:: SpeciesDtb_LnkToVec
  public:: LnkSpc_Build
  public:: Species_Read_AquSize
  public:: Species_Write_AquSize
  !
  public:: T_LnkSpc
  !
  type:: T_LnkSpc
    type(T_SpeciesDtb)      :: Value
    type(T_LnkSpc), pointer :: Next
  end type T_LnkSpc
  !
contains

subroutine Species_RemoveDoublons
  use M_T_Species,  only: T_Species,Species_Index
  use M_Global_Vars,only: vSpc
  !
  ! type(T_Species),allocatable:: vSpcNew(:)
  ! logical,allocatable:: vSkip(:)
  integer:: N,NNew,I

  N= size(vSpc)
  NNew= N
  do while(N>0)
    I= Species_Index(vSpc(N)%NamSp,vSpc)
    if(I<N) then
      print *,"replace species ",trim(vSpc(N)%NamSp), trim(vSpc(I)%NamSp)
      call Pause_
      vSpc(I)= vSpc(N)
      NNew= NNew -1
    end if
    N= N-1
  end do

  if(NNew<N) then
    print *,"Old/NewDim",N,NNew
    call Pause_
  end if

  return
end subroutine Species_RemoveDoublons

subroutine LnkSpc_Build(Init,Inputt,L,P)
  logical,           intent(in):: init
  type(t_SpeciesDtb),intent(in):: inputt
  type(t_Lnkspc),    pointer   :: l,p

  if(init) then
    nullify(l)
    allocate(l)
    nullify(l%next)
    l%Value=inputt
    p=>l
  else
    allocate(p%next)
    nullify(p%next%next)
    p%next%Value=inputt
    p=>p%next
  end if

  return
end subroutine LnkSpc_Build

subroutine SpeciesDtb_BuildLnk(LnkSpc,N)
!--
!-- build vSpc directly from HSV files
!-- implements the S%Model and S%Indx
!-- which link the species S with the database models
!--
  use M_Dtb_Vars, only: vDtbMinHkf,vDtbMinThr,vDtbAquHkf,vDtbLogKtbl,vDtbLogKAnl
  !
  type(T_LnkSpc),pointer    :: LnkSpc
  integer,       intent(out):: N
  !
  type(T_LnkSpc),pointer:: pCur
  type(T_SpeciesDtb):: S
  integer:: I
  !
  N= 0
  if(size(vDtbLogKtbl)>0) then
    do I=1,size(vDtbLogKtbl)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbLogKtbl
      S%DtbModel= "LOGKTBL"
      S%DtbTrace= trim(vDtbLogKtbl(I)%Num)
      call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    end do
    ! stop here,
    ! because logK discrete databases can't be mixed
    ! with bases of other type
    return !------------------------------------------------------return
  end if
  !
  N= 0
  if(size(vDtbLogKAnl)>0) then
    do I=1,size(vDtbLogKAnl)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbLogKAnl
      S%DtbModel= "LOGKANL"
      S%DtbTrace= trim(vDtbLogKAnl(I)%Num)
      call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    end do
    ! stop here,
    ! because logK analytic databases can't be mixed
    ! with bases of other type
    return !------------------------------------------------------return
  end if
  !
  N= 0
  !------------------------------------------------------aqueous species 
  !-------------------------------------place H2O as species Nr1 in list
  if(size(vDtbAquHkf)>0) then
    ! if there are aqueous species following the HKF equation of state,
    ! then the eos for H2O is necessarily the Haar-Gallagher-Kell model
    N=N+1
    S%Indx= 1
    S%DtbModel= "H2O_HGK" != Haar-Gallagher-Kell
    S%DtbTrace= "HAAR_etal"
    call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    !
    do I=1,size(vDtbAquHkf)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbAquHkf
      S%DtbModel= "AQU_HKF"
      S%DtbTrace= trim(vDtbAquHkf(I)%Num)
      call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    end do
  end if
  !
  !---------------------------------------then mineral (& gas ?) species
  if(size(vDtbMinHkf)>0) then
    do I=1,size(vDtbMinHkf)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbMinHkf
      S%DtbModel= trim(vDtbMinHkf(I)%Typ)//"_HKF"
      S%DtbTrace= trim(vDtbMinHkf(I)%Num)
      call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    end do
  end if
  !
  !---------------------------------------idem, for Theriak-like formats
  if(size(vDtbMinThr)>0) then
    do I=1,size(vDtbMinThr)
      N=N+1
      S%Indx=     I
      S%DtbModel= trim(vDtbMinThr(I)%Typ)//"_THR"
      S%DtbTrace= trim(vDtbMinThr(I)%Num)
      call LnkSpc_Build(N==1,S,LnkSpc,pCur)
    end do
  end if

  return
end subroutine SpeciesDtb_BuildLnk

subroutine SpeciesDtb_LnkToVec(Lnk,vSpc)
!--
!-- transfer the content of a linked list of species
!-- to an array of species,
!-- then deallocate the list
!--
  type(T_LnkSpc),    pointer    :: Lnk
  type(T_SpeciesDtb),intent(out):: vSpc(:)
  !
  type(T_LnkSpc),pointer:: pCur, pPrev
  type(T_SpeciesDtb):: S
  integer:: N

  N= 0
  pCur=> Lnk
  do while (associateD(pCur))
    S= pCur%Value
    N= N+1
    vSpc(N)= S
    pPrev=>  pCur
    pCur=>   pCur%next
    deallocate(pPrev)
  end do

  return
end subroutine SpeciesDtb_LnkToVec

subroutine Species_Read_AquSize(vSpc)
!--
!-- read size parameters of (some) aqueous species
!-- -> overwrite current value in vSpc
!--
  use M_IOTools
  use M_Files,    only: NamFInn
  use M_T_Species,only: T_Species,Species_Index
  !
  type(T_Species),intent(inout):: vSpc(:)
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios,i
  !
  logical,allocatable:: vIsPresent(:)

  if(count(vSpc(:)%Typ=="AQU")==0) return

  if(idebug>1) write(fTrc,'(/,A)') "< Species_Read_AquSize"

  call GetUnit(F)
  open(F,file=trim(NamFInn))

  allocate(vIsPresent(size(vSpc)))
  vIsPresent(:)= vSpc(:)%Z==0

  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    !
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)

    select case(W)
    !
    case("SPECIES.size")

      DoBlock: do

        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        call AppendToEnd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines

        select case(W)
          case("END","ENDSPECIES.size")
            exit DoFile
        end select

        i= Species_Index(W,vSpc)

        if(i>0) then
          if(vSpc(i)%Typ/="AQU") then
            call Warning_("Species "//trim(W)//" is not aqueous")
          else
            call LinToWrd(L,W,EoL)
            call WrdToReal(W,vSpc(i)%AquSize)
            vIsPresent(i)= .true.
            !print '(A,I3,1X,F7.2,1X,A)', &
            !& "AquSize ... ",i,vSpc(i)%AquSize,vSpc(I)%NamSp
            !pause_
          end if
        else
          call Warning_("Species "//trim(W)//" not in base")
        end if

      end do DoBlock

    end select

  end do DoFile
  close(F)

  if(idebug>1) then
    write(fTrc,'(A,/)') "aqu'size parameters"
    do i=1,size(vSpc)
      if(vSpc(i)%Typ=="AQU" .and. vSpc(i)%Z/=0) then
        if(vSpc(I)%AquSize >1.E-6) then
          write(fTrc,"(A,A1,F7.2)") &
          & vSpc(I)%NamSp, T_, vSpc(I)%AquSize
        else
          write(fTrc,"(A,A1,A)") &
          & vSpc(I)%NamSp, T_, " <- size parameter not defined"
        end if
      end if
    end do
  end if

  deallocate(vIsPresent)

  if(idebug>1) write(fTrc,'(A,/)') "</ Species_Read_AquSize"

  return
end subroutine Species_Read_AquSize

subroutine Species_Write_AquSize(vSpc)
!--
!-- write size parameters of charged aqueous species
!--
  use M_IOTools
  use M_T_Species,only: T_Species
  !
  type(T_Species),intent(in):: vSpc(:)
  !
  integer:: F,i
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Species_Write_AquSize"
  !
  call GetUnit(F)
  open(F,file="aqusize.tab")
  !
  write(F,'(A)') "SPECIES.size"
  !
  do i=1,size(vSpc)
    if(vSpc(i)%Typ=="AQU" .and. vSpc(i)%Z/=0) &
    & write(F,'(A,A1,F12.2)') vSpc(I)%NamSp, T_, vSpc(i)%AquSize
  end do
  !
  write(F,'(A)') "END"
  !
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Species_Write_AquSize"

  return
end subroutine Species_Write_AquSize

subroutine Canevas
!--
!-- template for file reading
!--
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Canevas"
  !
  call GetUnit(F)
  open(F,file="XXX")
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT")
      exit DoFile
    !
    case("block") !!!!!!canevas for reading one "block"
      !... I=0
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        call AppendToEnd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines
        !
        select case(W)
        case("ENDINPUT")
          exit DoFile
        case("END","ENDblock")
          exit DoBlock
          !-> in case reading data from several block..end blocks
          !case("END"); exit DoFile
          !!-> in case reading only the first available block
        case("XXX")
          !....
        case("YYY")
          !....
        end select
        !...
        !
      end do DoBlock
    !
    end select
    !
  end do DoFile
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Canevas"

  return
end subroutine Canevas

end module M_SpeciesDtb_Read

!! integer function SpeciesDtb_Index(Str,V)
!! !.position of species with %NamSp==Str in vSpc_(1:nSpc_)
!!   character(len=*),intent(in):: Str
!!   type(T_SpeciesDtb), intent(in):: V(:)
!!   !
!!   integer     ::I
!!   SpeciesDtb_Index=0
!!   if(size(V)==0) return
!!   I=0
!!   do
!!     I=I+1 !; if(idebug>1) write(fTrc,'(A)') vCpn(I)%SpName
!!     if(trim(Str)==trim(V(I)%NamSp)) then
!!       SpeciesDtb_Index=I
!!       exit
!!     end if
!!     if(I==size(V)) exit
!!   end do !if Str not found -> I=0
!!   return
!! end function SpeciesDtb_Index
!! 
!! logical function LnkSpc_Found(L,W,Spc)
!! !.find index of string W in linked list LnkSp
!! !.return corresponding species Spc
!!   !
!!   !integer::I_LnkSpc
!!   type(T_LnkSpc),    pointer    :: L
!!   character(len=*),  intent(in) :: W
!!   type(T_SpeciesDtb),intent(out):: Spc
!!   type(T_LnkSpc),    pointer    :: P,pPrev
!!   integer::I
!!   !
!!   P=>L; I=0
!!   LnkSpc_Found=.false.
!!   do while (associateD(P))
!!     I=I+1
!!     if(trim(w)==trim(P%Value%NamSp)) then
!!       LnkSpc_Found=.true.
!!       Spc=P%Value
!!       exit;
!!     end if
!!     pPrev=>P
!!     P=> P%next
!!   end do
!! end function LnkSpc_Found

