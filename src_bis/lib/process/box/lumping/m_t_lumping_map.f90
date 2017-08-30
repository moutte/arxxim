module M_T_Lumping_Map
  
  use M_Kinds
  use M_Trace
  
  implicit none
  !--
  private 

  type:: T_Lumping_Map 
     character(len=80):: cPhase
     character(len=80):: cCons
     character(len=80):: cSpecies
     !---
     integer :: iSpecies
     integer :: iPhase
     integer :: iCons
  end type T_Lumping_Map

  integer, parameter ::     P_WATER = 1 
  integer, parameter ::     P_GAS   = 2
  integer, parameter ::     P_OIL   = 3

  !// Public Functions
  public:: T_Lumping_Map 

  public:: Lumping_Map_Zero 
  public:: Lumping_Map_Index
  public:: Lumping_Map_Add
  
  public:: Lumping_Map_Print
  public:: Lumping_Map_Print_Header
  public:: Find_Lumping_Map_Species
  public:: Lumping_Map_IsPhaseWater
  
contains

  function Lumping_Map_IsPhaseWater(M) result (Ok)
    implicit none
    type(T_Lumping_Map), intent(in) :: M
    logical :: Ok
    !---
    Ok = ( M%iPhase .eq. P_WATER ) 
    
  end function Lumping_Map_IsPhaseWater
  
  !---
  
  subroutine Lumping_Map_Zero(M)
    !========================================
    ! Create an empty Lumping_Map Line
    !========================================
    implicit none
    type(T_Lumping_Map) :: M    
    !---    
    M%cPhase   = "none"
    M%cCons    = "none" 
    M%cSpecies = "none"
    !--
    M%iPhase   = 0
    M%iCons    = 0
    M%iSpecies = 0

  end subroutine Lumping_Map_Zero

  !---

  subroutine Find_Lumping_Map_Species(iSpc, IdxSpc, NumSpc, vLumpMap, sizeLumpMap)
    implicit none

    type(T_Lumping_Map) :: vLumpMap(:)
    integer, intent(in) :: sizeLumpMap
    integer             :: IdxSpc
    integer             :: iSpc
    !--
    integer             :: I
    logical             :: OkSpc
    integer             :: NumSpc
    !--
    IdxSpc=0
    NumSpc=0
    do I=1,sizeLumpMap
       OkSpc = ( iSpc == vLumpMap(I)%iSpecies)
       if ( OkSpc ) then
          IdxSpc=I
          NumSpc= NumSpc +1
          exit
       end if
    enddo

  end subroutine Find_Lumping_Map_Species

  !---

  function Lumping_Map_Index(cPhase, cCons, vLumpMap, sizeLumpMap)  result(Idx)
    implicit none

    character(len=*),             intent(in):: cPhase, cCons
    type(T_Lumping_Map) :: vLumpMap(:)
    integer, intent(in) :: sizeLumpMap
    integer             :: Idx
    !--
    integer             :: I
    logical             :: OkPhase, OkCons
    !--
    Idx=0
    do I=1,sizeLumpMap
       OkPhase = ( trim(cPhase) == trim( vLumpMap(I)% cPhase ) )
       OkCons  = ( trim(cCons)  == trim( vLumpMap(I)% cCons  ) )
       if ( OkPhase .and. OkCons ) then
          Idx=I
          exit
       end if
    enddo

  end function Lumping_Map_Index

  !---

  subroutine Lumping_Map_Add( sPhase, sCons, sSpecies, vLumpMap, Size_LumpMap, vSpc)
    use M_T_Species
    use M_StringUtils
    implicit none

    type(T_Lumping_Map) :: vLumpMap(:)
    type(T_Species) :: vSpc(:)
    
    integer :: Size_LumpMap
    character(len=*), intent(in)  :: sPhase, sCONS, sSpecies
    character(len=80) :: cPhase, cCONS, cSpecies
    integer :: iPhase, iSpecies, iCons
    integer :: iMap
    integer :: ZSpecies
    !---    

    !// Copy and Clean Fields    
    cPhase   = trim(adjustl(sPhase))
    cCons    = trim(adjustl(sCons))
    cSpecies = trim(adjustl(sSpecies))

    !// Check Phase-Cons         
    select case(trim(cPhase))
    case('WATER') ; iPhase = P_WATER  
    case('OIL')   ; iPhase = P_OIL
    case('GAS')   ; iPhase = P_GAS
    case default
       call FATAL_("Wrong Lumping Phase = " //trim(cPhase))
    end select

    !// Get Map Index    
    iMap  = Lumping_Map_Index(cPhase, cCons, vLumpMap, Size_LumpMap)
    if ( (trim(cCons) == 'C-LUMP')) then
       ! mutiple definition allowed only for C-LUMP
    else
       if (iMap>0) then
          call Fatal_("Lumping already defined for ["//trim(cPhase)//'|'//trim(cCons)//"]")
       end if
    end if
    
    Size_LumpMap = Size_LumpMap + 1
    if ( Size_LumpMap  > size(vLumpMap) ) then 
       call Fatal_("Maximal size for Lumping Map")
    end if
    iMap = Size_LumpMap

    !// Get Species Name   
    iSpecies = Species_Index(sSpecies, vSpc)
        
    if (iSpecies == 0) call FATAL_("Species Unknown = " //trim(sSpecies))
    
    ZSpecies = vSpc(iSpecies)%Z
    if ( (.not.(trim(cCons) .eq. 'C-LUMP')) .and. (.not.(ZSpecies .eq. 0))) then    
      call Fatal_("Charged Species not Allowed in Lumping : "//trim(sSpecies))
    end if
    
    !// Fill LumpingMap    
    vLumpMap(iMap) % cPhase    = cPhase
    vLumpMap(iMap) % cCons     = cCons
    vLumpMap(iMap) % cSpecies  = cSpecies
    vLumpMap(iMap) % iSpecies  = iSpecies
    vLumpMap(iMap) % iPhase    = iPhase

    !// Reference Index for Cons
    iCons = String_Index(vLumpMap(1:Size_LumpMap)%cCons , cCons)
    vLumpMap(iMap) % iCons     = iCons    

  end subroutine Lumping_Map_Add

  !---

  subroutine Lumping_Map_Print_Header(F)
    !========================================
    ! Print a Lumping Map Header
    !========================================
    implicit none
    integer :: F

    Write(*,'(A6, 3A16, 3A9)' ) &
         'iMap    ',&
         'cPhase          ', 'cCons            ', 'cSpecies        ' , &
         'iPhase  ', 'iCons     ', 'iSpecies '

    Write(*,'(A6, 3A16, 3A9)' ) &
         '----- ', &
         '--------------  ', '--------------  ', '--------------  ', &
         '-------- ', '-------- ', '-------- '

  end subroutine Lumping_Map_Print_Header

  !---

  subroutine Lumping_Map_Print(M,iMap,F)
    !========================================
    ! Print a Lumping Map Line Values
    !========================================
    implicit none
    integer :: F
    type(T_Lumping_Map) :: M
    integer :: iMap
    !---    
    Write(*,'(I4,2X, 3A16, 3(I6,3X))' ) &
         iMap, &
         M%cPhase, M%cCons, M%cSpecies , &
         M%iPhase, M%iCons, M%iSpecies

  end subroutine Lumping_Map_Print

  !---

end module M_T_Lumping_Map
