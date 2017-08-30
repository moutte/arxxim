module m_T_EosIfp_vState
  
  use M_Kinds,only: dp

  implicit none

  !// Computational TPZ Vectors 
  type:: T_EosIfp_vState

    integer :: NCONS, LVTH
    
    real(dp), pointer :: T(:), P(:) ! (LVTH) ! TP conditions
    real(dp), pointer :: X(:,:)     ! (LVTH, NCONS) ! composition 

    !------ Computed Properties 
    real(dp), pointer :: Z(:)         ! (LVTH) ! covolumes (LVTH)
    real(dp), pointer :: LnGamma(:,:) ! (LVTH, NCONS) ! Ln of Fugacity Coefficient

  end type T_EosIfp_vState

  !----

  public :: EosIfp_vState_new
  public :: EosIfp_vState_delete
  public :: EosIfp_vState_ComputeFugacity
  public :: EosIfp_vState_set

contains 

subroutine EosIfp_vState_new( T, NCONS, LVTH)

  implicit none
  type(T_EosIfp_vState), intent(inout) :: T
  integer, intent(in) :: NCONS, LVTH

  !-------------

  ! sizes
  T% NCONS = NCONS
  T% LVTH  = LVTH

  ! arrays
  allocate( T% T(LVTH) )
  allocate( T% P(LVTH) )
  allocate( T% X(LVTH, NCONS) )
  allocate( T% Z(LVTH) )
  allocate( T% LnGamma(LVTH, NCONS) )

  T% T = 0.D0
  T% P = 0.D0
  T% X = 0.D0
  T% Z = 0.D0
  T% LnGamma = 0.D0

end subroutine EosIfp_vState_new

!----

subroutine EosIfp_vState_delete( T )

  implicit none
  type(T_EosIfp_vState), intent(inout) :: T

  ! arrays
  deallocate( T% T )
  deallocate( T% P )
  deallocate( T% X )
  deallocate( T% Z )
  deallocate( T% LnGamma )

end subroutine EosIfp_vState_delete

!---

subroutine EosIfp_vState_Set( T, vTemperature, vPressure, vComposition )

  implicit none
  type(T_EosIfp_vState), intent(inout) :: T
  real(dp), intent(in) :: vTemperature(:)
  real(dp), intent(in) :: vPressure(:)
  real(dp), intent(in) :: vComposition(:,:)
  integer :: NCONS, LVTH

  !---
  NCONS = T% NCONS
  LVTH  = T% LVTH
  !---
  T% T(1:LVTH) = vTemperature
  T% P(1:LVTH) = vPressure
  T% X(1:LVTH, 1:NCONS) = vComposition

end subroutine EosIfp_vState_Set

!-----
subroutine EosIfp_vState_Print( T, funit )

  integer :: funit
  type(T_EosIfp_vState), intent(in) :: T
  integer :: NCONS, LVTH
  integer :: I, J

  !---
  NCONS = T% NCONS
  LVTH  = T% LVTH

  write(funit,*) ">> vState Print"
  do i=1, LVTH
    write(funit,*) "------------------------"
    write(funit,*) "Point : ", i
    write(funit,*) "P =", T%P(i)
    write(funit,*) "T =", T%T(i)
    do j=1, NCONS
      write(funit, *) "X (",j,") = ", T%X(i,j)
    end do
    write(funit,*) "------------------------"
    write(funit, *) "Z =", T%Z(i)
    do j=1, NCONS
      write(funit, *) "LnGamma (" , j," ) =", T%LnGamma(i,j)
    end do
  end do

  write(funit, *)
  
end subroutine EosIfp_vState_Print

subroutine EosIfp_vState_ComputeFugacity( T, T_Const, T_Param, NCONS, LVTH )

  use M_T_EosIfp_Const
  use M_T_EosIfp_Param
  use M_EosIfp_Calc, only :  EosIfp_VCUBZ, EosIfp_VFUGA
  implicit none
  type(T_EosIfp_Param), intent(in) :: T_Param
  type(T_EosIfp_Const), intent(in) :: T_Const
  type(T_EosIfp_vState), intent(inout) :: T

  integer, intent(in) :: NCONS
  integer, intent(in) :: LVTH
  !!!--------------- computation tables
  !----------------------------- intermediate parameters
  real(dp) :: COEFG(NCONS, NCONS)
  real(dp) :: COEF1(NCONS, NCONS)
  real(dp) :: COEF2(NCONS, NCONS)
  real(dp) :: PBI(NCONS)
  real(dp) :: POM(NCONS)
  real(dp) :: XVAL, OME
  !---------------------------- local copy of T_Param 
  real(dp) :: TC(NCONS), PC(NCONS)
  real(dp) :: OMEGA(NCONS)
  real(dp) :: COINT(NCONS,NCONS)
  real(dp) :: OMGA(NCONS)
  real(dp) :: OMGB(NCONS)
  real(dp) :: AIC(NCONS)
  real(dp) :: DLTA1, DLTA2
  !----------------------------- outputs 
  real(dp) :: PB(LVTH), PA(LVTH), GB(LVTH), AKJGA(LVTH),  CZ1(LVTH)
  real(dp) :: SAIJXJ(NCONS, LVTH)
  real(dp) :: RPB(LVTH), RBBRTD(LVTH)
  logical :: TRAINV(LVTH)
  real(dp) :: CZBAS2(LVTH), CZALS3(LVTH), CZ2S3(LVTH), DISCRI(LVTH)
  real(dp) :: Z1(LVTH), Z2(LVTH), Z3(LVTH)
  !----------------------------- outputs
  real(dp) :: ZPD2GB(LVTH),RZPDGB(LVTH),ZD12BL(LVTH) 
  real(dp) :: ZD12SB(LVTH),RZMGB(LVTH),ZM1SB(LVTH)

  !------------------ intermediate outputs
  real(dp) :: FUG(LVTH, NCONS)
  real(dp) :: Z(LVTH)
  real(dp) :: XAF(0:3,4)  
  !!!!!
  integer :: I, J, IT, JT
  !---
  real(dp), parameter :: R = 8.31439D0
  !------------------ remplissage de DLTA1 et DLTA2
  !  + valeurs par defaut de OMGA0 et OMGB0
  !------------------

  !------------------ copy parameters to standard arrays

  TC    = T_Param%TC(:)
  PC    = T_Param%PC(:)
  OMGA  = T_Param%OMGA(:)
  OMGB  = T_Param%OMGB(:)
  OMEGA = T_Param%OMEGA(:)
  COINT = T_Param%COINT(:,:)
  XAF   = T_Const%XAF(:,:)
  DLTA1 = T_Param%DLTA1
  DLTA2 = T_Param%DLTA2

  !// Model to index !!! special select case action ... see IT,JT
  select case( T_Param%MODEL)
  case ("PR1")  ; IT = 1
  case ("PR2")  ; IT = 2
  case ("SRK1") ; IT = 3
  case ("SRK2") ; IT = 4
  case default  
    write(*,*) "Error Bad Model <",  T_Param%MODEL, ">"
    stop "Fatal Error"
  end select

  do I=1,NCONS 
    AIC(I)= OMGA(I)*(R*TC(I))**2/PC(I)
    PBI(I)= OMGB(I)*R*TC(I)/PC(I)
    OME=OMEGA(I)
    JT=IT
    if ((IT==2.or.IT==4).and.OME.LE..49) JT=IT-1
    POM(I)=XAF(0,JT)+(XAF(1,JT)+(XAF(2,JT)+XAF(3,JT)*OME)*OME)*OME
  end do

  do I=1, NCONS
    do J=1, NCONS
      XVAL=(1.- COINT(I,J))*SQRT(AIC(I)*AIC(J))
      COEF1(I,J)= XVAL*(1.+ POM(I))*(1.+POM(J))
      COEF2(I,J)=-XVAL*( (1.+ POM(I)) * POM(J) / SQRT(TC(J)) &
        + (1.+ POM(J)) * POM(I) / SQRT(TC(I)) )
      COEFG(I,J)= XVAL*POM(I)*POM(J) / SQRT(TC(I) * TC(J))
    end do
  end do

  !-------- Cubic Resolution for Z

  call EosIfp_VCUBZ ( &
      LVTH, NCONS    &                               ! taille
    , T%T, T%P, T%X  &                               ! variables
    , Z  &                                           ! sortie principale
    , COEFG, COEF1, COEF2, PBI, DLTA1, DLTA2 &       ! parametres
    , PB, PA, GB, AKJGA, CZ1, SAIJXJ, RPB, RBBRTD &  ! sorties masquees 
    , CZBAS2,CZALS3,CZ2S3,DISCRI,TRAINV,Z1,Z2,Z3  &  ! variables de calcul
    )                                               

  ! store Z result
  T %Z = Z(:)
  !!!!write(*,*) "Z=", Z

  !-------- Fugacity Computation


  call  EosIfp_VFUGA ( LVTH, NCONS &                   ! taille 
    , Z &                                              ! variables
    , FUG &                                            ! sortie principale
    , PBI, DLTA1, DLTA2 &                              ! parametres
    , PB, PA, GB, AKJGA, CZ1, SAIJXJ, RPB, RBBRTD &    ! parametres
    , ZPD2GB, RZPDGB, ZD12BL, ZD12SB, RZMGB, ZM1SB &   ! variables de calcul
    )

  ! store Fuga result
  T% LnGamma = FUG(:,:)

end subroutine EosIfp_vState_ComputeFugacity

end module m_T_EosIfp_vState
