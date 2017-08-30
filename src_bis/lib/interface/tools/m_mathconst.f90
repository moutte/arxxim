module M_mathconst

  !! *******************************************************************
  !! * File name:     M_mathconst
  !! *
  !! * Purpose:       constantes mathematiques et physiques
  !! *                fondamentales
  !! *
  !! * Author:        Anthony Michel (anthony.michel@ifp.fr)
  !! *
  !! * Created:       2005
  !! *
  !! * Modification:  2005
  !! *
  !! *
  !! ********************************************************************

  implicit none
  public

  ! date and time
  ! ---
  real(kind=8), parameter :: second_in_day = 86400.d0               ! (3600*24)
  real(kind=8), parameter :: day_in_second = 1.157407407407407D-005 ! 1/second_in_day
  ! ---
  real (kind=8), parameter :: second_in_year = 31536000.D0            ! (3600*24*365)
  real (kind=8), parameter :: year_in_second = 3.170979198376459D-008 ! (1/second_in_year)
  ! ---
  real(kind=8), parameter :: day_in_year = 365.D0             ! 365
  real(kind=8), parameter :: year_in_day = 0.2739726027397260D-02 ! 1/day_in_year

  ! volumes
  ! ---
  real (kind=8), parameter :: liter_in_m3 = 1.d+3
  real (kind=8), parameter :: m3_in_liter = 1.d-3

  ! permeabilite
  real(kind=8), parameter :: Darcy_in_SI   = 0.10132499658281448750D+16
  real(kind=8), parameter :: SI_in_Darcy   = 9.869233D-016

  ! temperature
  real(kind=8), parameter :: K_origin_degC = 273.14d0

  ! pressions
  real(kind=8), parameter :: bar_in_Pa     = 1.d-5
  real(kind=8), parameter :: Pa_in_bar     = 1.d5



end module M_mathconst

  ! temperature
