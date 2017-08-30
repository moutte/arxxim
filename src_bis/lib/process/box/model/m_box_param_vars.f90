module M_Box_Param_Vars

  !------------------------------------------------------------
  ! Size Parameters for Module Box
  !------------------------------------------------------------

  implicit none
  public

  !// System size and structure  

  integer :: nSpc ! number of species
  
  integer :: nAq  ! number of aqueous species
  integer :: nMk  ! number of kinetic minerals      

  integer :: nCp  ! number of components
  integer :: nPr  ! number of primary aqueous species
  integer :: nAs  ! number of secondary aqueous species
 
end module M_Box_Param_Vars
