module M_FILES_VARS
  implicit none
  private

  character(len=80),public:: cTitle
  
  ! Input File Name
  character(len=80),public:: NamFInn
  character(len=80),public:: NamFLogK,NamFEle,NamFKin,NamFSol,NamFPtz
  character(len=80),public:: NamDtbAqu,NamDtbMin,NamDtbMlt

  character(len=80),public:: DirOut,DirLog
  character(len=80),public:: DirDtbOut,DirDtbLog 

end module M_FILES_VARS
