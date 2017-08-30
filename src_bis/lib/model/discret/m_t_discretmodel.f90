module M_T_DiscretModel
  !.structures for describing the "discretization" of a mixture phase
  !.to an array of phases of fixed composition (so-called "pure phase")
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_
  !
  private
  !
  public:: T_DiscretModel
  public:: T_DiscretParam
  !
  integer,public:: MaxDiscret= 100
  integer,public:: MinDiscret= 1
  !
  type:: T_DiscretModel
    character(len=15):: Name
    integer:: Dim1,Dim2,DimTot
    ! for binary discretization: Dim2=0
    integer:: iMix !-> index in vMixModel
    integer:: P1,P2,P3 !-> indexes of end members in mixing model
  end type T_DiscretModel
  !
  type:: T_DiscretParam
    integer:: iModel !-> points to a model in vDiscretModel
    integer:: I,J,K  !-> coordinates of discretization point
    integer:: iSpc   !-> points back to index in current vSpc
  end type T_DiscretParam
  !
end module M_T_DiscretModel
