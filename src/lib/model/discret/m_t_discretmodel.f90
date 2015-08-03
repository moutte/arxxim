MODULE M_T_DiscretModel
  !.structures for describing the "discretization" of a mixture phase
  !.to an array of phases of fixed composition (so-called "pure phase")
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_
  !
  PRIVATE
  !
  PUBLIC:: T_DiscretModel
  PUBLIC:: T_DiscretParam
  !
  INTEGER,PUBLIC:: MaxDiscret= 100
  INTEGER,PUBLIC:: MinDiscret= 1
  !
  TYPE:: T_DiscretModel
    CHARACTER(LEN=15):: Name
    INTEGER:: Dim1,Dim2,DimTot
    ! for binary discretization: Dim2=0
    INTEGER:: iMix !-> index in vMixModel
    INTEGER:: P1,P2,P3 !-> indexes of end members in mixing model
  ENDTYPE T_DiscretModel
  !
  TYPE:: T_DiscretParam
    INTEGER:: iModel !-> points to a model in vDiscretModel
    INTEGER:: I,J,K  !-> coordinates of discretization point
    INTEGER:: iSpc   !-> points back to index in current vSpc
  ENDTYPE T_DiscretParam
  !
ENDMODULE M_T_DiscretModel
