BeginPackage["DetectorModel`"]

DetectorModel::usage = "DetectorModel[s1, s2, cc] - Requires user to input three arrays that follow the structure: s1 = {{power1(mW), singles1(raw counts)}, {power2(mW), singles2(raw counts)}, ... }. This is also the same for the coincidence data structure.";

Begin["`Private`"]

pEm[n_, \[Gamma]_, P_] := 
    (1 - \[Gamma] P) (\[Gamma] P)^n;

pDet[\[Eta]_, n_, k_: 1] := Block[{}, 
    (1 - Sum[Binomial[n, idet] \[Eta]^idet (1 - \[Eta])^(n - idet), {idet, 0, k - 1}])];

pReady[\[Gamma]_, P_, rep_, \[Eta]_, rt_, kready_]:= 
    ((-1 + P \[Gamma])/(-1 + P \[Gamma] (1 - \[Eta])^kready))^Floor[rep rt];

Singles[\[Gamma]_, P_, rep_, \[Eta]_, rt_]:= 
    (P rep \[Gamma] ((1 - P \[Gamma])/(1 + P \[Gamma] (-1 + \[Eta])))^Floor[rep rt] \[Eta])/(1 + P \[Gamma] (-1 + \[Eta]));    

DetProbsMultiSourceAll[n_,\[Gamma]_,P_,rep_,\[Eta]1_,\[Eta]2_,rt1_,rt2_]:=
    rep ((P \[Gamma] ((1-P \[Gamma])/(1+P \[Gamma] (-1+\[Eta]1)))^Floor[rep rt1] \[Eta]1 ((1-P \[Gamma])/(1+P \[Gamma] (-1+\[Eta]2)))^Floor[rep rt2] (-1+P^2 \[Gamma]^2 (-1+\[Eta]1) (-1+\[Eta]2)) \[Eta]2)/((1+P \[Gamma] (-1+\[Eta]1)) (1+P \[Gamma] (-1+\[Eta]2)) (-1+P \[Gamma] (-1+\[Eta]1) (-1+\[Eta]2))))^n;

deadTime = Total@{0.0000000503125,0.0000000378125,0.0000000253125,0.0000000378125}/Length@{0.0000000503125,0.0000000378125,0.0000000253125,0.0000000378125};

DetectorModel[singles1_,singles2_,cc_,rep0_]:=

    Block[{\[Tau],\[Gamma],\[Eta]1,\[Eta]2,detectorModel},

    detectorModel = NMinimize[{(Total@Table[Abs[Singles[\[Gamma],singles1[[i,1]]*10^-3,rep0,\[Eta]1,deadTime]-singles1[[i,2]]]+Abs[Singles[\[Gamma],singles2[[i,1]]*10^-3,rep0,\[Eta]2,deadTime]-singles2[[i,2]]]+
    Abs[DetProbsMultiSourceAll[1,\[Gamma],cc[[i,1]]*10^-3,rep0,\[Eta]1,\[Eta]2,deadTime,deadTime]-cc[[i,2]]],{i,1,Length[singles1]}])
    ,0<=\[Gamma]<=0.4,0.1<=\[Eta]1<=0.8,0.1<=\[Eta]2<=0.8}
    ,{\[Gamma],\[Eta]1,\[Eta]2}];

    {\[Tau]=detectorModel[[2,1,2]],
    \[Eta]1=detectorModel[[2,2,2]],
    \[Eta]2=detectorModel[[2,3,2]]}]


End[]
EndPackage[]