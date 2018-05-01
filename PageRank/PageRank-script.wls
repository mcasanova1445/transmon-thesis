#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["~/Documents/transmon-thesis/PageRank"]


(*g=Graph[{1->2,1->3,1->4,2->1,2->4,3->4,4->3}];*) (*Grafo*)
(*g=Graph[{1->2,1->3,2->4}];*) (*Grafo \[AAcute]rbol*)
(*g=Graph[{1->2,2->3,3->4,4\[Rule]1}];*) (*Grafo c\[IAcute]clico*)
(*g=Graph[{1->2,1->3,1->4}];*) (*Grafo estrella*)
(*g=Graph[{1->2,1->3,1->4,2->1,2->3,2->4,3->1,3->2,3->4,4->1,4->2,4->3}];*) (*Grafo completo*)
g=Graph[{1->2,2->3,3->1, 1->3,3->2,2->1, 1->4,2->4,3->4}] (*Loke directed*)
(*g=Graph[{1->2,2->3,3\[Rule]1, 1\[Rule]3,3\[Rule]2,2\[Rule]1, 1->4,2\[Rule]4,3\[Rule]4, 4\[Rule]1,4\[Rule]2,4\[Rule]3}]*) (*Loke unidirected*)

Export["grafo.png",g]


Em=AdjacencyMatrix[g]; 
For[j=1,j<=4,j++,
	OutDeg=Sum[Em[[j,i]],{i,1,4}];
	If[OutDeg!=0, Em[[j]]=Em[[j]]/OutDeg, Em[[j]]={1/4,1/4,1/4,1/4};]]
Em=Transpose[Em]; (*Matriz de adyacencia ponderada del grafo*)
"Matriz de adyacencia ponderada: " >> "output.txt"
Em >>> "output.txt"



\[Alpha]=0.85;
G=\[Alpha] Em + (1-\[Alpha])/Dimensions[Em][[1]] Table[1,{Dimensions[Em][[1]]},{Dimensions[Em][[1]]}]; (*Matriz de Google*)
"Matriz de Google: " >>> "output.txt"
G >>> "output.txt"



Ivl1={1/4}; (*PageRank inicial de cada nodo*)
Ivl2={1/4};
Ivl3={1/4};
Ivl4={1/4};
Iv={{1/4},{1/4},{1/4},{1/4}}; (*Vector de PageRank inicial*)
For[i=1,i<=30,i++,
	Iv=G.Iv; (*Iterar con la matriz de Google*)
	(*Iv=Iv/Norm[Iv];*)
	AppendTo[Ivl1,Iv[[1,1]]];
	AppendTo[Ivl2,Iv[[2,1]]];
	AppendTo[Ivl3,Iv[[3,1]]];
	AppendTo[Ivl4,Iv[[4,1]]];]
"PageRank obtenido recorriendo la caminata cl\[AAcute]sica: " >>> "output.txt"
Iv >>> "output.txt" (*Vector de PageRank final*)



Export["Nodo1-PRC.png",ListLinePlot[Ivl1,PlotRange->{0,1},DataRange->{0,30}]] (*Graficar evoluci\[OAcute]n del PageRank de cada nodo*)
Export["Nodo2-PRC.png",ListLinePlot[Ivl2,PlotRange->{0,1},DataRange->{0,30}]]
Export["Nodo3-PRC.png",ListLinePlot[Ivl3,PlotRange->{0,1},DataRange->{0,30}]]
Export["Nodo4-PRC.png",ListLinePlot[Ivl4,PlotRange->{0,1},DataRange->{0,30}]]



"PageRank dado por la funci\[OAcute]n de Mathematica:" >>> "output.txt"
PageRankCentrality[g,\[Alpha]] >>> "output.txt" (*Comparar con el PageRank dado por la funci\[OAcute]n de Mathematica*)



ket[k_,n_]:=Table[{i==k},{i,0,n-1}]/.{True->1,False->0}; (*Generador de kets: k \[Rule] estado, n \[Rule] dimensi\[OAcute]n*)



Subscript[\[Phi]x, j_]:=KroneckerProduct[ket[j-1,4],Sum[Sqrt[G[[k,j]]] ket[k-1,4],{k,1,4}]]; (*Estados A \[Rule] B*)
Subscript[\[Psi]y, j_]:=KroneckerProduct[Sum[Sqrt[G[[k,j]]] ket[k-1,4],{k,1,4}],ket[j-1,4]]; (*Estados B \[Rule] A*)
\[Phi]x0=Sum[Subscript[\[Phi]x, j],{j,1,4}]/2; (*Superposici\[OAcute]n de los estados A \[Rule] B*)
\[Psi]y0=Sum[Subscript[\[Psi]y, j],{j,1,4}]/2; (*Superposici\[OAcute]n de los estados B \[Rule] A*)
\[CapitalPi]x=Sum[Subscript[\[Phi]x, j].Subscript[\[Phi]x, j]\[ConjugateTranspose],{j,1,4}]; (*Proyector A \[Rule] B*)
\[CapitalPi]y=Sum[Subscript[\[Psi]y, j].Subscript[\[Psi]y, j]\[ConjugateTranspose],{j,1,4}]; (*Proyector B \[Rule] A*)
S=Sum[KroneckerProduct[ket[j,4],ket[k,4]].KroneckerProduct[ket[k,4],ket[j,4]]\[ConjugateTranspose],{j,0,3},{k,0,3}]; (*Operador SWAP*)
Ux=S.(2\[CapitalPi]x-IdentityMatrix[16]); (*Operador de difusi\[OAcute]n 1/2*)
Ax=Sum[Subscript[\[Phi]x, j].ket[j-1,4]\[ConjugateTranspose],{j,1,4}]; (*A \[Rule] B*)
By=Sum[Subscript[\[Psi]y, j].ket[j-1,4]\[ConjugateTranspose],{j,1,4}]; (*B \[Rule] A*)
H=HadamardMatrix[2]; (*Hadamard*)
s=KroneckerProduct[H,H,H,H].KroneckerProduct[ket[0,2],ket[0,2],ket[0,2],ket[0,2]]//N; (*Estado de superposici\[OAcute]n uniforme*)
(*W=(IdentityMatrix[16]-2 Ax.Ax\[ConjugateTranspose]).(2 s.s\[ConjugateTranspose]-IdentityMatrix[16]);*) (*Operador de difusi\[OAcute]n m\[AAcute]s Grover-like*)
(*W=(2 By.By\[ConjugateTranspose]-IdentityMatrix[16]).(2 Ax.Ax\[ConjugateTranspose]-IdentityMatrix[16]);*) (*Operador de difusi\[OAcute]n 2*)
W=MatrixPower[Ux,2]; (*Operador de difusi\[OAcute]n 1*)
"Operador de difusi\[OAcute]n cu\[AAcute]ntico" >>> output.txt
W >>> output.txt
(*W=MatrixPower[W/.\[Alpha]\[Rule]0.99,1/10];*)
(*\[Phi]x0=s;*)



a=0.85;
NSteps = 50;
qIvl1=Table[((\[Phi]x0\[ConjugateTranspose]).MatrixPower[(W\[ConjugateTranspose]),m].KroneckerProduct[IdentityMatrix[4],ket[0,4].ket[0,4]\[ConjugateTranspose]].MatrixPower[(W),m].(\[Phi]x0))[[1,1]],{m,0,NSteps}];
qIvl2=Table[((\[Phi]x0\[ConjugateTranspose]).MatrixPower[(W\[ConjugateTranspose]),m].KroneckerProduct[IdentityMatrix[4],ket[1,4].ket[1,4]\[ConjugateTranspose]].MatrixPower[(W),m].(\[Phi]x0))[[1,1]],{m,0,NSteps}];
qIvl3=Table[((\[Phi]x0\[ConjugateTranspose]).MatrixPower[(W\[ConjugateTranspose]),m].KroneckerProduct[IdentityMatrix[4],ket[2,4].ket[2,4]\[ConjugateTranspose]].MatrixPower[(W),m].(\[Phi]x0))[[1,1]],{m,0,NSteps}];
qIvl4=Table[((\[Phi]x0\[ConjugateTranspose]).MatrixPower[(W\[ConjugateTranspose]),m].KroneckerProduct[IdentityMatrix[4],ket[3,4].ket[3,4]\[ConjugateTranspose]].MatrixPower[(W),m].(\[Phi]x0))[[1,1]],{m,0,NSteps}];



Export["Nodo1-PRQ.png",ListLinePlot[qIvl1,PlotRange->{0,1},DataRange->{0,NSteps}]]
Export["Nodo2-PRQ.png",ListLinePlot[qIvl2,PlotRange->{0,1},DataRange->{0,NSteps}]]
Export["Nodo3-PRQ.png",ListLinePlot[qIvl3,PlotRange->{0,1},DataRange->{0,NSteps}]]
Export["Nodo4-PRQ.png",ListLinePlot[qIvl4,PlotRange->{0,1},DataRange->{0,NSteps}]]

Export["Nodo1-PRQ.csv",qIvl1]
Export["Nodo2-PRQ.csv",qIvl2]
Export["Nodo3-PRQ.csv",qIvl3]
Export["Nodo4-PRQ.csv",qIvl4]

"PageRank cu\[AAcute]ntico: " >>> "output.txt"
{Mean[qIvl1],Mean[qIvl2],Mean[qIvl3],Mean[qIvl4]} >>> "output.txt"











