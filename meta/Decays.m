(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["Decays`",
   {"SARAH`", "CConversion`", "CXXDiagrams`", "TreeMasses`", "TextFormatting`", "Utils`", "Vertices`"
      }];

FSParticleDecay::usage="head used for storing details of an particle decay,
in the format
   FSParticleDecay[particle, {final state particles}, {{loopOrder, {{topology, {diagrams}}, {topology, {diagrams}}}}, ...}]
where the topology is encoded as the relevant adjacency matrix.
";

IsSupportedDecayParticle::usage="returns True if decays for the given particle are supported.";

CreateCompleteParticleList::usage="";
GetDecaysForParticle::usage = "Creates 'objects' FSParticleDecay";
GetVerticesForDecays::usage="gets required vertices for a list of decays";

CreateSMParticleAliases::usage="creates aliases for SM particles present in model.";
CallDecaysCalculationFunctions::usage="creates calls to functions calculating
decays of the given particles.";
CreateDecaysCalculationPrototypes::usage="creates prototypes for convenience
functions calculating all decays for each particle.";
CreateDecaysCalculationFunctions::usage="creates definitions for convenience
functions calculating all decays for each particle.";
CreatePartialWidthCalculationPrototypes::usage="create prototypes for
functions computing partial widths of all decays.";
CreatePartialWidthCalculationFunctions::usage="creates definitions for
functions computing partial widths of all decays.";
CreateDecaysGetterFunctions::usage="create getters for specific particle decays";
CreateDecayTableGetterPrototypes::usage="create getter prototypes for C++ decay table.";
CreateDecayTableGetterFunctions::usage="create getter definitions for C++ decay table.";
CreateDecayTableInitialization::usage="create C++ initializer for decay table."
CreateTotalAmplitudeSpecializations::usage="creates specialized functions for higher-order decays.";
CreatePartialWidthSpecializations::usage="creates specialized functions for particular decays.";

WrapCodeInLoopOverInternalVertices::usage=""; (* @todo remove *)
Begin["`Private`"];

StripDiagramColorFactor[fields_, colorFactor_] :=
   Module[{fieldIn = fields[[1]], fieldsOut = Drop[fields, 1], initialStateRep, finalStateReps},
      initialStateRep = TreeMasses`GetColorRepresentation@fieldsIn;
      finalStateReps = TreeMasses`GetColorRepresentation /@ fieldsOut;
      Switch[initialStateRep,
         S, 
            Switch[finalStateReps,
               {S,S}, colorFactor,
               {T,-T}, Coefficient[colorFactor, Superscript[ColorMath`CM\[Delta], ColorMathInterface`GetFieldColorIndex/@fieldsOut]],
               {O,O}, Coefficient[colorFactor, Superscript[ColorMath`CM\[CapitalDelta], ColorMathInterface`GetFieldColorIndex/@fieldsOut]], 
               _, Print["Unhandled final state color representation in ", initialStateRep, " -> ", finalStateReps, " decay"]; Quit[1];
            ],
         _, Print["Unhandled initial state color representation ", initialStateRep, " -> ", finalStateReps, " decay"]; Quit[1];
      ]
   ];

GetInitialState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := particle;

GetFinalState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := finalState;

(*
  Returns a list of the form

     {{integer loop order, {{topology one, {diagram one, diagram two, ...}}, ...}},
      {integer loop order, {{topology one, {diagram one, diagram two, ...}}, ...}}}

  where each sublist contains the topologies and corresponding diagrams contributing
  to the decay at that particular loop order.
*)
GetDecayTopologiesAndDiagrams[FSParticleDecay[particle_, finalState_List, diagrams_List]] := diagrams;

(*
  Returns a list of the form

     {{integer loop order, {diagram one, diagram two, ...}},
      {integer loop order, {diagram one, diagram two, ...}}, ...}

  where each sublist contains the diagrams contributing to the decay at that
  particular loop-order.
*)
GetDecayDiagramsOnly[FSParticleDecay[particle_, finalState_List, diagrams_List]] :=
    {#[[1]], Flatten[Last /@ #[[2]], 1]}& /@ diagrams;

(*
   Returns a list of the form

      {{topology one, {diagram one, diagram two, ...}},
       {topology two, {diagram one, diagram two, ...}}, ...}

   containing all topologies and diagrams contributing to the decay
   at the given loop order.
*)
GetDecayTopologiesAndDiagramsAtLoopOrder[decay_FSParticleDecay, loopOrder_Integer] :=
    Module[{toposAndDiags},
           toposAndDiags = Select[GetDecayTopologiesAndDiagrams[decay], (First[#] == loopOrder)&];
           If[toposAndDiags =!= {},
              toposAndDiags = Last[Flatten[toposAndDiags, 1]];
             ];
           toposAndDiags
          ];

(*
  Returns a list of the form

     {diagram one, diagram two, ...}

  containing all diagrams contributing to the decay at the given loop order.
*)
GetDecayDiagramsOnlyAtLoopOrder[decay_FSParticleDecay, loopOrder_Integer] :=
    Module[{diags},
           diags = Select[GetDecayDiagramsOnly[decay], (First[#] == loopOrder)&];
           If[diags =!= {},
              diags = Last[Flatten[diags, 1]];
             ];
           diags
          ];

GetPossibleDecayTopologies[nProducts_, nLoops_] :=
    (
     Print["Error: decay topology with ", nProducts, " particles and ", nLoops, " loops not supported"];
     Quit[1]
    );

(* tree-level two-body decay, with
   vertex 1 = incoming state
   vertex 2 = outgoing state
   vertex 3 = outgoing state
   vertex 4 = internal vertex *)
GetPossibleDecayTopologies[2, 0] :=
    {
     {{0,0,0,1},
      {0,0,0,1},
      {0,0,0,1},
      {1,1,1,0}}
    };

GetPossibleDecayTopologies[2, 1] :=
    {
     {{0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1},
      {1,0,0,0,1,1},
      {0,1,0,1,0,1},
      {0,0,1,1,1,0}}
     ,
     {{0,0,0,1,0},
      {0,0,0,0,1},
      {0,0,0,0,1},
      {1,0,0,0,2},
      {0,1,1,2,0}}
     ,
     {{0,0,0,1,0},
      {0,0,0,1,0},
      {0,0,0,0,1},
      {1,1,0,0,2},
      {0,0,1,2,0}}
     ,
     {{0,0,0,1,0},
      {0,0,0,0,1},
      {0,0,0,1,0},
      {1,0,1,0,2},
      {0,1,0,2,0}}
    };

IsTreeLevelDecayTopology[t_] := MemberQ[GetPossibleDecayTopologies[2, 0], t];
IsTreeLevelTwoBodyDecayTopology[t_] := IsTreeLevelDecayTopology[t];

IsOneLoopDecayTopology[t_] := MemberQ[GetPossibleDecayTopologies[2, 1], t];
IsOneLoopTwoBodyDecayTopology[t_] := IsOneLoopDecayTopology[t];

GetPossibleDecayTopologies[nProducts_] := Join @@ (GetPossibleDecayTopologies[nProducts, #]& /@ {0, 1});

GetTreeLevelDecayTopologyName[nFinalParticles_Integer] := "e" <> ToString[nFinalParticles + 1] <> "_l0_t1";

GetOneLoopDecayTopologyName[nFinalParticles_Integer, topology_] :=
    Module[{pos},
           pos = Position[GetPossibleDecayTopologies[nFinalParticles, 1], topology];
           If[pos === {},
              Print["Error: unknown one-loop topology."];
              Quit[1];
             ];
           pos = First[First[pos]];
           "e" <> ToString[nFinalParticles + 1] <> "_l1_t" <> ToString[pos]
          ];

GetDecayTopologyName[t_] :=
    Which[MemberQ[GetPossibleDecayTopologies[2, 0], t],
          GetTreeLevelDecayTopologyName[2],
          MemberQ[GetPossibleDecayTopologies[2, 1], t],
          GetOneLoopDecayTopologyName[2, t],
          True,
          Print["Error: unknown decay topology provided"];
          Quit[1];
         ];

CreateCompleteParticleList[particles_List] := DeleteDuplicates[Join[particles, SARAH`AntiField[#]& /@ particles]];

GenericScalarName[] := "scalar";
GenericVectorName[] := "vector";
GenericFermionName[] := "fermion";
GenericGhostName[] := "ghost";

SimplifiedName[Susyno`LieGroups`conj[particle_]] :=
    Susyno`LieGroups`conj[SimplifiedName[particle]];
SimplifiedName[SARAH`bar[particle_]] :=
    SARAH`bar[SimplifiedName[particle]];

SimplifiedName[particle_?TreeMasses`IsSMLepton] := "lep";
SimplifiedName[particle_?TreeMasses`IsSMDownQuark] := "dq";
SimplifiedName[particle_?TreeMasses`IsSMUpQuark] := "uq";
SimplifiedName[particle_ /; TreeMasses`GetHiggsBoson[] =!= Null && particle === TreeMasses`GetHiggsBoson[]] := "H";
SimplifiedName[particle_ /; TreeMasses`GetPseudoscalarHiggsBoson[] =!= Null && particle === TreeMasses`GetPseudoscalarHiggsBoson[]] := "AH";
SimplifiedName[particle_ /; TreeMasses`GetWBoson[] =!= Null && particle === TreeMasses`GetWBoson[]] := "W";
SimplifiedName[particle_ /; TreeMasses`GetZBoson[] =!= Null && particle === TreeMasses`GetZBoson[]] := "Z";
SimplifiedName[particle_ /; TreeMasses`GetPhoton[] =!= Null && particle === TreeMasses`GetPhoton[]] := "A";
SimplifiedName[particle_ /; TreeMasses`GetGluon[] =!= Null && particle === TreeMasses`GetGluon[]] := "G";
SimplifiedName[particle_] := particle;

CreateParticleAlias[particle_, namespace_String] :=
    "using " <> SimplifiedName[particle] <> " = " <>
    CXXDiagrams`CXXNameOfField[particle, prefixNamespace -> namespace] <> ";";

CreateParticleAliases[particles_, namespace_:""] :=
    Utils`StringJoinWithSeparator[CreateParticleAlias[#, namespace]& /@ particles, "\n"];

CreateSMParticleAliases[namespace_:""] :=
    Module[{smParticlesToAlias},
           smParticlesToAlias = Select[{TreeMasses`GetHiggsBoson[],
                                        TreeMasses`GetPseudoscalarHiggsBoson[],
                                        TreeMasses`GetWBoson[], TreeMasses`GetZBoson[],
                                        TreeMasses`GetGluon[], TreeMasses`GetPhoton[],
                                        TreeMasses`GetDownLepton[1] /. field_[generation_] :> field,
                                        TreeMasses`GetUpQuark[1] /. field_[generation_] :> field,
                                        TreeMasses`GetDownQuark[1] /.field_[generation_] :> field
                                       }, (# =!= Null)&];
           CreateParticleAliases[smParticlesToAlias, namespace]
          ];

GetGenericTypeName[p_?TreeMasses`IsScalar] := GenericScalarName[];
GetGenericTypeName[p_?TreeMasses`IsVector] := GenericVectorName[];
GetGenericTypeName[p_?TreeMasses`IsFermion] := GenericFermionName[];
GetGenericTypeName[p_?TreeMasses`IsGhost] := GenericGhostName[];

GetFeynArtsTypeName[p_?TreeMasses`IsScalar] := S;
GetFeynArtsTypeName[p_?TreeMasses`IsVector] := V;
GetFeynArtsTypeName[p_?TreeMasses`IsFermion] := F;
GetFeynArtsTypeName[p_?TreeMasses`IsGhost] := U;

IsSupportedDecayParticle[particle_] :=
    (
     !TreeMasses`IsGhost[particle] &&
     !TreeMasses`IsVector[particle] &&
     !TreeMasses`IsMassless[particle] &&
     TreeMasses`GetDimensionWithoutGoldstones[particle] > 0
    );

(* returns False if state consists only of Goldstones or ghosts, for any set of generation indices *)
IsPhysicalFinalState[finalState_List] :=
    Module[{goldstones, onlyGoldstones},
           goldstones = Select[finalState, TreeMasses`IsGoldstone];
           isAlwaysGoldstone = (TreeMasses`GetDimensionWithoutGoldstones[#] == 0)& /@ goldstones;
           (Select[finalState, TreeMasses`IsGhost] === {}) && !(Or @@ isAlwaysGoldstone)
          ];

IsElectricChargeConservingDecay[initialParticle_, finalState_List] :=
    Module[{chargeSum},
           chargeSum = Simplify[Plus @@ (Join[{-TreeMasses`GetElectricCharge[initialParticle]}, TreeMasses`GetElectricCharge /@ finalState])];
           PossibleZeroQ[chargeSum]
          ];

(* @todo handle more than 2 particles in final state and non-SM color representations *)
IsColorInvariantDecay[initialParticle_, finalState_List] :=
    Module[{initialStateRep, finalStateReps, result = True},
           If[Length[finalState] == 2,
              initialStateRep = TreeMasses`GetColorRepresentation[initialParticle];
              finalStateReps = Sort[TreeMasses`GetColorRepresentation /@ finalState];
              Switch[initialStateRep,
                     S, result = ((finalStateReps === {S, S}) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}]) ||
                                  (finalStateReps === {O, O}));,
                     T|-T, result = ((finalStateReps === Sort[{T, S}]) ||
                                     (finalStateReps === Sort[{-T, S}]) ||
                                     (finalStateReps === Sort[{O, S}]));,
                     O, result = ((finalStateReps === Sort[{O, S}]) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}]));,
                     _, result = True; (* unhandled case *)
                    ];
             ];
           result
          ];

FinalStateContainsInitialState[initialParticle_, finalState_List] :=
    Module[{containsInitialMultiplet, dim},
           containsInitialMultiplet = !FreeQ[finalState, initialParticle];
           dim = TreeMasses`GetDimension[initialParticle];
           containsInitialMultiplet && dim == 1
          ];

IsPossibleNonZeroVertex[fields_List, useDependences_:False] :=
    Module[{numFields, cachedVertices = {}},
           numFields = Length[fields];
           Vertices`IsNonZeroVertex[fields, Vertices`GetCachedVertices[numFields, useDependences], useDependences]
         ];

IsPossibleNonZeroDiagram[diagram_, useDependences_:False] :=
    Module[{vertices},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           And @@ (IsPossibleNonZeroVertex[#, useDependences]& /@ vertices)
          ];

IsPossibleTreeLevelDecay[decay_FSParticleDecay, useDependences_:False] :=
    Module[{treeLevelDiags = GetDecayDiagramsOnlyAtLoopOrder[decay, 0]},
           treeLevelDiags =!= {} && (And @@ (IsPossibleNonZeroDiagram[#, useDependences]& /@ treeLevelDiags))
          ];

IsPossibleOneLoopDecay[decay_FSParticleDecay] :=
    GetDecayDiagramsOnlyAtLoopOrder[decay, 1] =!= {};

ContainsOnlySupportedVertices[diagram_] :=
    Module[{vertices, vertexTypes, unsupportedVertices},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           vertexTypes = CXXDiagrams`VertexTypeForFields /@ vertices;
           unsupportedVertices = Complement[vertexTypes, CXXDiagrams`VertexTypes[]];
           If[unsupportedVertices =!= {},
              MapIndexed[(If[!MemberQ[CXXDiagrams`VertexTypes[], vertexTypes[[First[#2]]]],
                             Print["Warning: vertex with fields ", #1, " is not currently supported."];
                             Print["    Diagrams involving this vertex will be discarded."];
                            ];)&, vertices];
             ];
           unsupportedVertices === {}
          ];

IsSupportedDiagram[diagram_] := ContainsOnlySupportedVertices[diagram];

GetFinalStateExternalField[particle_] := SARAH`AntiField[particle];

GetContributingDiagramsForDecayGraph[initialField_, finalFields_List, graph_] :=
    Module[{externalFields, diagrams},
           externalFields = Join[{1 -> initialField}, MapIndexed[(First[#2] + 1 -> #1)&, finalFields]];
           diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph, externalFields];
           Select[diagrams, IsPossibleNonZeroDiagram]
          ];

GetContributingGraphsForDecay[initialParticle_, finalParticles_List, maxLoops_Integer] :=
    Module[{i, nFinalParticles = Length[finalParticles], topologies, diagrams},
           topologies = Join[Table[{i, GetPossibleDecayTopologies[nFinalParticles, i]}, {i, 0, maxLoops}]];
           diagrams = {#[[1]], {#, GetContributingDiagramsForDecayGraph[initialParticle, GetFinalStateExternalField /@ finalParticles, #]}& /@ #[[2]]}&
                      /@ topologies;
           diagrams = {#[[1]], With[{toposAndDiags = #[[2]]}, Select[toposAndDiags, #[[2]] =!= {}&]]}& /@ diagrams;
           diagrams = DeleteCases[diagrams, {_Integer, {}}];
           {#[[1]], With[{toposAndDiags = #[[2]]}, {#[[1]], Select[#[[2]], IsSupportedDiagram]}& /@ toposAndDiags]}& /@ diagrams
          ];

GetContributingGraphsForDecay[initialParticle_, finalParticles_List] :=
    GetContributingGraphsForDecay[initialParticle, finalParticles, 1];

(* defines a fixed ordering for final state particles  *)
(* @todo decide on what this ordering actually will be *)
OrderFinalState[initialParticle_?TreeMasses`IsScalar, finalParticles_List] :=
    Module[{orderedFinalState},
           orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
           orderedFinalState = Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]];
           If[Length[orderedFinalState] === 2,
              (* re-order to SSV *)
              If[TreeMasses`IsVector[orderedFinalState[[1]]] && TreeMasses`IsScalar[orderedFinalState[[2]]],
                 orderedFinalState = Reverse[orderedFinalState]
                ];
              (* re-order to S bar[F] F *)
              If[TreeMasses`IsFermion[orderedFinalState[[1]]] && TreeMasses`IsFermion[orderedFinalState[[2]]],
                 If[Head[orderedFinalState[[2]]] === SARAH`bar && !Head[orderedFinalState[[1]]] === SARAH`bar,
                    orderedFinalState = Reverse[orderedFinalState];
                   ];
                ];
             ];
            orderedFinalState
          ];

OrderFinalState[initialParticle_, finalParticles_List] :=
    Module[{orderedFinalState},
           orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
           Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]]
          ];

GetDecaysForParticle[particle_, Infinity, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {n_, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity, n_}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2,
                     allowedFinalStateParticles_List] :=
    GetDecaysForParticle[particle, {2, maxNumberOfProducts}, allowedFinalStateParticles];

GetDecaysForParticle[particle_, {minNumberOfProducts_Integer /; minNumberOfProducts >= 2,
                                 maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2},
                                allowedFinalStateParticles_List] :=
    Module[{i, finalStateSizes},
           finalStateSizes = Table[{i}, {i, minNumberOfProducts, maxNumberOfProducts}];
           Flatten[GetDecaysForParticle[particle, #, allowedFinalStateParticles]& /@ finalStateSizes, 1]
          ];

GetDecaysForParticle[particle_, {exactNumberOfProducts_Integer}, allowedFinalStateParticles_List] :=
    Module[{genericFinalStates, finalStateParticlesClassified,
            isPossibleDecay, concreteFinalStates, decays},

           If[exactNumberOfProducts > 2,
              Print["Error: decays with ", exactNumberOfProducts,
                    " final particles are not currently supported."];
              Quit[1];
           ];

           genericFinalStates = GetAllowedGenericFinalStates[particle, exactNumberOfProducts];

           (* @todo checks on colour and Lorentz structure *)
           isPossibleDecay[finalState_] := (IsPhysicalFinalState[finalState] &&
                                            IsElectricChargeConservingDecay[particle, finalState] &&
                                            IsColorInvariantDecay[particle, finalState] &&
                                            !FinalStateContainsInitialState[particle, finalState]);
           concreteFinalStates = Join @@ (GetParticleCombinationsOfType[#, allowedFinalStateParticles, isPossibleDecay]& /@ genericFinalStates);
           concreteFinalStates = OrderFinalState[particle, #] & /@ concreteFinalStates;
           decays = FSParticleDecay[particle, #, GetContributingGraphsForDecay[particle, #]]& /@ concreteFinalStates;
           Select[decays, GetDecayTopologiesAndDiagrams[#] =!= {}&]
          ];

GetDecaysForParticle[particle_, n_, allowedFinalStateParticles_List] :=
    (
     Print["Error: invalid number of final state particles: ", n];
     Quit[1];
    )

GatherParticlesByType[particles_List] :=
    Module[{areSameType},
           areSameType[p1_, p2_] := Or @@ ((#[p1] && #[p2])& /@ { TreeMasses`IsScalar,
                                                                  TreeMasses`IsVector,
                                                                  TreeMasses`IsFermion,
                                                                  TreeMasses`IsGhost });
           Gather[particles, areSameType]
          ];

(* returns a list of lists of the form {{particle type 1, {particles of type 1}}, {particle type 2, {particles of type 2}}, ...} *)
GetClassifiedParticleLists[particles_List] :=
    Module[{classified, foundTypes},
           classified = {GetGenericTypeName[First[#]], #}& /@ GatherParticlesByType[particles];
           foundTypes = First /@ classified;
           If[Length[Union[foundTypes]] != Length[foundTypes],
              Print["Error: particles incorrectly classified: ", classified];
              Quit[1];
             ];
           classified
          ];

BaseMulticombination[k_] := Module[{i}, Table[1, {i, 1, k}]];

(* see, e.g., http://www.martinbroadhurst.com/multicombinations.html *)
NextMulticombination[n_, combination_] :=
    Module[{k = Length[combination], i, incrementable, pos, val},
           incrementable = Position[combination, x_ /; x < n];
           If[Length[incrementable] == 0,
              {},
              pos = First[Last[incrementable]];
              val = combination[[pos]] + 1;
              Join[Take[combination, pos - 1], {val}, Table[val, {i, 1, k - pos}]]
             ]
          ];

NextMulticombinationsList[setSizes_List, combinations_List] :=
    Module[{numCombinations = Length[combinations], next},
           next = combinations;
           For[i = numCombinations, i > 0, i--,
               nextCombination = NextMulticombination[setSizes[[i]], combinations[[i]]];
               If[nextCombination =!= {},
                  next[[i]] = nextCombination;
                  Return[next];,
                  next[[i]] = BaseMulticombination[Length[next[[i]]]];
                 ];
              ];
           {}
          ];

GetParticleCombinationsOfType[genericState_List, particles_List, isValidTest_:Function[True]] :=
    Module[{genericTypeCounts, classifiedParticles, indexLists, candidate, combinations},
           genericTypeCounts = {#, Count[genericState, #]}& /@ DeleteDuplicates[genericState];
           classifiedParticles = Select[GetClassifiedParticleLists[particles], MemberQ[genericState, First[#]]&];
           genericTypeCounts = genericTypeCounts[[Flatten[Position[First /@ classifiedParticles, First[#]]& /@ genericTypeCounts]]];
           indexLists = BaseMulticombination[Last[#]]& /@ genericTypeCounts;
           combinations = Reap[
               While[indexLists =!= {},
                     candidate = Flatten[MapIndexed[With[{pos = First[#2], indices = #1},
                                                    Last[classifiedParticles[[pos]]][[#]]& /@ indices]&, indexLists]];
                     If[isValidTest[candidate],
                        Sow[candidate];
                       ];
                     indexLists = NextMulticombinationsList[Length[Last[#]]& /@ classifiedParticles, indexLists];
                    ];
               ];
           Flatten[Last[combinations], 1]
          ];

GetAllowedGenericFinalStates[particle_ /; (TreeMasses`IsScalar[particle] || TreeMasses`IsVector[particle]),
                             n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericScalarName[]},
               {GenericScalarName[], GenericVectorName[]},
               {GenericVectorName[], GenericVectorName[]},
               {GenericFermionName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetAllowedGenericFinalStates[particle_?TreeMasses`IsFermion, n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericFermionName[]},
               {GenericVectorName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetVerticesForDecay[decay_FSParticleDecay] :=
    Module[{diagrams = GetDecayDiagramsOnly[decay]},
           DeleteDuplicates[Join @@ (Flatten[(CXXDiagrams`VerticesForDiagram /@ Last[#]), 1]& /@ diagrams)]
          ];

GetVerticesForDecays[particleDecays_List] :=
    DeleteDuplicates[Flatten[GetVerticesForDecay /@ particleDecays, 1]]

LoopOverIndex[loopBody_String, index_, start_, stop_, type_:CConversion`ScalarType[CConversion`integerScalarCType]] :=
    Module[{idxStr, startStr, stopStr},
           idxStr = CConversion`ToValidCSymbolString[index];
           startStr = CConversion`ToValidCSymbolString[start];
           stopStr = CConversion`ToValidCSymbolString[stop];
           "for (" <> CConversion`CreateCType[type] <> " " <> idxStr <>
           " = " <> startStr <> "; " <> idxStr <> " < " <> stopStr <> "; ++" <>
           idxStr <> ") {\n" <> TextFormatting`IndentText[loopBody] <> "\n}"
          ];

(* generates a loop over the given indices, in the form
   {{idx1, start, stop}, {idx2, start, stop}, ...} with
   the first list entry being the innermost loop *)
LoopOverIndexCollection[loopBody_String, indices_List] :=
    Fold[LoopOverIndex[#1, Sequence @@ #2]&, loopBody, indices];

CreateGenericGetPartialWidthFunctionName[] := "get_partial_width";

CreateSpecializedPartialWidthCalculationName[initialState_, finalState_List, fieldsNamespace_] :=
    CreateGenericGetPartialWidthFunctionName[] <> "<" <>
    CXXDiagrams`CXXNameOfField[initialState, prefixNamespace -> fieldsNamespace] <> "," <>
    Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#, prefixNamespace -> fieldsNamespace]& /@ finalState, ","] <> " >";

CreatePartialWidthCalculationName[decay_FSParticleDecay, scope_:""] :=
    Module[{initialState, initialStateName,
            finalState, finalStateName},
           initialState = GetInitialState[decay];
           initialStateName = CConversion`ToValidCSymbolString[initialState];
           finalState = GetFinalState[decay];
           finalStateName = StringJoin[CConversion`ToValidCSymbolString /@ finalState];
           scope <> If[scope != "", "::", ""] <> "partial_width_" <> initialStateName <> "_to_" <> finalStateName
          ];

CreatePartialWidthCalculationPrototype[decay_FSParticleDecay] :=
    Module[{returnType = "", functionName = "", functionArgs = "",
            initialStateDim, finalStateDims},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = FlexibleSUSY`FSModelName <> "_mass_eigenstates&" <>
                          If[initialStateDim > 1, ", int", ""] <>
                          StringJoin[If[# > 1, ", int", ""]& /@ finalStateDims];
           functionName = CreatePartialWidthCalculationName[decay];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const;"
          ];

CreatePartialWidthCalculationFunction[decay_FSParticleDecay, fieldsNamespace_] :=
    Module[{i, returnType = "", functionName = "", functionArgs = "",
            initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, setFieldIndices, body = ""},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = FlexibleSUSY`FSModelName <> "_mass_eigenstates& model" <>
                          If[initialStateDim > 1, ", int gI1", ""] <>
                          StringJoin[MapIndexed[If[#1 > 1, ", int gO" <> ToString[First[#2]], ""]&, finalStateDims]];
           functionName = CreatePartialWidthCalculationName[decay, "CLASSNAME"];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           setFieldIndices[field_, indicesName_, indexVal_] :=
               Module[{i, dim, numIndices, result = ""},
                      dim = TreeMasses`GetDimension[field];
                      numIndices = CXXDiagrams`NumberOfFieldIndices[field];
                      result = "const typename field_indices<" <> CXXDiagrams`CXXNameOfField[field] <> " >::type " <> indicesName;
                      If[numIndices == 0 || dim <= 1,
                         result = result <> "{};\n";,
                         result = result <> "{{" <> ToString[indexVal] <>
                                  StringJoin[Table[", 0", {i, 1, numIndices - 1}]] <> "}};\n";
                        ];
                      result
                     ];
           body = "context_base context{model};\n" <>
                  (* @todo ask Jobst which usings are necessary *)
                  "using namespace " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams;\n" <>
                  "using namespace " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields;\n" <>
                  "using fields::bar;" <>
                  "using fields::conj;" <>
                  StringJoin[setFieldIndices[#[[1]], #[[2]], #[[3]]]& /@
                                 Join[{{initialState, "in_indices", If[initialStateDim > 1, "gI1", ""]}},
                                      MapIndexed[{#1, "out_" <> ToString[First[#2]] <> "_indices",
                                                  If[finalStateDims[[First[#2]]] > 1, "gO" <> ToString[First[#2]], ""]}&, finalState]]];
           body = body <> "\nreturn " <> CreateSpecializedPartialWidthCalculationName[initialState, finalState, fieldsNamespace] <>
                  "(context, in_indices" <> StringJoin[Table[", out_" <> ToString[i] <> "_indices", {i, 1, Length[finalState]}]] <> ");\n";
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const\n{\n" <>
                  TextFormatting`IndentText[body] <> "}\n"
          ];

CreatePartialWidthCalculationPrototypes[particleDecays_List] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationPrototype /@ allDecays, "\n"]
          ];

CreatePartialWidthCalculationFunctions[particleDecays_List, fieldsNamespace_] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationFunction[#, fieldsNamespace]& /@ allDecays, "\n"]
          ];

CallPDGCodeGetter[SARAH`bar[particle_], args__] :=
    "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[Susyno`LieGroups`conj[particle_], args__] :=
   "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[particle_, idx_String, namespace_] :=
    Module[{dim = TreeMasses`GetDimension[particle], particleStr, result = ""},
           particleStr = namespace <> If[namespace != "", "::", ""] <> CConversion`ToValidCSymbolString[particle];
           result = namespace <> If[namespace != "", "::", ""] <> "get_pdg_code_for_particle(" <>
                    particleStr;
           If[dim > 1,
              result = result <> ", " <> idx;
             ];
           result <> ")"
          ];

CallPartialWidthCalculation[decay_FSParticleDecay] :=
    Module[{i, initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, functionArgs = "",
            pdgsList = "", loopIndices, body = ""},
           initialStateDim = TreeMasses`GetDimension[initialState];
           finalStateDims = TreeMasses`GetDimension /@ finalState;
           finalStateStarts = TreeMasses`GetDimensionStartSkippingGoldstones /@ finalState;
           functionArgs = "model_" <> If[initialStateDim > 1, ", gI1", ""] <>
                          MapIndexed[If[#1 > 1, ", gO" <> ToString[First[#2]], ""]&, finalStateDims];
           pdgsList = MapIndexed[With[{idx = First[#2]},
                                      CallPDGCodeGetter[#1, If[finalStateDims[[idx]] > 1, "gO" <> ToString[idx], ""], FlexibleSUSY`FSModelName <> "_info"]]&,
                                 finalState];
           pdgsList = "{" <> Utils`StringJoinWithSeparator[pdgsList, ", "] <> "}";
           body = "decays.set_decay(" <> CreatePartialWidthCalculationName[decay] <> "(" <> functionArgs <> "), " <> pdgsList <> ");";
           loopIndices = Reverse[Select[MapIndexed[With[{idx = First[#2]},
                                                        If[#1 > 1,
                                                           {"gO" <> ToString[idx], Utils`MathIndexToCPP@finalStateStarts[[idx]], #1},
                                                           {}
                                                          ]
                                                       ]&, finalStateDims], (# =!= {})&]];
           If[loopIndices =!= {},
              body = LoopOverIndexCollection[body, loopIndices];
             ];
           body <> "\n"
          ];

CreateDecaysCalculationFunctionName[particle_, scope_:""] :=
    scope <> If[scope =!= "", "::", ""] <> "calculate_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysCalculationPrototype[particle_] :=
    "void " <> CreateDecaysCalculationFunctionName[particle] <> "();";

CreateDecaysCalculationPrototypes[decayLists_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationPrototype[First[#]]& /@ decayLists, "\n"];

CreateLocalScope[body_] := "{\n" <> TextFormatting`IndentText[body] <> "}\n";

CreateDecaysCalculationFunction[decaysList_] :=
    Module[{particle = First[decaysList], particleDim, particleStart,
            decayChannels = Last[decaysList],
            runToScale = "", body = ""},
           particleDim = TreeMasses`GetDimension[particle];
           particleStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
           runToScale = "const auto& decay_mass = PHYSICAL(" <>
                        CConversion`ToValidCSymbolString[TreeMasses`GetMass[particle]] <>
                        ");\nmodel_.run_to(decay_mass" <> If[particleDim > 1, "(gI1)", ""] <> ");\n" <>
                        "model_.calculate_DRbar_masses();\n";
           body = StringJoin[CallPartialWidthCalculation /@ decayChannels];
           body = "auto& decays = decay_table.get_" <> CConversion`ToValidCSymbolString[particle] <>
                  "_decays(" <> If[particleDim > 1, "gI1", ""] <> ");\n\n" <> body;
           body = "auto model_ = model;\n\nif (run_to_decay_particle_scale) {\n" <>
                  TextFormatting`IndentText[runToScale] <> "}\n\n" <> body;
           If[particleDim > 1,
              body = LoopOverIndexCollection[body, {{"gI1", Utils`MathIndexToCPP@particleStart, particleDim}}] <> "\n";
             ];
           "void " <> CreateDecaysCalculationFunctionName[particle, "CLASSNAME"] <>
           "()\n{\n"
           <> TextFormatting`IndentText[body] <> "}\n"
          ];

CreateDecaysCalculationFunctions[particleDecays_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationFunction /@ particleDecays, "\n"];

CallDecaysFunction[particle_, arg_:"model", obj_:""] :=
    obj <> CreateDecaysCalculationFunctionName[particle] <> "();\n"

CallThreadedDecaysFunction[particle_, ptr_:"this", pool_:"tp"] :=
    pool <> ".run_task([" <> ptr <> "] () { " <>
    If[ptr === "this", "", ptr <> "->"] <>
    CreateDecaysCalculationFunctionName[particle] <> "(); });\n";

CallDecaysCalculationFunctions[particles_List, enableDecaysThreads_] :=
    Module[{result = ""},
           If[enableDecaysThreads,
              result = "Thread_pool tp(std::min(std::thread::hardware_concurrency(), " <>
                       ToString[Length[particles]] <> "u));\n\n" <>
                       StringJoin[CallThreadedDecaysFunction /@ particles];
              ,
              result = StringJoin[CallDecaysFunction /@ particles];
             ];
           result
          ];

GetDecayAmplitudeType[initialParticle_?TreeMasses`IsScalar, finalState_List] :=
    Module[{vertexType},
           vertexType = CXXDiagrams`VertexTypeForFields[Join[{initialParticle}, finalState]];
           Switch[vertexType,
                  ScalarVertex, "Decay_amplitude_SSS",
                  ChiralVertex, "Decay_amplitude_SFF",
                  MomentumDifferenceVertex, "Decay_amplitude_SSV",
                  InverseMetricVertex, "Decay_amplitude_SVV",
                  _, Print["Warning: decay ", initialParticle, " -> ", finalState, " is not supported."];
                     "Unknown_amplitude_type"
                 ]
          ];

GetDecayAmplitudeType[initialParticle_?TreeMasses`IsFermion, finalState_List] :=
    Module[{vertexType},
           vertexType = CXXDiagrams`VertexTypeForFields[Join[{initialParticle}, finalState]];
           Switch[vertexType,
                  ChiralVertex, "Decay_amplitude_FFS",
                  ChiralVertex, "Decay_amplitude_FFV",
                  _, Print["Warning: decay ", initialParticle, " -> ", finalState, " is not supported."];
                     "Unknown_amplitude_type"
                 ]
          ];

GetDecayAmplitudeType[decay_FSParticleDecay] :=
    GetDecayAmplitudeType[GetInitialState[decay], GetFinalState[decay]];

CreateFieldIndices[particle_String] :=
    "typename " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::field_indices<" <> particle <> " >::type";

CreateFieldIndices[particle_, fieldsNamespace_] :=
    CreateFieldIndices[CXXDiagrams`CXXNameOfField[particle, prefixNamespace -> fieldsNamespace]];

CreateTotalAmplitudeFunctionName[] := "calculate_amplitude";

CreateTotalAmplitudeSpecializationDecl[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            returnType = "", fieldsNamespace, fieldsList, templatePars = "", args = ""},
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           fieldsNamespace = FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields";
           fieldsList = Join[{initialParticle}, finalState];
           templatePars = "<" <> Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#, prefixNamespace -> modelName <> "_cxx_diagrams::fields"]& /@
                                                               fieldsList, ", "] <> ">";
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[#, fieldsNamespace] <> "&")& /@ fieldsList, ", "];
           "template<>\n" <> returnType <> " " <> modelName <> "_decays::" <>
           CreateTotalAmplitudeFunctionName[] <> templatePars <> "(" <> args <> ") const;\n"
          ];

FillSSSDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_out_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_2);\n";
           assignments = assignments <> structName <> ".m_out_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_3);\n";
           assignments
          ];

FillSFFDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_out_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_2);\n";
           assignments = assignments <> structName <> ".m_out_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_3);\n";
           assignments
          ];

FillSSVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, finalState, scalar, scalarPos, vector, vectorPos, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           finalState = GetFinalState[decay];
           scalar = First[Select[finalState, TreeMasses`IsScalar]];
           scalarPos = First[First[Position[finalState, scalar]]];
           vector = First[Select[finalState, TreeMasses`IsVector]];
           vectorPos = First[First[Position[finalState, vector]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_scalar = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[scalar, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[scalarPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_vector = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[vector, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[vectorPos + 1] <> ");\n";
           assignments
          ];

FillSVVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_out_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_2);\n";
           assignments = assignments <> structName <> ".m_out_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_3);\n";
           assignments
          ];

FillFFSDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, finalState, scalar, scalarPos, fermion, fermionPos, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           finalState = GetFinalState[decay];
           fermion = First[Select[finalState, TreeMasses`IsFermion]];
           fermionPos = First[First[Position[finalState, fermion]]];
           scalar = First[Select[finalState, TreeMasses`IsScalar]];
           scalarPos = First[First[Position[finalState, scalar]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_fermion = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[fermion, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[fermionPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_scalar = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[scalar, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[scalarPos + 1] <> ");\n";
           assignments
          ];

FillFFVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsNamespace, finalState, fermion, fermionPos, vector, vectorPos, assignments = ""},
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           finalState = GetFinalState[decay];
           fermion = First[Select[finalState, TreeMasses`IsFermion]];
           fermionPos = First[First[Position[finalState, fermion]]];
           vector = First[Select[finalState, TreeMasses`IsVector]];
           vectorPos = First[First[Position[finalState, vector]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay], prefixNamespace -> fieldsNamespace] <>
                         " >(idx_1);\n";
           assignments = assignments <> structName <> ".m_fermion = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[fermion, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[fermionPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_vector = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[vector, prefixNamespace -> fieldsNamespace] <>
                         " >(idx_" <> ToString[vectorPos + 1] <> ");\n";
           assignments
          ];

FillDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Switch[GetDecayAmplitudeType[decay],
           "Decay_amplitude_SSS", FillSSSDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SFF", FillSFFDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SSV", FillSSVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SVV", FillSVVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFS", FillFFSDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFV", FillFFVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           _, ""
          ];

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SSS",
                              structName_] :=
    structName <> ".form_factor = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SSV",
                              structName_] :=
    structName <> ".form_factor = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SVV",
                              structName_] :=
    structName <> ".form_factor_g = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_11 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_12 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_21 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_22 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_eps = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SFF",
                              structName_] :=
    structName <> ".form_factor_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_right = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_FFS",
                              structName_] :=
    structName <> ".form_factor_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_right = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_FFV",
                              structName_] :=
    structName <> ".form_factor_gam_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_gam_right = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_p_1 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_p_2 = std::complex<double>(0., 0.);\n";

GetTreeLevelTwoBodyDecayVertex[decay_FSParticleDecay] :=
    Module[{treeLevelDiags, vertices = {}},
           treeLevelDiags = GetDecayDiagramsOnlyAtLoopOrder[decay, 0];
           If[treeLevelDiags =!= {},
              vertices = Flatten[CXXDiagrams`VerticesForDiagram /@ treeLevelDiags, 1];
             ];
           vertices
          ];

EvaluateTreeLevelTwoBodyDecayVertex[decay_FSParticleDecay, modelName_, indicesName_, paramsStruct_, resultName_:"vertex"] :=
    Module[{vertexFields, fieldsNamespace, templatePars},
           vertexFields = GetTreeLevelTwoBodyDecayVertex[decay];
           If[Length[vertexFields] > 1,
              Print["Error: more than a single vertex in tree-level decays."];
              Quit[1];
             ];
           If[vertexFields =!= {},
              vertexFields = First[vertexFields];
              fieldsNamespace = modelName <> "_cxx_diagrams::fields";
              templatePars = "<" <>
                              Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#, prefixNamespace -> fieldsNamespace]&
                                                            /@ vertexFields, ", "] <> " >";
              "const auto " <> resultName <> " = " <> modelName <> "_cxx_diagrams::Vertex" <> templatePars <> "::evaluate(" <>
              indicesName <> ", " <> paramsStruct <> ");\n",
              ""
             ]
          ];

FillSSSTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor += vertex.value();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSFFTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_left += vertex.left();\n" <>
                         structName <> ".form_factor_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSSVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor += vertex.value(0, 1);\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSVVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_g += vertex.value();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillFFSTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_left += vertex.left();\n" <>
                         structName <> ".form_factor_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillFFVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{i, fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_gam_left += vertex.left();\n" <>
                         structName <> ".form_factor_gam_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Switch[GetDecayAmplitudeType[decay],
           "Decay_amplitude_SSS", FillSSSTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SFF", FillSFFTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SSV", FillSSVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SVV", FillSVVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFS", FillFFSTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFV", FillFFVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           _, ""
          ];
Compare[a_, b_] := Module[{},

   For[i = 1, i <= Length[a], i++,
      If[Select[b, (a[[i,2]] === #[[2]] && Sort[a[[i,1]]] === Sort[#[[1]]])&] === {},
         Return[False]
      ]
   ];

   Return[True];
];

TranslationForDiagram[topology_, diagram_] := Module[{
   diagramsWithCorrectTopology,
   file = Get["utils/loop_decays/output/generic_loop_decay_diagram_classes.m"], temp,res
},

   diagramsWithCorrectTopology = Select[file, MemberQ[#, topology]&];
   Utils`AssertWithMessage[diagramsWithCorrectTopology =!= {},
      "Can't find topology " <> ToString@topology <> " for diagram " <> ToString@diagram
   ];
   temp = (#[[1]] -> GetFeynArtsTypeName[#[[2]]])& /@ InsertionsOnEdgesForDiagram[topology, diagram];

   res = Select[diagramsWithCorrectTopology, Compare[#[[3]]/.#[[4]], temp]&];

   Utils`AssertWithMessage[Length[res] === 1,
      "Error! Couldn't find translation for a diagram"
   ];

   First@res
];

EvaluateOneLoopTwoBodyDecayDiagramWithTopology[decay_, topology_, diagram_] :=
    Module[{
       file = Get["utils/loop_decays/output/generic_loop_decay_diagram_classes.m"],
       diagramsWithCorrectTopology, cos3, particles, particles2, con1, con2, con3, p1, p2, p3,
      numberOfVertices = Count[diagram, el_ /; Head[el] === List],
      externalParticles = Join[{GetInitialState[decay]}, GetFinalState[decay]],
      internalParticles, permutedParticles, vertices, temp},

      temp = (#[[1]] -> GetFeynArtsTypeName[#[[2]]])& /@ InsertionsOnEdgesForDiagram[topology, diagram];
      Switch[numberOfVertices,
         (* half-candy diagram *)
         2, 
            con1 = CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[4, 5, diagram, topology];
            p1 = diagram[[4, con1[[1,1]]]];
            p2 = diagram[[4, con1[[2,1]]]];,
         3, 
            con1 = CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[4, 5, diagram, topology];
            p1 = diagram[[4, con1[[1,1]]]];
            con2 = CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[4, 6, diagram, topology];
            p2 = diagram[[4, con2[[1,1]]]];
            con3 = CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[5, 6, diagram, topology];
            p3 = diagram[[5, con3[[1,1]]]];
      ];

      diagramsWithCorrectTopology = Select[file, MemberQ[#, topology]&];
      Utils`AssertWithMessage[diagramsWithCorrectTopology =!= {},
         "Can't find topology " <> ToString@topology <> " for diagram " <> ToString@diagram
      ];

      internalParticles = If[numberOfVertices === 2, {p1,p2}, {p1, p2, p3}];
      particles = Join[externalParticles, internalParticles];

      If[GetFeynArtsTypeName /@ particles === {S,V,V,V,S},
        particles = Join[Drop[particles, -2], {particles[[5]], particles[[4]]}];
      ];
      particles2 = MapIndexed[Field[#2[[1]]] -> #1&, GetFeynArtsTypeName /@ particles];

      cos3 = Select[diagramsWithCorrectTopology, MemberQ[#, particles2]&];
      If[!(cos3 === Select[diagramsWithCorrectTopology, Compare[#[[3]]/.#[[4]], temp]&]),
         Print["WHERE???"];
         Print[cos3];
         Print[InsertionsOnEdgesForDiagram[topology, diagram]];
         Print[temp];
         Print[Select[diagramsWithCorrectTopology, ((#[[3]]/.#[[4]])===temp)&]];
         Print[Select[diagramsWithCorrectTopology, Compare[#[[3]]/.#[[4]], temp]&]];
         Quit[1]
         ];
      cos3 = Select[diagramsWithCorrectTopology, Compare[#[[3]]/.#[[4]], temp]&];
      If[cos3 === {},
         Print["Topology found, but no diagram for ", diagram];
         Print[particles2];
         Print[diagramsWithCorrectTopology];
         Quit[1];
      ];

      vertices = cos3[[1, 5]];

      {"std::complex<double> " <> ToString@N[Utils`FSReIm[StripDiagramColorFactor[externalParticles, ColorFactorForDiagram[topology, diagram]]], 16] <> " * " <>
         "calculate_" <> cos3[[1, 1]],
         Length[cos3[[1, 5]]],
         cos3[[1, -1]],
         Drop[cos3[[1, 3]]/. Rule[a_, b_] :> a, 3],
         vertices,
         cos3[[1,3]]
      }

   ];


FieldFromDiagram[diagram_, {i_, j_}] :=
   If[!ListQ[diagram[[i]]],
      diagram[[i]],
      diagram[[i,j]]
   ];

(* Returns a list of particles on edges of a diagram in form of list of elements as
   {vertex1, vertex2} -> particle connecting vertex1 and 2 *)
(* @todo: not clear about the conugation of particle *)
InsertionsOnEdgesForDiagram[topology_, diagram_] := Module[{sortedVertexCombinations},

   sortedVertexCombinations =
      Select[
         Tuples[Range[CXXDiagrams`NumberOfExternalParticlesInTopology[topology] + CXXDiagrams`NumberOfPropagatorsInTopology[topology]], 2],
         (OrderedQ[#] && !SameQ@@#)&
      ];

   Print[
"THIS", diagram, " i? ",
      Flatten[
         DeleteCases[
            Switch[Length[CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology]],
               1, {{#1, #2} -> (Print[{#1, #2, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology]}];FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[1,1]]}])},
               2, {
               {#1, #2} -> FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[1,1]]}],
               {#1, #2} -> FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[2,1]]}]
            },
               _, {}
            ]& @@@ sortedVertexCombinations,
            {}
         ],
         1
      ]
   ];
      Flatten[
         DeleteCases[
            Switch[Length[CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology]],
               1, {{#1, #2} -> FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[1,1]]}]},
               2, {
                  {#1, #2} -> FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[1,1]]}],
                  {#1, #2} -> FieldFromDiagram[diagram, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, diagram, topology][[2,1]]}]
               },
               _, {}
            ]& @@@ sortedVertexCombinations,
            {}
         ],
         1
      ]
];


(* In generic_loop_decay_diagram_classes.m the couplings are written in terms of edges as
   Cp[list of fields with edge numbers as arguments of Fields.
   This function translates vertex given in terms of edges to the list of vertices connected by those edges.
   The order is preserved. *)
EdgesToVertexConnections[edges_, list_] := Module[{temp = {}, vertexNumber},

   (* for list of edges {n1, n2,... , nn} find list of pair of vertices that the edges correspond to as
      {pair of vertices connected by edge n1, pair of vertices connected by edge 2, etc} *)
   For[i = 1, i <= Length[edges], i++,
      AppendTo[temp, Select[list, MatchQ[#, {_, _} -> Field[edges[[i]]]]&]];
   ];
   Utils`AssertWithMessage[Dimensions[temp] === {Length[edges], 1}, ""];

   First /@ First  /@ temp
];

ConvertCouplingToCPP[Cp[particles__][lor_], vertices_, indices_] := Module[{
      vertexEdges, res, quit= False, pos
      },

   vertexEdges = (List[particles] /. Index[Generic, n_] :> n /. Index[Generic, n_] :> n);
   pos = First@First@Position[vertices, vertexEdges];
   res = Replace[lor,
      {LorentzProduct[_, PL] :> "left()",
         LorentzProduct[_, PR] :> "right()",
         PL :> "left()",
         PR :> "right()",
         1 :> "value()",
         Mom[U[Index[Generic, n_]]] :> (
            quit=True;
            "value(" <> ToString[Utils`MathIndexToCPP[Position[{particles}, U[Index[Generic, n]]][[1,1]]]] <> ")"),
         Mom[-U[Index[Generic, n_]]] :> (quit=True;
            "value(" <> ToString[Utils`MathIndexToCPP[Position[{particles}, -U[Index[Generic, n]]][[1,1]]]] <> ")"),
         Mom[f_[Index[Generic, n_]]] - Mom[-f_[Index[Generic, m_]]] :> "value(1,0)",
         Mom[f_[Index[Generic, n_]]] - Mom[f_[Index[Generic, m_]]] :> "value(1,0)",
         g[_, _] :> "value()",
         g[lt1_, lt2_] (-Mom[V[Index[Generic, 3]]] + Mom[V[Index[Generic, 5]]])
            + g[lt1_, lt3_] (Mom[V[Index[Generic, 3]]] - Mom[V[Index[Generic, 6]]])
            + g[lt2_, lt3_] (-Mom[V[Index[Generic, 5]]] + Mom[V[Index[Generic, 6]]]) :> "value(TripleVectorVertex::odd_permutation {})",
         g[lt1, lt2] (-Mom[V[Index[Generic, 2]]] + Mom[V[Index[Generic, 4]]])
            + g[lt1, lt3] (Mom[V[Index[Generic, 2]]] - Mom[-V[Index[Generic, 6]]])
            + g[lt2, lt3] (-Mom[V[Index[Generic, 4]]] + Mom[-V[Index[Generic, 6]]]) :> "value(TripleVectorVertex::odd_permutation {})",
         g[lt1_, lt2_] g[lt3_, lt4_] :> "value1()",
         Mom[S[n_]] - Mom[-S[Index[Generic, m_]]] :> "value(1,0)",
         lor :> (Print["Unidentifuied lorentz struct ", lor]; Quit[1])
      }
   ];

   "vertex" <> ToString@indices[[pos]] <> "::evaluate(index" <> ToString@indices[[pos]] <> ", context)." <> res
]

(* Fields for vertex given by connections *)
pppp[diagram_, topology_, connections_] := Module[{vertexNumber, temp},

   temp = CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[Sequence@@#, diagram, topology]& /@ connections;
      (
         If[ ListQ[diagram[[#2[[1]]]]],
            diagram[[#2[[1]]]][[#1[[1,1]]]], diagram[[#2[[1]]]]
   ]

      )& @@@ Transpose[{temp, connections}]
]

(* Returns translation of the form
   {Field[1] -> hh, Field[2] -> VG, Field[3] -> VG, Field[4] -> Fd, Field[5] -> bar[Fd], Field[6] -> Fd}

   A: {{1, 4} -> hh, {2, 5} -> VP, {3, 5} -> VP, {4, 5} -> Hp,

>    {4, 5} -> conj[Hp]}
B: {{1, 4} -> Field[1], {2, 5} -> Field[2], {3, 5} -> Field[3],

>    {4, 5} -> Field[4], {4, 5} -> Field[5]}
{Field[1] -> hh, Field[2] -> VP, Field[3] -> VP, Field[4] -> Hp,

>   Field[5] -> Hp}

*)
GetFieldsAssociations[concreteFieldOnEdgeBetweenVertices_, fieldNumberOnEdgeBetweenVertices_] :=
   Module[{temp = {}},

      Print[concreteFieldOnEdgeBetweenVertices];
      Print[fieldNumberOnEdgeBetweenVertices];
      temp = (Reverse /@ fieldNumberOnEdgeBetweenVertices);

      (* the numbers of vertices in Dylna might be unsorted *)
      temp = (#1 -> Sort[#2])& @@@ temp;
      Print[temp];

      For[i = 1, i <= Length[temp], i++,
         temp[[i]] = temp[[i]] /. concreteFieldOnEdgeBetweenVertices[[i]]
      ];

      Print[temp];
      temp
];

WrapCodeInLoopOverInternalVertices[topology_, diagram_, code_, internalMasses_, verticesTranslation_, fromDylan_] :=
   Module[{vertices, vertices2, indices, cppVertices, loop,
      con = {}, temp, temp2 = "", temp3, masses = "", mass = {}, temp9, kupa, fuck, translation, fieldAssociation,
      externalEdges, internalEdges,
   externalFieldsLocationsInVertices,
      internalFieldsLocationsInVertices, verticesInFieldTypes, matchExternalFieldIndicesCode, matchInternalFieldIndicesCode
   },

      translation = TranslationForDiagram[topology, diagram];

      (* {Field[1] -> concrete field, ...} *)
      fieldAssociation = GetFieldsAssociations[InsertionsOnEdgesForDiagram[topology, diagram], translation[[3]]];

      (* vertex in terms of field types (S,V,F,...) and indices 1, 2 *)
      verticesInFieldTypes =
         List @@@ (DeleteDuplicates[translation[[-3]] /. Index[Generic, n_Integer] -> n /. Cp[x___][__] -> Cp[x]]);

      (* vertices in an orientation as required by Cp *)
      Print["1 ", verticesInFieldTypes];
(*      Print["1.5", InsertionsOnEdgesForDiagram[topology, diagram]];*)
(*      Print["1.8 ", fieldAssociation, " ", ((#1 -> #2@@#1)& @@@ translation[[4]])];*)
(*      Print["2 ", fieldAssociation /. ((#1 -> #2@@#1)& @@@ translation[[4]])];*)
      vertices = verticesInFieldTypes /. (fieldAssociation /. ((#1 -> #2@@#1)& @@@ translation[[4]])) /. - e_ :> AntiField[e];

      (* set of unique indices used in names of vertices and indices *)
      indices = Table[Unique["Id"], {Length@vertices}];

      (* create using declarations for vertices *)
      cppVertices =
         "using vertex" <> ToString@#1 <> " = Vertex<" <>
            (StringJoin@Riffle[CXXDiagrams`CXXNameOfField /@ #2  ,", "] <> ">;\n")& @@@ Transpose[{indices, vertices}];

      (* loop over indices *)
      loop = "for(const auto& index" <> ToString@# <> ": index_range<vertex" <> ToString@# <> ">()) {\n"& /@ indices;

      (* List of {integer, integer} -> Field[integer] *)
      externalEdges =
         Select[
            translation[[3]],
            (MatchQ[#, ({i_Integer, j_Integer} -> Field[_Integer]) /; (First@Sort[{i,j}] >= 1 && First@Sort[{i,j}] <= 3 && Last@Sort[{i,j}] > 3)])&
         ];

      matchExternalFieldIndicesCode =
         "auto externalFieldIndicesIn" <> ToString[#] <>
            " = vertex" <>
            ToString@indices[[ First@First@Position[verticesInFieldTypes, _[#]] ]] <>
         "::template indices_of_field<" <>
            ToString@Utils`MathIndexToCPP@Last@First@Position[verticesInFieldTypes, _[#]] <>
            ">(index" <>
            ToString@indices[[ First@First@Position[verticesInFieldTypes, _[#]] ]] <> ");\n" & /@ Range[3];

      matchInternalFieldIndicesCode =
         ("if(vertex" <> ToString@indices[[  First@First@Position[verticesInFieldTypes, _[#1]]   ]] <>
         "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@ Last@First@Position[verticesInFieldTypes, _[#1]] <> ">(index" <>
         ToString@indices[[  First@First@Position[verticesInFieldTypes, _[#1]]   ]] <> ") != " <>
      "vertex" <> ToString@indices[[  First@First@Position[verticesInFieldTypes, _[#2]]   ]] <>
         "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@ Last@First@Position[verticesInFieldTypes, _[#2]] <> ">(index" <>
         ToString@indices[[  First@First@Position[verticesInFieldTypes, _[#2]]   ]] <> ")) {\n" <> TextFormatting`IndentText["continue;\n"] <> "}\n")&
      @@@
      DeleteCases[DeleteDuplicates[Sort /@ Tuples[Range[4, 3+Length[vertices]], 2]], {n_Integer, n_Integer}];

      mass =(
      "const auto mInternal" <> ToString[#-3] <> " = context.mass<" <>
         CXXNameOfField[
            vertices[[  Sequence@@First@Position[verticesInFieldTypes /. -field_ -> field, _[#]]    ]]
            ] <> ">(" <>
         "vertex" <> ToString@indices[[First@First@Position[verticesInFieldTypes/. -field_->field, _[#]]]] <> "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP[
         Last@First@Position[verticesInFieldTypes, _[#]]
         ] <> ">(index" <> ToString@indices[[First@First@Position[verticesInFieldTypes/.-field_->field, _[#]]]] <> "));\n"
)&/@
      Range[4, 3+Length@vertices];

      mass = StringJoin@@mass;

      internalEdges =
         Select[
            translation[[3]],
            (MatchQ[#, ({i_Integer, j_Integer} -> Field[_Integer]) /; (First@Sort[{i,j}] > 3 && Last@Sort[{i,j}] > 3)])&
         ];
      internalFieldsLocationsInVertices =
         {
            {First@#1, First@First@(Position[translation[[-2, First@#1]] /. Susyno`LieGroups`conj -> Identity, #2])},
            {Last@#1,  First@First@(Position[translation[[-2, Last@#1 ]] /. Susyno`LieGroups`conj -> Identity, #2])}

}& @@@ internalEdges;


            (* write number of propagators in topology *)
            "// topology with " <> ToString@CXXDiagrams`NumberOfPropagatorsInTopology[topology] <> " propagators\n" <>
               cppVertices <>
               loop <> "\n" <>
               TextFormatting`IndentText[

                  (* skip indices that don't match external indices *)
                  matchExternalFieldIndicesCode <>
                     "if(externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3) {\n" <>
                     TextFormatting`IndentText["continue;"] <>
                     "\n}\n\n" <>

                  matchInternalFieldIndicesCode <>
                  mass <>

                  "result += calculate_" <> translation[[1]] <> "(\n" <>
                  TextFormatting`IndentText[
                  "result.m_decay, result.m_out_1, result.m_out_2,\n" <>
                    StringJoin @@ Riffle[("mInternal" <> ToString@#)& /@ Range@CXXDiagrams`NumberOfPropagatorsInTopology[topology], ", "] <> ", " <> "\n"
                     ] <>
                     (* couplings *)
                     StringJoin @@ Riffle[ToString /@ ConvertCouplingToCPP[#, verticesInFieldTypes, indices]& /@ translation[[-3]], ", "] <> ",\n" <>



(*                  translation[[-3]] <>*)
(*                     StringJoin@Riffle[("vertex" <> ToString@indices[[First[Intersection@@#2]-3]] <> "::evaluate(index" <> ToString@indices[[First[Intersection@@#2]-3]] <> ", context)."
 <> ToString[#1])& @@@ kupa, ", "] <> ",\n" <>*)
                    code[[2]]

               ] <>
               "\n" <>
               StringJoin@@ConstantArray["}\n", Length@vertices] <>
               "\n"
   ];

EvaluateDecayDiagramWithTopology[decay_, topology_, diagram_] :=
    Which[IsOneLoopTwoBodyDecayTopology[topology],
          EvaluateOneLoopTwoBodyDecayDiagramWithTopology[decay, topology, diagram],
          True,
          Print["Error: requested evaluation of unsupported topology."];
          Quit[1];
         ];

ColorFactorForDiagram[topology_, diagram_] :=
   ColourFactorForIndexedDiagramFromGraph[
      CXXDiagrams`IndexDiagramFromGraph[diagram, topology], topology
   ];

FillOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{oneLoopTopAndInsertion, body = "", temp},

       (* list of elements like {topology, insertion (diagram)} *)
       oneLoopTopAndInsertion = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];

           (
           temp = TranslationForDiagram[Sequence @@ #];
              body = body <>
                 WrapCodeInLoopOverInternalVertices[
                    Sequence @@ #,
                    {"result += calculate_" <> temp[[1]] <>
                       "(result.m_decay, result.m_out_1, result.m_out_2,\n" <>
                          "// number of internal masses " <> ToString@CXXDiagrams`NumberOfPropagatorsInTopology[#[[1]]] <> "\n" <>
                           "mInternal1, mInternal2," <> If[CXXDiagrams`NumberOfPropagatorsInTopology[#[[1]]] === 3, " mInternal3,", ""] <> "\n" <>

                    "// number of couplings " <> ToString@EvaluateDecayDiagramWithTopology[decay, Sequence @@ #][[2]] <> "\n",

                     (* scale *)
                    "// renormalization scale\n" <>
                    "result.m_decay\n" <>
                       (* finite or not *)
                    "// finite part\n" <>
                       If[!EvaluateDecayDiagramWithTopology[decay, Sequence @@ #][[3]], ", 1.", ""] <>
                     ");"},
               EvaluateDecayDiagramWithTopology[decay, Sequence @@ #][[4]],
               EvaluateDecayDiagramWithTopology[decay, Sequence @@ #][[5]],
                    EvaluateDecayDiagramWithTopology[decay, Sequence @@ #][[6]]

                  ]
            )& /@ oneLoopTopAndInsertion;

           "// 1-loop amplitude(s)\n" <> body
          ];

(*FillSFFOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)
    (*Module[{oneLoopDiags, body = ""},*)
           (*oneLoopDiags = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];*)
           (*(body = body <>*)
               (*"result.form_factor_left += " <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #] <> ";\n" <>*)
               (*"result.form_factor_right += " <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #] <> ";\n"*)
            (*)& /@ oneLoopDiags;*)
           (*"// 1-loop amplitude(s)\n" <> body*)
          (*];*)

(*FillSSVOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)
    (*Module[{oneLoopDiags, body = ""},*)
           (*oneLoopDiags = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];*)
           (*(body = body <>*)
              (*"result.form_factor += " <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #] <> ";\n"*)
            (*)& /@ oneLoopDiags;*)
           (*"// 1-loop amplitude(s)\n" <> body*)
          (*];*)

(*FillSVVOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)
    (*Module[{oneLoopDiags, body = ""},*)
           (*oneLoopDiags = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];*)
           (*(body = body <> *)
               (*"result.form_factor_g += " <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #] <> ";\n"*)

            (*)& /@ oneLoopDiags;*)
           (*"// 1-loop amplitude(s)\n" <> body*)
          (*];*)

(*FillFFSOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)
    (*Module[{oneLoopDiags, body = ""},*)
           (*oneLoopDiags = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];*)
           (*(body = body <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #])& /@ oneLoopDiags;*)
           (*"// 1-loop amplitude(s)\n" <> body*)
          (*];*)

(*FillFFVOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)
    (*Module[{oneLoopDiags, body = ""},*)
           (*oneLoopDiags = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];*)
           (*(body = body <> EvaluateDecayDiagramWithTopology[decay, Sequence @@ #])& /@ oneLoopDiags;*)
           (*"// 1-loop amplitude(s)\n" <> body*)
          (*];*)

(*FillOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=*)

(*
    Switch[GetDecayAmplitudeType[decay],
           "Decay_amplitude_SSS", FillSSSOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SFF", FillSFFOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SSV", FillSSVOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SVV", FillSVVOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFS", FillFFSOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFV", FillFFVOneLoopDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           _, ""
          ];
          *)

(* creates `calculate_amplitude` function
   that returns a total (sumed over internal insertions) 1-loop amplitude for a given external particles *)
CreateTotalAmplitudeSpecializationDef[decay_FSParticleDecay, modelName_] :=
   Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsNamespace = modelName <> "_cxx_diagrams::fields",
            returnVar = "result", paramsStruct = "context", returnType = "",
            externalFieldsList, templatePars = "", args = "",
            body = ""},

           Print["Entry point: creating amplitude for ", initialParticle, " -> ", finalState]; (* @todo: remove *)

           (* Decay_amplitude_XXX *)
           returnType = GetDecayAmplitudeType[decay];

           (* template arguments *)
           externalFieldsList = Join[{initialParticle}, finalState];
           templatePars = "<" <> Utils`StringJoinWithSeparator[
              CXXDiagrams`CXXNameOfField[#, prefixNamespace -> fieldsNamespace]& /@ externalFieldsList, ", "] <> ">";

           (* function arguments *)
           args =
              "const " <> modelName <> "_cxx_diagrams::context_base& " <> paramsStruct <> ", " <>
                 Utils`StringJoinWithSeparator[
                     MapIndexed[(CreateFieldIndices[#1, fieldsNamespace] <> " const& idx_" <> ToString[First[#2]])&, externalFieldsList],
                     ", "
                 ];

           (* body *)
           body = returnType <> " " <> returnVar <> ";\n";
           body = body <> FillDecayAmplitudeMasses[decay, modelName, returnVar, paramsStruct] <> "\n"; (* external particle masses *)
           body = body <> ZeroDecayAmplitudeFormFactors[decay, returnVar] <> "\n";
           If[IsPossibleTreeLevelDecay[decay, True],
              body = body <> "// @todo correct prefactors\n" <> FillTreeLevelDecayAmplitudeFormFactors[decay, modelName, returnVar, paramsStruct] <> "\n";
             ];
           If[!IsPossibleTreeLevelDecay[decay, True] && IsPossibleOneLoopDecay[decay],
              body = body <> FillOneLoopDecayAmplitudeFormFactors[decay, modelName, returnVar, paramsStruct] <> "\n";
             ];
           body = body <> "return " <> returnVar <> ";\n";

            "// " <> ToString@initialParticle <>  " -> " <> ToString@finalState <> "\n" <>
           "template<>\n" <>
              returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
               templatePars <>
                  "(" <> args <> ") const\n{\n" <>
                     TextFormatting`IndentText[body] <>
                  "}\n"
   ];

GetHiggsBosonDecays[particleDecays_List] :=
    If[TreeMasses`GetHiggsBoson =!= Null, Select[particleDecays, (First[#] === TreeMasses`GetHiggsBoson[])&], {}];

SelectDecayByFinalState[finalState_List, decays_List] :=
    Select[decays, (Sort[GetFinalState[#]] === Sort[finalState])&];

SelectDownQuarkDownQuarkFinalState[decays_List] :=
    Module[{downQuarkSymbol, result = {}},
           downQuarkSymbol = TreeMasses`GetDownQuark[1] /. field_[generation_] :> field;
           If[downQuarkSymbol =!= Null,
              result = SelectDecayByFinalState[{downQuarkSymbol, SARAH`AntiField[downQuarkSymbol]}, decays];
             ];
           result
          ];

SelectChargedLeptonChargedLeptonFinalState[decays_List] :=
    Module[{chargedLeptonSymbol, result = {}},
       chargedLeptonSymbol = TreeMasses`GetDownLepton[1] /. field_[generation_] :> field;
       If[chargedLeptonSymbol =!= Null,
          result = SelectDecayByFinalState[{chargedLeptonSymbol, SARAH`AntiField[chargedLeptonSymbol]}, decays];
       ];
       result
    ];

SelectGluonGluonFinalState[decays_List] :=
    Module[{gluonSymbol = TreeMasses`GetGluon[], result = {}},
           If[gluonSymbol =!= Null,
              result = SelectDecayByFinalState[{gluonSymbol, gluonSymbol}, decays];
             ];
           result
          ];

SelectHiggsHiggsFinalState[decays_List] :=
    Module[{higgsSymbol = TreeMasses`GetHiggsBoson[], result = {}},
           If[higgsSymbol =!= Null,
              result = SelectDecayByFinalState[{higgsSymbol, higgsSymbol}, decays];
             ];
           result
          ];

SelectPhotonPhotonFinalState[decays_List] :=
    Module[{photonSymbol = TreeMasses`GetPhoton[], result = {}},
           If[photonSymbol =!= Null,
              result = SelectDecayByFinalState[{photonSymbol, photonSymbol}, decays];
             ];
           result
          ];

SelectHiggsHiggsFinalState[decays_List] :=
    Module[{psSymbol = TreeMasses`GetPseudoscalarHiggsBoson[], result = {}},
           If[psSymbol =!= Null,
              result = SelectDecayByFinalState[{psSymbol, psSymbol}, decays];
             ];
           result
          ];

SelectUpQuarkUpQuarkFinalState[decays_List] :=
    Module[{upQuarkSymbol, result = {}},
           upQuarkSymbol = TreeMasses`GetUpQuark[1] /. field_[generation_] :> field;
           If[upQuarkSymbol =!= Null,
              result = SelectDecayByFinalState[{upQuarkSymbol, SARAH`AntiField[upQuarkSymbol]}, decays];
             ];
           result
          ];

SelectWWFinalState[decays_List] :=
    Module[{wBosonSymbol = TreeMasses`GetWBoson[], result = {}},
           If[wBosonSymbol =!= Null,
              result = SelectDecayByFinalState[{wBosonSymbol, SARAH`AntiField[wBosonSymbol]}, decays];
             ];
           result
          ];

SelectZPhotonFinalState[decays_List] :=
    Module[{zBosonSymbol = TreeMasses`GetZBoson[], photonSymbol = TreeMasses`GetPhoton[], result = {}},
           If[zBosonSymbol =!= Null && photonSymbol =!= Null,
              result = SelectDecayByFinalState[{zBosonSymbol, photonSymbol}, decays];
             ];
           result
          ];

SelectZZFinalState[decays_List] :=
    Module[{zBosonSymbol = TreeMasses`GetZBoson[], result = {}},
           If[zBosonSymbol =!= Null,
              result = SelectDecayByFinalState[{zBosonSymbol, zBosonSymbol}, decays];
             ];
           result
          ];

CreateHiggsToGluonGluonTotalAmplitudeFunction[hggDecay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[hggDecay], finalState = GetFinalState[hggDecay],
            fieldsList, returnType = "", args = "", templatePars = "", body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#]] <> "&")& /@ fieldsList, ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = returnType <> " result;\nreturn result;\n";
           "template<>\n" <> returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreateHiggsToGluonGluonTotalAmplitude[particleDecays_List, modelName_] :=
    Module[{higgsDecays, hggDecay, prototype = "", function = ""},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              hggDecay = SelectGluonGluonFinalState[Last[higgsDecays]];
              If[hggDecay =!= {},
                 hggDecay = First[hggDecay];
                 prototype = CreateTotalAmplitudeSpecializationDecl[hggDecay, modelName];
                 function  = CreateHiggsToGluonGluonTotalAmplitudeFunction[hggDecay, modelName]
                ];
             ];
           {prototype, function}
          ];

CreateHiggsToPhotonPhotonTotalAmplitudeFunction[hgamgamDecay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[hgamgamDecay], finalState = GetFinalState[hgamgamDecay],
            fieldsList, returnType = "", args = "", templatePars = "", body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#]] <> "&")& /@ fieldsList, ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = returnType <> " result;\nreturn result;\n";
           "template<>\n" <> returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreateHiggsToPhotonPhotonTotalAmplitude[particleDecays_List, modelName_] :=
    Module[{higgsDecays, hgamgamDecay, prototype = "", function = ""},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              hgamgamDecay = SelectPhotonPhotonFinalState[Last[higgsDecays]];
              If[hgamgamDecay =!= {},
                 hgamgamDecay = First[hgamgamDecay];
                 prototype = CreateTotalAmplitudeSpecializationDecl[hgamgamDecay, modelName];
                 function  = CreateHiggsToPhotonPhotonTotalAmplitudeFunction[hgamgamDecay, modelName]
                ];
             ];
           {prototype, function}
          ];

CreateTotalAmplitudeSpecialization[decay_FSParticleDecay, modelName_] :=
    Module[{decl = "", def = ""},
           decl = CreateTotalAmplitudeSpecializationDecl[decay, modelName];
           def = CreateTotalAmplitudeSpecializationDef[decay, modelName];
           {decl, def}
          ];

CreateTotalAmplitudeSpecializations[particleDecays_List, modelName_] :=
    Module[{specializations},
           specializations = Flatten[(CreateTotalAmplitudeSpecialization[#, modelName]& /@ Last[#])& /@ particleDecays, 1];
           specializations = Select[specializations, (# =!= {} && # =!= {"", ""})&];
           Utils`StringJoinWithSeparator[#, "\n"]& /@ Transpose[specializations]
          ];

CreatePartialWidthSpecializationDecl[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, fieldsNamespace, args},
           fieldsList = Join[{initialParticle}, finalState];
           fieldsNamespace = If[modelName != "", modelName <> "_cxx_diagrams::fields", False];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[#, fieldsNamespace] <> "&")& /@ fieldsList, ", "];
           "template <>\n" <>
           "double " <> modelName <> "_decays::" <>
           CreateSpecializedPartialWidthCalculationName[initialParticle, finalState, fieldsNamespace] <>
           "(" <> args <> ") const;"
          ];

CreateIncludedPartialWidthSpecialization[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            declaration = "", includeStatement = ""},
           declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
           includeStatement = "#include \"templates/sm_h_decays/decay_" <>
                              SimplifiedName[initialParticle] <> "_to_" <>
                              StringJoin[SimplifiedName[# /. SARAH`bar|Susyno`LieGroups`conj -> Identity]& /@ finalState] <>
                              ".cpp\"";
           {declaration, includeStatement}
          ];

CreateHiggsToZZPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectZZFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToWWPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectWWFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToGluonGluonPartialWidthFunction[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, args, templatePars, body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           fieldIndices = {"in_idx", "out1_idx", "out2_idx"};
           args = "const " <> modelName <> "_cxx_diagrams::context_base& context, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#1]] <> "& " <> #2)& @@@ Transpose[{fieldsList, fieldIndices}], ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = 
              "const auto amp = calculate_amplitude<H, G, G>(context, in_idx, out1_idx, out2_idx);\n" <>
              "return amp.square();\n";
           "template <>\ndouble CLASSNAME::" <> CreateGenericGetPartialWidthFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateHiggsToPhotonPhotonPartialWidthFunction[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, args, templatePars, body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           fieldIndices = {"in_idx", "out1_idx", "out2_idx"};
           args = "const " <> modelName <> "_cxx_diagrams::context_base& context, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#1]] <> "& " <> #2)& @@@ Transpose[{fieldsList, fieldIndices}], ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = 
              "const auto amp = calculate_amplitude<H, A, A>(context, in_idx, out1_idx, out2_idx);\n" <>
              "return amp.square();\n";
           "template <>\ndouble CLASSNAME::" <> CreateGenericGetPartialWidthFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateHiggsToGluonGluonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectGluonGluonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
              function = CreateHiggsToGluonGluonPartialWidthFunction[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToPhotonPhotonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectPhotonPhotonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
              function = CreateHiggsToPhotonPhotonPartialWidthFunction[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToZPhotonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, orderZPhotonFinalState, declaration = "", function = ""},
           decay = SelectZPhotonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              orderZPhotonFinalState[finalState_] :=
                  Module[{zBosonSymbol = TreeMasses`GetZBoson[], photonSymbol = TreeMasses`GetPhoton[]},
                         If[finalState === {zBosonSymbol, photonSymbol},
                            finalState,
                            Reverse[finalState]
                           ]
                        ];
              decay = ReplacePart[decay, {{2}} -> orderZPhotonFinalState[GetFinalState[decay]]];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToHiggsHiggsPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectHiggsHiggsFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToPseudoscalarPseudoscalarPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectPseudoscalarPseudoscalarFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToUpQuarkUpQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectUpQuarkUpQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToDownQuarkDownQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectDownQuarkDownQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];
CreateHiggsToChargedLeptonChargedLeptonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
       decay = SelectChargedLeptonChargedLeptonFinalState[decaysList];
       If[decay =!= {},
          decay = First[decay];
          {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
       ];
       {declaration, function}
    ];

CreateHiggsDecayPartialWidthSpecializations[particleDecays_, modelName_] :=
    Module[{higgsDecays, specializations = {}},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              specializations = {CreateHiggsToZZPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToWWPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToGluonGluonPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToPhotonPhotonPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToZPhotonPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToHiggsHiggsPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToUpQuarkUpQuarkPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToDownQuarkDownQuarkPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToChargedLeptonChargedLeptonPartialWidth[higgsDecays, modelName]};
             ];
           specializations
          ];

CreatePartialWidthSpecializations[particleDecays_List, modelName_] :=
    Module[{specializations},
           specializations = CreateHiggsDecayPartialWidthSpecializations[particleDecays, modelName];
           specializations = Select[specializations, (# =!= {} && # =!= {"", ""})&];
           Utils`StringJoinWithSeparator[#, "\n\n"]& /@ Transpose[specializations]
          ];

CreateDecaysGetterFunctionName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysGetterFunction[particle_] :=
    Module[{dim, body = ""},
           dim = TreeMasses`GetDimension[particle];
           body = "return decay_table." <> CreateDecayTableEntryGetterName[particle] <>
                  "(" <> If[dim > 1, "i", ""] <> ");";
           "const Decays_list& " <> CreateDecaysGetterFunctionName[particle] <> "(" <>
           If[dim > 1, "int i", ""] <> ") const { " <> body <> " }"
          ];

CreateDecaysGetterFunctions[particles_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysGetterFunction /@ particles, "\n"];

CreateDecayTableEntryGetterName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecayTableEntryGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           "Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ");"
          ];

CreateDecayTableEntryConstGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           "const Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ") const;"
          ];

CreateDecayTableEntryGetterFunctionBody[particle_, rows_List] :=
    Module[{i, dim, idxName = "gI1", errMsg = "", body = ""},
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           If[dim != Length[rows],
              Print["Error: number of rows (", Length[rows], ") does not match size of"];
              Print["    ", particle, " multiplet."];
              Quit[1];
             ];
           If[dim == 1,
              body = "return decay_table[" <> ToString[rows[[1]]] <> "];\n";
              ,
              body = "switch (" <> idxName <> ") {\n";
              For[i = 0, i < dim, i++,
                  body = body <> "case " <> ToString[i] <> ": return decay_table[" <> ToString[rows[[i + 1]]] <> "]; break;\n";
                 ];
              body = body <> "}\n\n";
              errMsg = "std::ostringstream sstr;\n" <>
                       "sstr << \"invalid particle index \" << std::to_string(" <> idxName <> ") << '\\n';\n\n" <>
                       "throw OutOfBoundsError(sstr.str());\n";
              body = body <> errMsg;
             ];
           body
          ];

CreateDecayTableEntryGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ")\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryConstGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "const Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryGetterFunction[particle, {row}, scope];

CreateDecayTableEntryConstGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryConstGetterFunction[particle, {row}, scope];

CreateDecayTableGetterPrototypes[decayParticles_List] :=
    Utils`StringJoinWithSeparator[(CreateDecayTableEntryGetterPrototype[#] <> "\n" <>
                                   CreateDecayTableEntryConstGetterPrototype[#])& /@ decayParticles, "\n"];

CreateDecayTableGetterFunctions[decayParticles_List, scope_:"CLASSNAME"] :=
    Module[{i, dims, offsets, rowAssignments, defs = ""},
           dims = TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles;
           offsets = If[Length[dims] == 1, {0}, Join[{0}, Accumulate[dims[[1;;-1]]]]];
           rowAssignments = MapIndexed[{decayParticles[[First[#2]]], Table[offsets[[First[#2]]] + i, {i, 0, #1 - 1}]}&, dims];
           defs = (CreateDecayTableEntryGetterFunction[#[[1]], #[[2]], scope] <> "\n\n" <>
                   CreateDecayTableEntryConstGetterFunction[#[[1]], #[[2]], scope])& /@ rowAssignments;
           Utils`StringJoinWithSeparator[defs, "\n"]
          ];

CreateDecayTableInitialization[decayParticles_List] :=
    Module[{i, dims, dimsWithoutGoldstones, starts, pdgCodes, initializerList = ""},
           dims = TreeMasses`GetDimension /@ decayParticles;
           dimsWithoutGoldstones = TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles;
           starts = TreeMasses`GetDimensionStartSkippingGoldstones /@ decayParticles;
           pdgCodes = Parameters`GetPDGCodesForParticle /@ decayParticles;
           For[i = 1, i <= Length[decayParticles], i++,
               If[dims[[i]] != Length[pdgCodes[[i]]],
                  Print["Error: number of PDG codes does not match size of ", decayParticles[[i]], " multiplet."];
                  Quit[1];
                 ];
               If[dimsWithoutGoldstones[[i]] > 0,
                  initializerList = initializerList <> If[initializerList == "", "", ", "] <>
                                    Utils`StringJoinWithSeparator[("Decays_list(" <> ToString[#] <> ")")& /@ pdgCodes[[i, starts[[i]] ;;]], ", "];
                 ];
              ];
           ": decay_table({" <> initializerList <> "})"
          ];

End[]

EndPackage[];
