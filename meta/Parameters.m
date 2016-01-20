
BeginPackage["Parameters`", {"SARAH`", "CConversion`", "Utils`", "Phases`"}];

FindSymbolDef::usage="";

CreateSetAssignment::usage="";
CreateDisplayAssignment::usage="";
CreateParameterNamesStr::usage="";
CreateParameterEnums::usage="";
CreateInputParameterEnum::usage="";
CreateInputParameterNames::usage="";
CreateStdVectorNamesAssignment::usage="";

CreateInputParameterArrayGetter::usage="";
CreateInputParameterArraySetter::usage="";

SetParameter::usage="set model parameter";
SetInputParameter::usage="set input parameter";
AddInputParameters::usage="add an input parameter";
SetPhases::usage="sets field phases";
GetPhases::usage="returns field phases";
SetPhase::usage="sets a phase to a value";

ApplyGUTNormalization::usage="Returns a list of repacement rules for
gauge couplings, which replace non-normalized gauge couplings
(e.g. gY) by normalized ones (e.g. g1).";

GetGUTNormalization::usage="Returns GUT normalization of the given
coupling";

CreateIndexReplacementRules::usage="";

AddRealParameter::usage="";
SetRealParameters::usage="";

SaveParameterLocally::usage="Save parameters in local variables";
RestoreParameter::usage="Restore parameters from local variables";

GetType::usage="";
GetPhase::usage="";
HasPhase::usage="";
GetRealTypeFromDimension::usage="";
GetParameterDimensions::usage="";
GetThirdGeneration::usage="returns parameter with third generation index";

IsRealParameter::usage="";
IsComplexParameter::usage="";
IsRealExpression::usage="";
IsMatrix::usage="returns true if parameter is a matrix";
IsSymmetricMatrixParameter::usage="returns true if parameter is a matrix";
IsModelParameter::usage="returns True if parameter is a model parameter";
IsInputParameter::usage="returns False if parameter is an input parameter";
IsOutputParameter::usage="returns True if parameter is a defined output parameter";
IsIndex::usage="returns true if given symbol is an index";
IsPhase::usage="returns True if given symbol is a phase";

GetIndices::usage="returns list of indices from a given parameter";

AllModelParametersAreReal::usage="returns True if all model parameters
are real, False otherwise";

SetInputParameters::usage="";
SetModelParameters::usage="";
SetOutputParameters::usage="";

GetInputParameters::usage="";
GetModelParameters::usage="";
GetOutputParameters::usage="";
GetModelParametersWithMassDimension::usage="Returns model parameters
with given mass dimension";

GetDependenceNumSymbols::usage="Returns symbols which have a
DependenceNum";

CreateLocalConstRefs::usage="creates local const references to model
parameters / input parameters.";

CreateLocalConstRefsForBetas::usage="";

CreateLocalConstRefsForInputParameters::usage="creates local const
references for all input parameters in the given expression.";

CreateLocalConstRefsForPhysicalParameters::usage="creates local const
references for all physical output parameters in the given
expression.";

FillInputParametersFromTuples::usage="";

DecreaseIndexLiterals::usage="";
IncreaseIndexLiterals::usage="";

DecreaseSumIdices::usage="";

GetEffectiveMu::usage="";
GetEffectiveMASqr::usage="";

GetParameterFromDescription::usage="Returns model parameter from a
given description string.";
GetParticleFromDescription::usage="Returns particle symbol from a
given description string.";
GetPDGCodesForParticle::usage="Returns the PDG codes for a particle."

NumberOfIndependentEntriesOfSymmetricMatrix::usage="Returns number of
independent parameters of a real symmetric nxn matrix";

ExpandExpressions::usage="";
AppendGenerationIndices::usage="";

StripIndices::usage="removes indices from model parameter";

FilterOutLinearDependentEqs::usage="returns linear independent equations";

FilterOutIndependentEqs::usage = "returns equations that depend on the
given list of parameters.  I.e. equations, that do not depend on the
given list of parameters are omitted from the output.";

FindAllParameters::usage = "returns list of all parameters contained
in the given expression";

Begin["`Private`"];

allInputParameters = {};
allModelParameters = {};
allOutputParameters = {};
allPhases = {};

SetInputParameters[pars_List] := allInputParameters = DeleteDuplicates[pars];
AddInputParameters[pars_List] := allInputParameters = DeleteDuplicates[Utils`ForceJoin[allInputParameters, pars]];
SetModelParameters[pars_List] := allModelParameters = DeleteDuplicates[pars];
SetOutputParameters[pars_List] := allOutputParameters = DeleteDuplicates[pars];
SetPhases[phases_List]        := allPhases = DeleteDuplicates[phases];

GetInputParameters[] := allInputParameters;
GetModelParameters[] := allModelParameters;
GetOutputParameters[] := allOutputParameters;
GetPhases[] := allPhases;

additionalRealParameters = {};

AddRealParameter[parameter_List] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, parameter]];

AddRealParameter[parameter_] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, {parameter}]];

SetRealParameters[parameters_List] :=
    additionalRealParameters = parameters;

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<Parameters>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

FindSymbolDef[sym_] :=
    Module[{symDef},
           symDef = Cases[SARAH`ParameterDefinitions,
                          {sym, {___, DependenceNum -> definition_, ___}} :> definition];
           If[Head[symDef] =!= List || symDef === {},
              Print["Error: Could not find definition of ",
                    sym, " in SARAH`ParameterDefinitions"];
              Return[0];
             ];
           If[Length[symDef] > 1,
              Print["Warning: ", sym, " defined multiple times"];
             ];
           symDef = symDef[[1]];
           Return[symDef];
          ];

(* Returns all parameters within an expression *)
FindAllParameters[expr_] :=
    Module[{symbols, compactExpr, allParameters},
           allParameters = Join[allModelParameters, allOutputParameters,
                                allInputParameters, Phases`GetArg /@ allPhases,
                                GetDependenceNumSymbols[]];
           compactExpr = RemoveProtectedHeads[expr];
           (* find all model parameters with SARAH head *)
           symbols = DeleteDuplicates[Flatten[
               { Cases[compactExpr, SARAH`L[a_][__] /; MemberQ[allParameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_][__] /; MemberQ[allParameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_][__] /; MemberQ[allParameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_][__] /; MemberQ[allParameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`L[a_]     /; MemberQ[allParameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_]     /; MemberQ[allParameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_]     /; MemberQ[allParameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_]     /; MemberQ[allParameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}]
               }]];
           (* remove parameters found from compactExpr *)
           compactExpr = compactExpr /. (RuleDelayed[#, CConversion`ToValidCSymbolString[#]]& /@ symbols);
           (* find all model parameters without SARAH head in compactExpr *)
           symbols = Join[symbols,
               { Cases[compactExpr, a_Symbol /; MemberQ[allParameters,a], {0,Infinity}],
                 Cases[compactExpr, a_[__] /; MemberQ[allParameters,a] :> a, {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]], {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
               }];
           DeleteDuplicates[Flatten[symbols]]
          ];

IsScalar[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 1 || Length[SARAH`getDimParameters[sym]] == 0;

IsMatrix[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsMatrix[sym];

IsMatrix[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 2;

IsSymmetricMatrixParameter[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsSymmetricMatrixParameter[sym];

IsSymmetricMatrixParameter[sym_] :=
    IsMatrix[sym] && MemberQ[SARAH`ListSoftBreakingScalarMasses, sym];

AllModelParametersAreReal[] := MemberQ[SARAH`RealParameters, All];

sarahIndices = {
    SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4,
    Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4
};

IsIndex[i_?NumberQ] := True;
IsIndex[i_ /; MemberQ[sarahIndices,i]] := True;
IsIndex[_] := False;
IsIndex[indices_List] := And @@ (IsIndex /@ indices);
IsIndex[indices__] := IsIndex[{indices}];

GetIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] := {indices};
GetIndices[parameter_] := {};

IsPhase[parameter_] := MemberQ[Phases`GetArg /@ allPhases, Phases`GetArg[parameter]];

IsModelParameter[parameter_] := MemberQ[allModelParameters, parameter];
IsModelParameter[Re[parameter_]] := IsModelParameter[parameter];
IsModelParameter[Im[parameter_]] := IsModelParameter[parameter];

IsModelParameter[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    IsModelParameter[parameter];

IsInputParameter[parameter_] := MemberQ[allInputParameters, parameter];

IsOutputParameter[lst_List] := And @@ (IsOutputParameter /@ lst);
IsOutputParameter[sym_]     := MemberQ[GetOutputParameters[],sym];

IsRealParameter[Re[sym_]] := True;
IsRealParameter[Im[sym_]] := True;

IsRealParameter[sym_] :=
    (IsModelParameter[sym] && AllModelParametersAreReal[]) ||
    MemberQ[Utils`ForceJoin[SARAH`realVar, additionalRealParameters, SARAH`RealParameters], sym];

IsComplexParameter[sym_] :=
    !IsRealParameter[sym];

IsRealExpression[parameter_ /; IsModelParameter[parameter]] :=
    IsRealParameter[parameter];

IsRealExpression[expr_?NumericQ] :=
    Element[expr, Reals];

IsRealExpression[_Complex] := False;

IsRealExpression[_Real] := True;

IsRealExpression[expr_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsRealExpression[expr];

IsRealExpression[HoldPattern[SARAH`Delta[_,_]]] := True;
IsRealExpression[HoldPattern[SARAH`ThetaStep[_,_]]] := True;
IsRealExpression[Cos[_]] := True;
IsRealExpression[Sin[_]] := True;
IsRealExpression[ArcCos[_]] := True;
IsRealExpression[ArcSin[_]] := True;
IsRealExpression[ArcTan[_]] := True;
IsRealExpression[Power[a_,b_]] := IsRealExpression[a] && IsRealExpression[b];
IsRealExpression[Susyno`LieGroups`conj[expr_]] :=
    IsRealExpression[expr];
IsRealExpression[SARAH`Conj[expr_]] := IsRealExpression[expr];
IsRealExpression[Conjugate[expr_]]  := IsRealExpression[expr];
IsRealExpression[Transpose[expr_]]  := IsRealExpression[expr];
IsRealExpression[SARAH`Tp[expr_]]   := IsRealExpression[expr];
IsRealExpression[SARAH`Adj[expr_]]  := IsRealExpression[expr];
IsRealExpression[bar[expr_]]        := IsRealExpression[expr];

IsRealExpression[expr_Symbol] := IsRealParameter[expr];

IsRealExpression[Times[Conjugate[a_],a_]] := True;
IsRealExpression[Times[a_,Conjugate[a_]]] := True;
IsRealExpression[Times[Susyno`LieGroups`conj[a_],a_]] := True;
IsRealExpression[Times[a_,Susyno`LieGroups`conj[a_]]] := True;
IsRealExpression[Times[SARAH`Conj[a_],a_]]            := True;
IsRealExpression[Times[a_,SARAH`Conj[a_]]]            := True;

IsRealExpression[factors_Times] :=
    And @@ (IsRealExpression[#]& /@ (List @@ factors));

IsRealExpression[terms_Plus] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsHermitian[a_] :=
    Or[IsScalar[a] && IsRealParameter[a],
       IsSymmetricMatrixParameter[a] && IsRealParameter[a],
       MemberQ[SARAH`ListSoftBreakingScalarMasses, a]
      ];

(* helper function which calculates the adjoint of an expression *)
FSAdj[Susyno`LieGroups`conj[a_]] := SARAH`Tp[a];
FSAdj[SARAH`Tp[a_]] := a /; IsRealParameter[a];
FSAdj[SARAH`Tp[a_]] := Susyno`LieGroups`conj[a];
FSAdj[a_] := a /; IsHermitian[a];
FSAdj[a_] := SARAH`Adj[a];
FSAdj[a__] := Sequence @@ (FSAdj /@ Reverse[{a}]);

TraceEquality[SARAH`trace[a__], SARAH`trace[b__]] :=
    Or @@ (({a} === #)& /@ NestList[RotateLeft, {b}, Length[{b}] - 1]);

(* if all parameters in the trace are real, the trace is real *)
IsRealExpression[SARAH`trace[a__]] :=
    True /; And @@ (IsRealParameter /@ FindAllParameters[{a}]);

IsRealExpression[SARAH`trace[a__]] :=
    TraceEquality[SARAH`trace[a] /. SARAH`Adj -> FSAdj, SARAH`trace[FSAdj[a]]];

IsRealExpression[terms_SARAH`MatMul] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsRealExpression[sum[index_, start_, stop_, expr_]] :=
    IsRealExpression[expr];

IsRealExpression[otherwise_] := False;

HasPhase[particle_] :=
    MemberQ[#[[1]]& /@ SARAH`ParticlePhases, particle];

GetPhase[particle_ /; HasPhase[particle]] :=
    Cases[SARAH`ParticlePhases, {particle, phase_} :> phase][[1]];

GetPhase[_] := Null;

GetTypeFromDimension[sym_, {}] :=
    If[IsRealParameter[sym],
       CConversion`ScalarType[CConversion`realScalarCType],
       CConversion`ScalarType[CConversion`complexScalarCType]
      ];

GetTypeFromDimension[sym_, {0|1}] :=
    GetTypeFromDimension[sym, {}];

GetTypeFromDimension[sym_, {num_?NumberQ}] :=
    Module[{scalarType},
           scalarType = If[IsRealParameter[sym],
                           CConversion`realScalarCType,
                           CConversion`complexScalarCType
                          ];
           CConversion`VectorType[scalarType, num]
          ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ}] :=
    If[IsRealParameter[sym],
       CConversion`MatrixType[CConversion`realScalarCType, num1, num2],
       CConversion`MatrixType[CConversion`complexScalarCType, num1, num2]
      ];

GetRealTypeFromDimension[{}] :=
    CConversion`ScalarType[CConversion`realScalarCType];

GetRealTypeFromDimension[{0}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{1}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{num_?NumberQ}] :=
    CConversion`VectorType[CConversion`realScalarCType, num];

GetRealTypeFromDimension[{num1_?NumberQ, num2_?NumberQ}] :=
    CConversion`MatrixType[CConversion`realScalarCType, num1, num2];

GetType[FlexibleSUSY`M[sym_]] :=
    GetTypeFromDimension[sym, {SARAH`getGen[sym, FlexibleSUSY`FSEigenstates]}];

GetType[sym_] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

GetType[sym_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

GetParameterDimensions[sym_] :=
    Module[{dim},
           dim = SARAH`getDimParameters[sym];
           Switch[dim,
                  {}, {1},
                  {0}, {1},
                  _, dim
                 ]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`ScalarType[_]}] := {};

CreateIndexReplacementRule[{parameter_, CConversion`VectorType[_,_] | CConversion`ArrayType[_,_]}] :=
    Module[{i},
           RuleDelayed @@ Rule[parameter[i_], parameter[i-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`MatrixType[_,_,_]}] :=
    Module[{i,j},
           RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]]
          ];

CreateIndexReplacementRule[parameter_] :=
    Module[{i,j,k,l, dim, rule},
           dim = SARAH`getDimParameters[parameter];
           rule = {};
           Switch[Length[dim],
                  1, rule = RuleDelayed @@ Rule[parameter[i_], parameter[i-1]];,
                  2, rule = RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]];,
                  3, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]];,
                  4, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]];
                 ];
           rule
          ];

CreateIndexReplacementRules[pars_List] :=
    Flatten[CreateIndexReplacementRule /@ pars];

GetGUTNormalization[coupling_Symbol] :=
    Module[{pos, norm},
           pos = Position[SARAH`Gauge, coupling];
           If[pos =!= {},
              norm = SARAH`GUTren[pos[[1,1]]];
              If[NumericQ[norm],
                 Return[norm];
                ];
             ];
           Return[1];
          ];

ApplyGUTNormalization[] :=
    Module[{i, rules = {}, coupling},
           For[i = 1, i <= Length[SARAH`Gauge], i++,
               If[NumericQ[SARAH`GUTren[i]],
                  coupling = SARAH`Gauge[[i,4]];
                  AppendTo[rules, Rule[coupling, coupling SARAH`GUTren[i]]];
                 ];
              ];
           Return[rules];
          ];

CreateSetAssignment[name_, startIndex_, parameterType_, struct_:"pars"] :=
    Block[{},
          Print["Error: CreateSetAssignment: unknown parameter type: ", ToString[parameterType]];
          Quit[1];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType], struct_:"pars"] :=
    Module[{ass = ""},
           ass = name <> " = " <> struct <> "(" <> ToString[startIndex] <> ");\n";
           Return[{ass, 1}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType], struct_:"pars"] :=
    Module[{ass = "", type},
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           ass = name <> " = " <> type <> "(" <> struct <> "(" <> ToString[startIndex] <>
                 "), " <> struct <> "(" <> ToString[startIndex + 1] <> "));\n";
           Return[{ass, 2}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`realScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <> struct <> "(" <>
                     ToString[startIndex + count] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`complexScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++; count+=2,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <>
                     type <> "(" <> struct <>
                     "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                     "(" <> ToString[startIndex + count + 1] <> "));\n";
              ];
           If[2 * rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> ") = " <> struct <> "(" <> ToString[startIndex + count] <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <>
                         ") = " <> type <> "(" <> struct <>
                         "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                         "(" <> ToString[startIndex + count + 1] <> "));\n";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateDisplayAssignment[name_, startIndex_, parameterType_, struct_:"pars"] :=
    Block[{},
          Print["Error: CreateDisplayAssignment: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType], struct_:"pars"] :=
    Module[{ass = ""},
           ass = struct <> "(" <> ToString[startIndex] <> ") = "
                 <> name <> ";\n";
           Return[{ass, 1}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType], struct_:"pars"] :=
    Module[{ass = ""},
           ass = struct <> "(" <> ToString[startIndex] <> ") = Re(" <> name <> ");\n" <>
                 struct <> "(" <> ToString[startIndex + 1] <> ") = Im(" <> name <> ");\n";
           Return[{ass, 2}];
          ];

CreateDisplayAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`realScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> struct <> "(" <> ToString[startIndex + count] <> ") = "
                     <> name <> "(" <> ToString[i] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateDisplayAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`complexScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count+=2,
               ass = ass <> struct <> "(" <> ToString[startIndex + count] <> ") = Re("
                     <> name <> "(" <> ToString[i] <> "));\n";
               ass = ass <> struct <> "(" <> ToString[startIndex + count + 1] <> ") = Im("
                     <> name <> "(" <> ToString[i] <> "));\n";
              ];
           If[2 * rows != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> struct <> "(" <> ToString[startIndex + count] <> ") = "
                          <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                          <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   ass = ass <> struct <> "(" <> ToString[startIndex + count] <> ") = Re("
                         <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> "));\n";
                   ass = ass <> struct <> "(" <> ToString[startIndex + count + 1] <> ") = Im("
                         <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> "));\n";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, parameterType_, struct_:"names"] :=
    Block[{},
          Print["Error: CreateStdVectorNamesAssignment: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType], struct_:"names"] :=
    Module[{ass = ""},
           ass = struct <> "[" <> ToString[startIndex] <> "] = \""
                 <> name <> "\";\n";
           Return[{ass, 1}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType], struct_:"names"] :=
    Module[{ass = ""},
           ass = struct <> "[" <> ToString[startIndex] <> "] = \"Re(" <> name <> ")\";\n" <>
                 struct <> "[" <> ToString[startIndex + 1] <> "] = \"Im(" <> name <> ")\";\n";
           Return[{ass, 2}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`realScalarCType, rows_], struct_:"names"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> struct <> "[" <> ToString[startIndex + count] <> "] = \""
                     <> name <> "(" <> ToString[i] <> ")\";\n";
              ];
           If[rows != count,
              Print["Error: CreateStdVectorNamesAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`complexScalarCType, rows_], struct_:"names"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count+=2,
               ass = ass <> struct <> "[" <> ToString[startIndex + count] <> "] = \"Re("
                     <> name <> "(" <> ToString[i] <> "))\";\n";
               ass = ass <> struct <> "[" <> ToString[startIndex + count + 1] <> "] = \"Im("
                     <> name <> "(" <> ToString[i] <> "))\";\n";
              ];
           If[2 * rows != count,
              Print["Error: CreateStdVectorNamesAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_], struct_:"names"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> struct <> "[" <> ToString[startIndex + count] <> "] = \""
                          <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                          <> ")\";\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateStdVectorAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateStdVectorNamesAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_], struct_:"names"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   ass = ass <> struct <> "[" <> ToString[startIndex + count] <> "] = \"Re("
                         <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> "))\";\n";
                   ass = ass <> struct <> "[" <> ToString[startIndex + count + 1] <> "] = \"Im("
                         <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> "))\";\n";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateStdVectorNamesAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateParameterNamesStr[name_, parameterType_] :=
    Block[{},
          Print["Error: CreateParameterNamesStr: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateParameterNamesStr[name_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]] :=
    "\"" <> CConversion`ToValidCSymbolString[name] <> "\"";

CreateParameterNamesStr[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    "\"Re(" <> CConversion`ToValidCSymbolString[name] <> ")\", " <>
    "\"Im(" <> CConversion`ToValidCSymbolString[name] <> ")\"";

CreateParameterNamesStr[name_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> "\"" <> CConversion`ToValidCSymbolString[name] <>
                     "(" <> ToString[i] <> ")\"";
              ];
           If[rows != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterNamesStr[name_, CConversion`VectorType[CConversion`complexScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count+=2,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> "\"Re(" <> CConversion`ToValidCSymbolString[name] <>
                     "(" <> ToString[i] <> "))\", ";
               ass = ass <> "\"Im(" <> CConversion`ToValidCSymbolString[name] <>
                     "(" <> ToString[i] <> "))\"";
              ];
           If[2 * rows != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[2 * rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterNamesStr[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> "\"" <> CConversion`ToValidCSymbolString[name] <>
                         "(" <> ToString[i] <> "," <> ToString[j] <> ")\"";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterNamesStr[name_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> "\"Re(" <> CConversion`ToValidCSymbolString[name] <>
                         "(" <> ToString[i] <> "," <> ToString[j] <> "))\", ";
                   ass = ass <> "\"Im(" <> CConversion`ToValidCSymbolString[name] <>
                         "(" <> ToString[i] <> "," <> ToString[j] <> "))\"";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[2 * rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, parameterType_] :=
    Block[{},
          Print["Error: CreateParameterEnums: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateParameterEnums[name_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]] :=
    CConversion`ToValidCSymbolString[name];

CreateParameterEnums[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    CConversion`ToValidCSymbolString[Re[name]] <> ", " <>
    CConversion`ToValidCSymbolString[Im[name]];

CreateParameterEnums[name_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> CConversion`ToValidCSymbolString[name[i]];
              ];
           If[rows != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, CConversion`VectorType[CConversion`complexScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count+=2,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> CConversion`ToValidCSymbolString[Re[name[i]]] <> ", "
                         <> CConversion`ToValidCSymbolString[Im[name[i]]];
              ];
           If[2 * rows != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[2 * rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> CConversion`ToValidCSymbolString[name[i,j]];
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> CConversion`ToValidCSymbolString[Re[name[i,j]]] <> ", "
                             <> CConversion`ToValidCSymbolString[Im[name[i,j]]];
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[2 * rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateInputParameterEnum[inputParameters_List] :=
    Module[{i, par, type, name, result = ""},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = Parameters`CreateParameterEnums[par, type];
               If[i > 1, result = result <> ", ";];
               result = result <> name;
              ];
           If[Length[inputParameters] > 0, result = result <> ", ";];
           result = result <> "NUMBER_OF_INPUT_PARAMETERS";
           result = "enum Input_parameters : unsigned {" <>
                    result <> "};\n";
           Return[result];
          ];

CreateInputParameterNames[inputParameters_List] :=
    Module[{i, par, type, name, result = ""},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = Parameters`CreateParameterNamesStr[par, type];
               If[i > 1, result = result <> ", ";];
               result = result <> name;
              ];
           result = "const char* input_parameter_names[NUMBER_OF_INPUT_PARAMETERS] = {" <>
                    result <> "};\n";
           Return[result];
          ];

SetInputParameter[parameter_, value_, wrapper_String, castToType_:None] :=
    Module[{parameterStr, valueStr},
           If[IsInputParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[parameter];
              valueStr = CConversion`RValueToCFormString[value];
              wrapper <> "(" <> parameterStr <> ") = " <> CConversion`CastTo[valueStr,castToType] <> ";\n",
              ""
             ]
          ];

SetPhase[phase_, value_, class_String] :=
    Module[{phaseStr, valueStr},
           If[IsPhase[phase],
              phaseStr = Phases`CreatePhaseName[phase];
              valueStr = CConversion`RValueToCFormString[value];
              class <> "->set_" <> phaseStr <> "(" <> valueStr <> ");\n",
              ""
             ]
          ];

ConcatIndices[indices_List] :=
    Utils`StringJoinWithSeparator[ToString /@ indices,","];

CreateIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    "(" <> ConcatIndices[{indices}] <> ")";

CreateIndices[parameter_] := "";

AppendIfNotEmpty[str_String, sym_String] := If[str == "", "", str <> sym];

SetParameter[Re[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, newValue},
           parameterStr = CConversion`RValueToCFormString[parameter];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(" <> value <> ",Im(MODELPARAMETER(" <> parameterStr <> ")))";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[Im[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, newValue},
           parameterStr = CConversion`RValueToCFormString[parameter];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(Re(MODELPARAMETER(" <> parameterStr <> "))," <> value <> ")";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[parameter_, value_String, class_String, castToType_:None] :=
    Module[{parameterStr, targetType = castToType},
           If[IsModelParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[StripIndices[parameter]];
              (* if the parameter indices, we need to cast to the element type *)
              If[GetIndices[parameter] =!= {} && targetType =!= None,
                 targetType = CConversion`GetScalarElementType[targetType];
                ];
              class <> "->set_" <> parameterStr <> "(" <>
              AppendIfNotEmpty[ConcatIndices[GetIndices[parameter]],","] <>
              CConversion`CastTo[value,targetType] <> ");\n",
              ""
             ]
          ];

SetParameter[parameter_, value_, class_String] :=
    SetParameter[parameter, CConversion`RValueToCFormString[value], class, GetType[parameter]]

SetParameter[Re[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".real(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[Im[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".imag(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[parameter_, value_, castToType_:None] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> CreateIndices[parameter] <> " = " <>
    CConversion`CastTo[CConversion`RValueToCFormString[value],castToType] <> ";\n";

SaveParameterLocally[parameters_List, prefix_String, caller_String] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> SaveParameterLocally[parameters[[i]], prefix, caller];
              ];
           Return[result];
          ];

SaveParameterLocally[parameter_, prefix_String, caller_String] :=
    Module[{ parStr, parStrSym },
           parStr = CConversion`RValueToCFormString[parameter];
           parStrSym = CConversion`ToValidCSymbolString[parameter];
           "const auto " <> prefix <> parStrSym <> " = " <>
           If[caller != "", caller <> "(" <> parStr <> ")", parStr] <> ";\n"
          ];

RestoreParameter[parameters_List, prefix_String, modelPtr_String] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> RestoreParameter[parameters[[i]], prefix, modelPtr];
              ];
           Return[result];
          ];

RestoreParameter[parameter_, prefix_String, modelPtr_String] :=
    Module[{ parStr, parStrSym },
           parStr = CConversion`RValueToCFormString[parameter];
           parStrSym = CConversion`ToValidCSymbolString[parameter];
           If[modelPtr != "",
              SetParameter[parameter, prefix <> parStrSym, modelPtr],
              parStr <> " = " <> prefix <> parStrSym <> ";\n"]
          ];

RemoveProtectedHeads[expr_] :=
    expr /. { FlexibleSUSY`LowEnergyConstant[__] -> FlexibleSUSY`LowEnergyConstant[],
              FlexibleSUSY`Pole[__]  -> FlexibleSUSY`Pole[] };

DefineLocalConstCopy[parameter_, macro_String, prefix_String:""] :=
    "const auto " <> prefix <> ToValidCSymbolString[parameter] <> " = " <>
    macro <> "(" <> ToValidCSymbolString[parameter] <> ");\n";

PrivateCallLoopMassFunction[FlexibleSUSY`M[particle_Symbol]] :=
    "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole();\n";

CalculateLocalPoleMasses[parameter_] :=
    "MODEL->" <> PrivateCallLoopMassFunction[parameter];

CreateLocalConstRefs[expr_] :=
    Module[{result = "", symbols, inputSymbols, modelPars, outputPars,
            poleMasses, phases, depNum},
           symbols = FindAllParameters[expr];
           poleMasses = {
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_]]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}],
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_[__]]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
                        };
           symbols = DeleteDuplicates[Flatten[symbols]];
           poleMasses = DeleteDuplicates[Flatten[poleMasses]];
           inputSymbols = DeleteDuplicates[Select[symbols, (MemberQ[allInputParameters,#])&]];
           modelPars    = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           outputPars   = DeleteDuplicates[Select[symbols, (MemberQ[allOutputParameters,#])&]];
           phases       = DeleteDuplicates[Select[symbols, (MemberQ[Phases`GetArg /@ allPhases,#])&]];
           depNum       = DeleteDuplicates[Select[symbols, (MemberQ[GetDependenceNumSymbols[],#])&]];
           (result = result <> DefineLocalConstCopy[#,"INPUTPARAMETER"])& /@ inputSymbols;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ modelPars;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ outputPars;
           (result = result <> DefineLocalConstCopy[#,"PHASE"         ])& /@ phases;
           (result = result <> DefineLocalConstCopy[#,"DERIVEDPARAMETER"])& /@ depNum;
           (result = result <> CalculateLocalPoleMasses[#])& /@ poleMasses;
           Return[result];
          ];

CreateLocalConstRefsForPhysicalParameters[expr_] :=
    Module[{result = "", symbols, outputPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, {0,Infinity}],
                       Cases[compactExpr, a_[__] /; MemberQ[allOutputParameters,a] :> a, {0,Infinity}],
                       Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]], {0,Infinity}],
                       Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
                     };
           symbols = DeleteDuplicates[Flatten[symbols]];
           outputPars = DeleteDuplicates[Select[symbols, (MemberQ[allOutputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#,"PHYSICAL"])& /@ outputPars;
           Return[result];
          ];

CreateLocalConstRefsForBetas[expr_] :=
    Module[{result = "", symbols, modelPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, {0,Infinity}],
                       Cases[compactExpr, a_[__] /; MemberQ[allModelParameters,a] :> a, {0,Infinity}] };
           symbols = DeleteDuplicates[Flatten[symbols]];
           modelPars = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#, "BETAPARAMETER", "beta_"])& /@ modelPars;
           Return[result];
          ];

CreateLocalConstRefsForInputParameters[expr_, head_String:"INPUT"] :=
    Module[{result = "", symbols, inputPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = Cases[compactExpr, _Symbol, {0,Infinity}];
           symbols = DeleteDuplicates[Flatten[symbols]];
           inputPars = DeleteDuplicates[Select[symbols, (MemberQ[allInputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#, head])& /@ inputPars;
           Return[result];
          ];

SetInputParameterTo[Re[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".real(" <> value <> ");";

SetInputParameterTo[Im[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".imag(" <> value <> ");";

SetInputParameterTo[par_, value_String] :=
    CConversion`ToValidCSymbolString[par] <> " = " <> value <> ";";

CreateCaseFromTuple[{key_?NumberQ, parameter_}] :=
    "case " <> ToString[key] <> ": input." <>
    SetInputParameterTo[parameter, "value"] <> " break;\n";

CreateCaseFromTuple[expr_] :=
    Block[{},
          Print["Error: not a valid {key,parameter} tuple: ", expr];
          ""
         ];

FillInputParametersFromTuples[minpar_List, blockName_String] :=
    Module[{result = ""},
           (result = result <> CreateCaseFromTuple[#])& /@ minpar;
           result = "switch (key) {\n" <> result <>
                    "default: WARNING(\"Unrecognized entry in block " <>
                    blockName <> ": \" << key); break;\n}\n";
           Return[result];
          ];

IncreaseIndex[ind_Integer, num_Integer] := ind + num;
IncreaseIndex[ind_, _]     := ind;
IncreaseIndices[a_[{ind__}], num_Integer] := a[IncreaseIndex[#,num]& /@ {ind}];
IncreaseIndices[a_[ind__], num_Integer] := a[Sequence @@ (IncreaseIndex[#,num]& /@ {ind})];
IncreaseIndices[a_, _]     := a;
IncreaseIndices[SARAH`Delta[a_, b_], num_Integer] :=
    CConversion`FSKroneckerDelta[IncreaseIndex[a,num], IncreaseIndex[b,num]];

IncreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, 1];

IncreaseIndexLiterals[expr_, num_Integer] :=
    IncreaseIndexLiterals[expr, num, Join[allInputParameters, allModelParameters,
                                          allOutputParameters]];

IncreaseIndexLiterals[expr_, num_Integer, heads_List] :=
    Module[{indexedSymbols, rules, decrExpr, allHeads},
           allHeads = Join[heads /. FlexibleSUSY`M -> Identity, {SARAH`Delta, SARAH`ThetaStep}];
           indexedSymbols = Cases[{expr}, s_[__] /; MemberQ[allHeads, s], Infinity];
           rules = Rule[#, IncreaseIndices[#,num]] & /@ indexedSymbols;
           decrExpr = expr /. rules;
           Return[decrExpr]
          ];

DecreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, -1];

DecreaseIndexLiterals[expr_, heads_List] :=
    IncreaseIndexLiterals[expr, -1, heads];

DecreaseSumIdices[expr_] :=
    expr //. SARAH`sum[idx_, start_, stop_, exp_] :> CConversion`IndexSum[idx, start - 1, stop - 1, exp];

GetEffectiveMu[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMu],
              Print["Error: effective Mu parameter not defined!"];
              Print["   Please set EffectiveMu to the expression of the"];
              Print["   effective Mu parameter."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMu
          ];

GetEffectiveMASqr[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMASqr],
              Print["Error: effective CP-odd Higgs mass not defined!"];
              Print["   Please set EffectiveMASqr to the expression of the"];
              Print["   effective CP-odd Higgs mass."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMASqr
          ];

GetParameterFromDescription[description_String] :=
    Module[{parameter},
           parameter =Cases[SARAH`ParameterDefinitions,
                            {parameter_,
                             {___, SARAH`Description -> description, ___}} :>
                            parameter];
           If[Length[parameter] == 0,
              Print["Error: Parameter with description \"", description,
                    "\" not found."];
              Return[Null];
             ];
           If[Length[parameter] > 1,
              Print["Warning: Parameter with description \"", description,
                    "\" not unique."];
             ];
           parameter[[1]]
          ];

GetParticleFromDescription[description_String, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{particle},
           particle =Cases[SARAH`ParticleDefinitions[eigenstates],
                            {particle_,
                             {___, SARAH`Description -> description, ___}} :>
                            particle];
           If[Length[particle] == 0,
              DebugPrint["Note: Particle with description \"", description,
                         "\" not found."];
              Return[Null];
             ];
           If[Length[particle] > 1,
              Print["Warning: Particle with description \"", description,
                    "\" not unique."];
             ];
           particle[[1]]
          ];

GetParticleFromDescription[multipletName_String, splitNames_List] :=
    Module[{result},
           result = GetParticleFromDescription[multipletName];
           If[result =!= Null, Return[{result}]];
           GetParticleFromDescription /@ splitNames
          ];

GetPDGCodesForParticle[particle_] :=
    Module[{pdgList},
            pdgList = SARAH`getPDGList[particle];
            If[pdgList === None,
               pdgList = {};
              ];
           pdgList
          ];

NumberOfIndependentEntriesOfSymmetricMatrix[n_] := (n^2 + n) / 2;

AppendGenerationIndices[expr_List] :=
    AppendGenerationIndices /@ expr;

AppendGenerationIndices[expr_Symbol] :=
    Switch[SARAH`getDimParameters[expr],
           {}                          , expr,
           {0}                         , expr,
           {1}                         , expr,
           {idx_}                      , expr[SARAH`gt1],
           {idx1_, idx2_}              , expr[SARAH`gt1, SARAH`gt2],
           {idx1_, idx2_, idx3_}       , expr[SARAH`gt1, SARAH`gt2, SARAH`gt3],
           {idx1_, idx2_, idx3_, idx4_}, expr[SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4]
          ];

AppendGenerationIndices[expr_] := expr;

(*
 * Expands a list of expressions of the form
 *
 *   { 1 + A[SARAH`gt1] }
 *
 * to
 *
 *   { 1 + A[1], 1 + A[2], 1 + A[3] }
 *
 * where the indices SARAH`gt1 and SARAH`gt2 are assumed to run from 1
 * to their maximum value.
 *)
ExpandExpressions[eqs_List] :=
    Join[Flatten[ExpandExpressions /@ eqs]];

GetIdxRange[{idx_,par_}] :=
    Module[{dim = GetParameterDimensions[par]},
           Switch[dim,
                  {}  , {idx,1,1},
                  {0} , {idx,1,1},
                  {n_?NumberQ}, {idx,1,dim[[1]]},
                  {__}, {idx,1,1}
                 ]
          ];

ExpandExpressions[eq_] :=
    Module[{par, indexSymbols, indexRanges, indices = {SARAH`gt1, SARAH`gt2}},
           indexSymbols = DeleteDuplicates[
               Join @@ (Cases[eq, par_[idx_ /; !FreeQ[idx,#]] /;
                                  MemberQ[Join[GetModelParameters[], GetOutputParameters[]], par] :> {#,par},
                              {0,Infinity}]& /@ indices)];
           indexRanges  = DeleteDuplicates[GetIdxRange /@ indexSymbols];
           If[indexRanges === {},
              eq,
              Table[eq, Evaluate[Sequence @@ indexRanges]]
             ]
          ];

(* removes indices from model Parameter, taking SARAH's L, B, T, Q
   into account.  *)
StripIndices[par_[idx___] /; And @@ (NumberQ /@ {idx})] := par;

StripIndices[par_[idx___] /; MemberQ[Join[allModelParameters,allOutputParameters],par]] := par;

StripIndices[par_] := par;

ExtractParametersFromSARAHBetaLists[beta_List] :=
    StripIndices[#[[1]]]& /@ beta;

GetModelParametersWithMassDimension[dim_?IntegerQ] :=
    Module[{dimPars},
           Switch[dim,
                  0, dimPars = Join[SARAH`BetaGauge, SARAH`BetaLijkl, SARAH`BetaYijk, SARAH`BetaQijkl];,
                  1, dimPars = Join[SARAH`BetaMuij, SARAH`BetaTijk, SARAH`BetaMi, SARAH`BetaDGi, SARAH`BetaVEV];,
                  2, dimPars = Join[SARAH`BetaLSi, SARAH`BetaBij, SARAH`Betam2ij];,
                  3, dimPars = Join[SARAH`BetaLi];,
                  _,
                  Print["Error: GetModelParametersWithMassDimension: ", dim,
                        " is not a valid dimension"];
                  Return[{}];
                 ];
           ExtractParametersFromSARAHBetaLists[dimPars]
          ];

AreLinearDependent[{eq1_, eq2_}, parameters_List] :=
    Module[{frac = Simplify[eq1/eq2 /. FlexibleSUSY`tadpole[_] -> 0],
            pars},
           (* ignore parameter heads Re[], Im[], Abs[], Phase[] *)
           pars = parameters /. { Re[p_] :> p, Im[p_] :> p,
                                  Abs[p_] :> p, FlexibleSUSY`Phase[p_] :> p };
           And @@ (FreeQ[frac,#]& /@ pars)
          ];

FilterOutLinearDependentEqs[{}, _List] := {};

FilterOutLinearDependentEqs[{eq_}, _List] := {eq};

FilterOutLinearDependentEqs[{eq_, rest__}, parameters_List] :=
    If[Or @@ (AreLinearDependent[#,parameters]& /@ ({eq,#}& /@ {rest})),
       (* leave out eq and check rest *)
       FilterOutLinearDependentEqs[{rest}, parameters],
       (* keep eq and check rest *)
       {eq, Sequence @@ FilterOutLinearDependentEqs[{rest}, parameters]}
      ];

FilterOutIndependentEqs[eqs_List, pars_List] :=
    DeleteDuplicates @ Flatten @ Join[FilterOutIndependentEqs[eqs,#]& /@ pars];

FilterOutIndependentEqs[eqs_List, p_] :=
    Select[eqs, (!FreeQ[#,p])&];

GetThirdGeneration[par_] :=
    Which[IsScalar[par], par,
          IsMatrix[par], par[2,2],
          True, Print["Warning: GetThirdGeneration[",par,"]: unknown type"]; par
         ];

GetDependenceNumSymbols[] :=
    DeleteDuplicates @ Flatten @
    Join[{SARAH`Weinberg},
         Cases[SARAH`ParameterDefinitions,
               {parameter_ /; !MemberQ[Parameters`GetModelParameters[], parameter] &&
                parameter =!= SARAH`Weinberg && parameter =!= SARAH`electricCharge,
                {___, SARAH`DependenceNum -> value:Except[None], ___}} :> parameter]
        ];

CreateInputParameterArrayGetter[inputParameters_List] :=
    Module[{get = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               get = get <> assignment;
               paramCount += nAssignments;
              ];
           get = "Eigen::ArrayXd pars(" <> ToString[paramCount] <> ");\n\n" <>
                 get <> "\n" <>
                 "return pars;";
           Return[get];
          ];

CreateInputParameterArraySetter[inputParameters_List] :=
    Module[{set = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           Return[set];
          ];

End[];

EndPackage[];
