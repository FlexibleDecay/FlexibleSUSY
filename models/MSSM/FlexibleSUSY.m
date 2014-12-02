
FSModelName = "MSSM";
FSEigenstates = SARAH`EWSB;

(* CMSSM input parameters *)

MINPAR = {
    {3, TanBeta},
    {4, Sign[\[Mu]]}
};

EXTPAR = {
    {0, Qin},
    {21, mHd2IN},
    {22, mHu2IN}
};

FSExtraInputParameters = {
    {Aeij, AeijIN, {3,3}}, (* 3x3 matrix *)
    {Adij, AdijIN, {3,3}}, (* 3x3 matrix *)
    {Auij, AuijIN, {3,3}}  (* 3x3 matrix *)
};

EWSBOutputParameters = { B[\[Mu]], \[Mu] };

SUSYScale = Sqrt[Product[M[Su[i]]^(Abs[ZU[i,3]]^2 + Abs[ZU[i,6]]^2), {i,6}]];

SUSYScaleFirstGuess = SM[MZ];

SUSYScaleInput = {};

HighScale == Qin;

HighScaleFirstGuess = 2.0 10^16;

HighScaleInput = {
   {T[Ye], Aeij*Ye},
   {T[Yd], Adij*Ye},
   {T[Yu], Auij*Ye},
   {mHd2, mHd2IN},
   {mHu2, mHu2IN},
   {mq2, LHInput[mq2]},
   {ml2, LHInput[ml2]},
   {md2, LHInput[md2]},
   {mu2, LHInput[mu2]},
   {me2, LHInput[me2]},
   {MassB, LHInput[MassB]},
   {MassWB,LHInput[MassWB]},
   {MassG, LHInput[MassG]}
};

LowScale = SM[MZ];

LowScaleFirstGuess = SM[MZ];

LowScaleInput = {
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic},
   {vd, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Cos[ArcTan[TanBeta]]},
   {vu, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Sin[ArcTan[TanBeta]]}
};

InitialGuessAtLowScale = {
   {vd, SM[vev] Cos[ArcTan[TanBeta]]},
   {vu, SM[vev] Sin[ArcTan[TanBeta]]},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

InitialGuessAtHighScale = {
   {\[Mu]   , 1.0},
   {B[\[Mu]], 0.0}
};

UseHiggs2LoopMSSM = True;
EffectiveMu = \[Mu];

PotentialLSPParticles = { Chi, Sv, Su, Sd, Se, Cha, Glu };

DefaultPoleMassPrecision = MediumPrecision;
HighPoleMassPrecision    = {hh, Ah, Hpm};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

ExtraSLHAOutputBlocks = {
   {ALPHA, {{ArcSin[Pole[ZH[2,2]]]}}},
   {HMIX , {{1, \[Mu]},
            {2, vu / vd},
            {3, Sqrt[vu^2 + vd^2]},
            {4, M[Ah[2]]^2},
            {101, B[\[Mu]]},
            {102, vd},
            {103, vu} } },
   {Au,    {{1, 1, T[Yu][1,1] / Yu[1,1]},
            {2, 2, T[Yu][2,2] / Yu[2,2]},
            {3, 3, T[Yu][3,3] / Yu[3,3]} } },
   {Ad,    {{1, 1, T[Yd][1,1] / Yd[1,1]},
            {2, 2, T[Yd][2,2] / Yd[2,2]},
            {3, 3, T[Yd][3,3] / Yd[3,3]} } },
   {Ae,    {{1, 1, T[Ye][1,1] / Ye[1,1]},
            {2, 2, T[Ye][2,2] / Ye[2,2]},
            {3, 3, T[Ye][3,3] / Ye[3,3]} } },
   {MSOFT, {{1, MassB},
            {2, MassWB},
            {3, MassG},
            {21, mHd2},
            {22, mHu2},
            {31, Sqrt[ml2[1,1]]},
            {32, Sqrt[ml2[2,2]]},
            {33, Sqrt[ml2[3,3]]},
            {34, Sqrt[me2[1,1]]},
            {35, Sqrt[me2[2,2]]},
            {36, Sqrt[me2[3,3]]},
            {41, Sqrt[mq2[1,1]]},
            {42, Sqrt[mq2[2,2]]},
            {43, Sqrt[mq2[3,3]]},
            {44, Sqrt[mu2[1,1]]},
            {45, Sqrt[mu2[2,2]]},
            {46, Sqrt[mu2[3,3]]},
            {47, Sqrt[md2[1,1]]},
            {48, Sqrt[md2[2,2]]},
            {49, Sqrt[md2[3,3]]} } }
};
