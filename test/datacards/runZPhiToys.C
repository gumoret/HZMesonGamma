{
  gROOT->ProcessLine(".L ../../../HiggsAnalysis/CombinedLimit/test/plotParametersFromToys.C");
  plotParametersFromToys("./fitDiagnosticsPhiZ_toys.root","./fitDiagnosticsPhiZ_total.root","datacard_STAT_Phi_Z_bdt0_2024.root","r<0")
}
