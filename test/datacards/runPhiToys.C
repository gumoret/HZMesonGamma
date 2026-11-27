{
  gROOT->ProcessLine(".L ../../../HiggsAnalysis/CombinedLimit/test/plotParametersFromToys.C");
  plotParametersFromToys("./fitDiagnosticsPhiGFcat_toys.root","./fitDiagnosticsPhiGFcat_total.root","datacard_STAT_Phi_GFcat_Comb_2024.root","r<0")
}
