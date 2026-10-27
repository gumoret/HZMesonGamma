{
  gROOT->ProcessLine(".L ../../../HiggsAnalysis/CombinedLimit/test/plotParametersFromToys.C");
  plotParametersFromToys("./fitDiagnosticsDstarGFcat_toys.root","./fitDiagnosticsDstarGFcat_total.root","datacard_STAT_Dstar_GFcat_Comb_2024.root","r<0")
}
