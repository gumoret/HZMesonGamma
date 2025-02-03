import ROOT
import argparse

# Suppress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  # Set to False to see the Canvas

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('outputfile_name', help='Provide output file name (to read histograms)')

args = p.parse_args()
output_filename = args.outputfile_name

# Open the output ROOT file
fInput = ROOT.TFile(output_filename, "READ")

# Retrieve the histograms
hist_pairs = [
    ("h_MesonMassTrkTrk", "h_MesonTreeMass", "h_MesonTreeMassKin", "Meson Mass [GeV/c^2]", "Meson Mass Comparison"),
    ("h_MesonPt", "h_MesonTreePt", None, "p_{T} [GeV/c]", "Meson pT Comparison"),
    ("h_MesonEta", "h_MesonTreeEta", None, "#eta", "Meson Eta Comparison"),
    ("h_MesonPhi", "h_MesonTreePhi", None, "#phi [rad]", "Meson Phi Comparison"),
]

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen]

for hist_names in hist_pairs:
    h1_name, h2_name, h3_name, x_title, plot_title = hist_names
    h1 = fInput.Get(h1_name)
    h2 = fInput.Get(h2_name)
    h3 = fInput.Get(h3_name) if h3_name else None
    
    h1.SetLineColor(colors[0])
    h1.SetLineWidth(2)
    h2.SetLineColor(colors[1])
    h2.SetLineWidth(2)
    if h3:
        h3.SetLineColor(colors[2])
        h3.SetLineWidth(2)
    
    c = ROOT.TCanvas(f"c_{h1_name}", plot_title, 800, 600)
    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    if h3:
        h3.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
    legend.AddEntry(h1, h1_name.replace("h_", ""), "l")
    legend.AddEntry(h2, h2_name.replace("h_", ""), "l")
    if h3:
        legend.AddEntry(h3, h3_name.replace("h_", ""), "l")
    legend.Draw()
    
    h1.GetXaxis().SetTitle(x_title)
    h1.GetYaxis().SetTitle("Events")
    h1.SetTitle(plot_title)
    
    c.SaveAs(f"/eos/user/g/gumoret/www/HZMesonGamma/var_comparisons/rho/{h1_name}_comparison.png")
    c.SaveAs(f"/eos/user/g/gumoret/www/HZMesonGamma/var_comparisons/rho/{h1_name}_comparison.pdf")

    c.Update()

# Close the ROOT file
fInput.Close()