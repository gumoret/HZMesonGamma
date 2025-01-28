import ROOT
import argparse

# Suppress the opening of many Canvas's
ROOT.gROOT.SetBatch(False)  # Set to False to see the Canvas

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('outputfile_name', help='Provide output file name (to read histograms)')

args = p.parse_args()
output_filename = args.outputfile_name

# Open the output ROOT file
fInput = ROOT.TFile(output_filename, "READ")

# Retrieve the histograms
h_MesonMassTrkTrk = fInput.Get("h_MesonMassTrkTrk")
h_MesonTreeMass = fInput.Get("h_MesonTreeMass")

# Set styles for the histograms
h_MesonMassTrkTrk.SetLineColor(ROOT.kRed)
h_MesonMassTrkTrk.SetLineWidth(2)
h_MesonTreeMass.SetLineColor(ROOT.kBlue)
h_MesonTreeMass.SetLineWidth(2)

# Create a Canvas
c1 = ROOT.TCanvas("c1", "Meson Mass Distributions", 800, 600)

# Draw the first histogram
h_MesonMassTrkTrk.Draw("HIST")  # Use HIST to draw as a histogram
h_MesonTreeMass.Draw("HIST SAME")  # Use SAME to overlay the second histogram

# Add a legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)  # Adjust the position as needed
legend.AddEntry(h_MesonMassTrkTrk, "Meson Mass (TrkTrk)", "l")
legend.AddEntry(h_MesonTreeMass, "Meson Mass (Tree)", "l")
legend.Draw()

# Add labels and title
h_MesonMassTrkTrk.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
h_MesonMassTrkTrk.GetYaxis().SetTitle("Events")
h_MesonMassTrkTrk.SetTitle("Comparison of Meson Mass Distributions")

# Save the Canvas
c1.SaveAs("meson_mass_comparison.png")

# Show the Canvas
c1.Update()

# Close the ROOT file
fInput.Close()
