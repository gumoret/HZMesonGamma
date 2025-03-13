import ROOT
import argparse

# Suppress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)  # Set to False to see the Canvas

isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma
isKAnalysis   = False # for Z -> K* Gamma
isD0Analysis  = False # for Z -> D0* Gamma

isHAnalysis   = False
isZAnalysis   = False


# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select plot options')
p.add_argument('boson_channel', help='type <<H>> or <<Z>>')
p.add_argument('meson_channel', help='type <<rho>> or <<phi>> or <<K>> or <<D0>>')
p.add_argument('rootfile_name', help='Type rootfile name to plot')


args = p.parse_args()

if args.boson_channel == "H": isHAnalysis = True 
elif args.boson_channel == "Z": isZAnalysis = True 
if args.meson_channel == "phi": isPhiAnalysis = True 
if args.meson_channel == "rho": isRhoAnalysis = True 
if args.meson_channel == "K": isKAnalysis = True 
if args.meson_channel == "D0": isD0Analysis = True 

# Open the output ROOT file
fInput = ROOT.TFile(args.rootfile_name, "READ")

h_EventsOR = fInput.Get("nEventsOR")
h_Events50 = fInput.Get("nEvents50")
h_Events35 = fInput.Get("nEvents35")


bin1contentOR  = h_EventsOR.GetBinContent(1)
bin2contentOR  = h_EventsOR.GetBinContent(2)
bin3contentOR  = h_EventsOR.GetBinContent(3)
bin4contentOR  = h_EventsOR.GetBinContent(4)
bin5contentOR  = h_EventsOR.GetBinContent(5)

bin1content50  = h_Events50.GetBinContent(1)
bin2content50  = h_Events50.GetBinContent(2)
bin3content50  = h_Events50.GetBinContent(3)
bin4content50  = h_Events50.GetBinContent(4)
bin5content50  = h_Events50.GetBinContent(5)

bin1content35  = h_Events35.GetBinContent(1)
bin2content35  = h_Events35.GetBinContent(2)
bin3content35  = h_Events35.GetBinContent(3)
bin4content35  = h_Events35.GetBinContent(4)
bin5content35  = h_Events35.GetBinContent(5)

nSignal      = bin1contentOR
scale_factor = 100./nSignal

h_EventsOR.Reset()
h_Events50.Reset()
h_Events35.Reset()

h_EventsOR.SetBinContent(1, bin1contentOR * scale_factor)  
h_Events50.SetBinContent(1, bin1content50 * scale_factor)  
h_Events35.SetBinContent(1, bin1content35 * scale_factor)

h_EventsOR.SetBinContent(2, bin2contentOR * scale_factor)  
h_Events50.SetBinContent(2, bin2content50 * scale_factor)  
h_Events35.SetBinContent(2, bin2content35 * scale_factor)  

h_EventsOR.SetBinContent(3, bin3contentOR * scale_factor)  
h_Events50.SetBinContent(3, bin3content50 * scale_factor)  
h_Events35.SetBinContent(3, bin3content35 * scale_factor)  

h_EventsOR.SetBinContent(4, bin4contentOR * scale_factor)  
h_Events50.SetBinContent(4, bin4content50 * scale_factor)  
h_Events35.SetBinContent(4, bin4content35 * scale_factor)  

h_EventsOR.SetBinContent(5, bin5contentOR * scale_factor)  
h_Events50.SetBinContent(5, bin5content50 * scale_factor)  
h_Events35.SetBinContent(5, bin5content35 * scale_factor)  

# Labels
h_EventsOR.GetXaxis().SetBinLabel(1, "Processed")
h_EventsOR.GetXaxis().SetBinLabel(2, "Triggered")
h_EventsOR.GetXaxis().SetBinLabel(3, "Photon")
h_EventsOR.GetXaxis().SetBinLabel(4, "Best meson")
h_EventsOR.GetXaxis().SetBinLabel(5, "Trks pT")

h_EventsOR.GetXaxis().SetTitle("")
h_EventsOR.GetYaxis().SetTitle("#epsilon (%)")


h_Events35.SetFillColor(ROOT.kBlue)
h_Events50.SetFillColor(ROOT.kRed)
h_EventsOR.SetFillColor(ROOT.kGreen)

h_Events35.SetFillStyle(3007)
h_Events50.SetFillStyle(3012)
h_EventsOR.SetFillStyle(3013)

# Canvas and plot
c = ROOT.TCanvas("c", "Trigger Comparison", 800, 600)
c.cd()

h_EventsOR.Draw("HIST TEXT0")
h_Events50.Draw("HIST SAME TEXT0")
h_Events35.Draw("HIST SAME TEXT0")


ROOT.gStyle.SetPaintTextFormat("4.2f ")
ROOT.gStyle.SetOptStat(0)
h_EventsOR.SetMarkerSize(1.4)
h_Events50.SetMarkerSize(1.4)
h_Events35.SetMarkerSize(1.4)
h_EventsOR.GetYaxis().SetRangeUser(0.,105)

# Legenda
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_Events35, "Photon35_TwoProngs", "f")
legend.AddEntry(h_Events50, "Photon50EB_TightID", "f")
legend.AddEntry(h_EventsOR, "Trigger OR", "f")
legend.Draw()


if isPhiAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HPhi/Signal/"
if isPhiAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZPhi_MadGraph/Signal/"

if isRhoAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HRho/Signal/"
if isRhoAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZRho_MadGraph/Signal/"

if isKAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HK/Signal/"
if isKAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZK/Signal/"

c.SaveAs(output_dir + "Triggers_comparison.pdf")
c.SaveAs(output_dir + "Triggers_comparison.png")