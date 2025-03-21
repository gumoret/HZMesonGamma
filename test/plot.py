import ROOT
import math
import copy
import sys
import argparse
import tdrstyle, CMS_lumi
from ROOT import gROOT


isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma
isKAnalysis   = False # for Z -> K* Gamma
isD0Analysis  = False # for Z -> D0* Gamma

isHAnalysis   = False
isZAnalysis   = False

plotOnlyData  = False

isWideRange   = False  #for plotting in the right directory

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   
# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select plot options')
p.add_argument('boson_channel', help='type <<H>> or <<Z>>')
p.add_argument('meson_channel', help='type <<rho>> or <<phi>> or <<K*>> or <<D0*>>')
p.add_argument('plot_only_data', help='type <<signal>> or <<data>>')
p.add_argument('SR_or_CR', help='type <<SR>> or <<CR>>')
p.add_argument('rootfile_name', help='Type rootfile name')


args = p.parse_args()

if args.boson_channel == "H": isHAnalysis = True 
elif args.boson_channel == "Z": isZAnalysis = True 
if args.meson_channel == "phi": isPhiAnalysis = True 
if args.meson_channel == "rho": isRhoAnalysis = True 
if args.meson_channel == "K": isKAnalysis = True 
if args.meson_channel == "D0": isD0Analysis = True 

if args.plot_only_data == "data": plotOnlyData = True
if args.plot_only_data == "signal": plotOnlyData = False

signal_magnify = 1

region_name = args.SR_or_CR #CR or SR

inputfile = args.rootfile_name


#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.9
CMS_lumi.cmsTextSize = 1.2
CMS_lumi.lumi_13TeV = ""

hstack  = dict()
hsignal = dict()
canvas  = dict()
histo_container = [] #just for memory management

#Get the list of histograms
list_histos = []
signalfile = ROOT.TFile(args.rootfile_name)
keylist = signalfile.GetListOfKeys()
key = ROOT.TKey()

for key in keylist :
    obj_class = ROOT.gROOT.GetClass(key.GetClassName())
    if not obj_class.InheritsFrom("TH1") :
        continue
    if not (key.ReadObj().GetName() == "h_efficiency" or key.ReadObj().GetName() == "h_cutOverflow"): #h_efficiency and h_cutOverflow is a plot plotted in other way   
        list_histos.append( key.ReadObj().GetName() )

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")


fileIn = ROOT.TFile(inputfile)

for histo_name in list_histos:
    histo = fileIn.Get(histo_name)

    print("histo_name = ",histo_name)

    #set stats box
    ROOT.gStyle.SetOptStat(111111) 
    histo.SetStats(True)

    # Set to 0 the bins containing negative values, due to negative weights
    hsize = histo.GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
    for bin in range(1,hsize+1): # The +1 is in order to get the last bin
        bincontent = histo.GetBinContent(bin)
        if bincontent < 0.:
            histo.SetBinContent(bin,0.)

    histo_container.append(copy.copy(histo))
    
    if not histo_name == "h_nMuons" and not histo_name == "h_nPhotons"  and not histo_name == "h_nElectrons" and not histo_name == "h_photonWP90" and not histo_name == "h_mesonMass":
        print(histo_name)
        
    if histo_name == "h_mesonMass": histo_container[-1].Rebin(2)

    histo_container[-1].SetLineStyle(1)   #continue line (2 for dashed)

    if isPhiAnalysis: histo_container[-1].SetLineColor(9)   #blue
    elif isRhoAnalysis: histo_container[-1].SetLineColor(46)   #red
    elif isKAnalysis: histo_container[-1].SetLineColor(8) #green

    histo_container[-1].SetLineWidth(4)   #kind of thick
    #histo_container[-1].Scale(1./histo_container[-1].GetEntries()) #normalize to 1
    hsignal[histo_name] = histo_container[-1]
    hstack[histo_name].Add(histo_container[-1])

fileIn.Close()

for histo_name in list_histos:

    canvas[histo_name] = ROOT.TCanvas("Canvas_" + histo_name,"",200,106,600,600)
    canvas[histo_name].cd()

    hstack[histo_name].SetTitle("")
    hsignal[histo_name].SetTitle("")


    #if plotOnlyData :

    hstack[histo_name].Draw("histo")
    hstack[histo_name].GetYaxis().SetTitleSize(0.04)
    hstack[histo_name].GetXaxis().SetTitleSize(0.045)
    hstack[histo_name].GetYaxis().SetTitleOffset(1.25)
    hstack[histo_name].GetXaxis().SetTitleOffset(1.15)
    hstack[histo_name].GetYaxis().SetTitle("Events")
    hstack[histo_name].GetYaxis().SetMaxDigits(3)
    hstack[histo_name].GetXaxis().SetLabelSize(0.04)
    hstack[histo_name].GetYaxis().SetLabelSize(0.04)
        
    hstack[histo_name].SetMaximum(1.5 * hsignal[histo_name].GetMaximum())

    #Legend ----------------------------------------
    leg1 = ROOT.TLegend(0.65, 0.68, 0.87, 0.91) #right positioning
    leg1.SetHeader(" ")
    leg1.SetNColumns(1)
    leg1.SetFillColorAlpha(0,0.)
    leg1.SetBorderSize(0)
    leg1.SetLineColor(1)
    leg1.SetLineStyle(1)
    leg1.SetLineWidth(1)
    leg1.SetFillStyle(1001)

    if histo_name == "h_bosonMass":
        #hstack[histo_name].Rebin(2)            
        hstack[histo_name].GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
        #hstack[histo_name].GetXaxis().SetLimits(115.,135.)

    
    if histo_name == "h_mesonMass" :
        hstack[histo_name].GetXaxis().SetTitle("m_{meson} [GeV]")
        if isPhiAnalysis:
            #hstack[histo_name].GetXaxis().SetLimits(1.00,1.042)
            leftLine  = ROOT.TLine(1.008,0.,1.008,hsignal[histo_name].GetMaximum()*1.1)
            rightLine = ROOT.TLine(1.032,0.,1.032,hsignal[histo_name].GetMaximum()*1.1)
            leftLine.SetLineColor(4)
            leftLine.SetLineStyle(2)
            leftLine.SetLineWidth(3)
            rightLine.SetLineColor(4)
            rightLine.SetLineStyle(2)
            rightLine.SetLineWidth(3)
            leftLine.Draw()
            rightLine.Draw()

        if isRhoAnalysis:
            #hstack[histo_name].GetXaxis().SetLimits(0.5,1.)
            leftLine  = ROOT.TLine(0.62,0.,0.62,hsignal[histo_name].GetMaximum()*1.1)
            rightLine = ROOT.TLine(0.92,0.,0.92,hsignal[histo_name].GetMaximum()*1.1)
            leftLine.SetLineColor(2)
            leftLine.SetLineStyle(2)
            leftLine.SetLineWidth(3)
            rightLine.SetLineColor(2)
            rightLine.SetLineStyle(2)
            rightLine.SetLineWidth(3)
            leftLine.Draw()
            rightLine.Draw()

        if isKAnalysis:
            hstack[histo_name].GetXaxis().SetLimits(0.65,1.1)
            leftLine  = ROOT.TLine(0.85,0.,0.85,hsignal[histo_name].GetMaximum()*1.1)
            rightLine = ROOT.TLine(0.96,0.,0.96,hsignal[histo_name].GetMaximum()*1.1)
            leftLine.SetLineColor(3)
            leftLine.SetLineStyle(2)
            leftLine.SetLineWidth(3)
            rightLine.SetLineColor(3)
            rightLine.SetLineStyle(2)
            rightLine.SetLineWidth(3)
            leftLine.Draw()
            rightLine.Draw()

    if histo_name == "h_firstTrkPt" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{1}} [GeV]")
        #hstack[histo_name].GetXaxis().SetLimits(15.,60.)

    if histo_name == "h_secondTrkPt" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{Trk_{2}} [GeV]")
        #hstack[histo_name].GetXaxis().SetLimits(5.,50.)

    if histo_name == "h_firstTrkEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{Trk_{1}}")
        hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

    if histo_name == "h_secondTrkEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{Trk_{2}}")
        hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

    if histo_name == "h_mesonEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{ditrk}")
        hstack[histo_name].GetXaxis().SetLimits(-2.5,2.5)

    if histo_name == "h_firstTrkPhi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi_{Trk_{1}} [rad]")

    if histo_name == "h_secondTrkPhi" :
        hstack[histo_name].GetXaxis().SetTitle("#phi_{Trk_{2}} [rad]")

    if histo_name == "h_mesonPt" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{meson} [GeV]")    

    if histo_name == "h_mesonIsoCh" :
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{meson}/(#Sigmap_{T}^{ch} + p_{T}^{meson})")    

    if histo_name == "h_trksDeltaR" :
        hstack[histo_name].GetXaxis().SetTitle("#DeltaR_{ditrk}")

    if histo_name == "h_nPhotons" :
        hstack[histo_name].GetXaxis().SetTitle("n.#gamma")
        hstack[histo_name].GetXaxis().SetLimits(-0.5,4.5)

    if histo_name == "h_photonEnergy" :
        hstack[histo_name].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
        #hstack[histo_name].GetXaxis().SetLimits(38.,160.)

    if histo_name == "h_photonEta" :
        hstack[histo_name].GetXaxis().SetTitle("#eta_{#gamma}")

    if histo_name == "h_nElectrons" :
        hstack[histo_name].GetXaxis().SetTitle("n.electrons")
        hstack[histo_name].GetXaxis().SetLimits(-0.5,2.5)

    if histo_name == "h_nMuons" :
        hstack[histo_name].GetXaxis().SetTitle("n.muons")
        hstack[histo_name].GetXaxis().SetLimits(-0.5,3.5)


    hstack[histo_name].Draw("SAMES,histo")
    if signal_magnify != 1: hsignal[histo_name].Scale(signal_magnify)  
    hsignal[histo_name].Draw("SAMES,histo")

    if plotOnlyData:
        if isPhiAnalysis: leg1.AddEntry(hsignal[histo_name],"#phi#gamma data","l")
        if isRhoAnalysis: leg1.AddEntry(hsignal[histo_name],"#rho#gamma data","l")
        if isKAnalysis: leg1.AddEntry(hsignal[histo_name],"K*#gamma data","l")

    elif not plotOnlyData:
        if isPhiAnalysis: leg1.AddEntry(hsignal[histo_name],"#phi#gamma signal","l")
        if isRhoAnalysis: leg1.AddEntry(hsignal[histo_name],"#rho#gamma signal","l")
        if isKAnalysis: leg1.AddEntry(hsignal[histo_name],"K*#gamma signal","l")

    #CMS_lumi.CMS_lumi(canvas[histo_name], iPeriod, iPos) #Print integrated lumi and energy information
    leg1.Draw()


    ################################################

    if plotOnlyData:
        if isPhiAnalysis and region_name == "SR": 
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Phi/Data/SR/"
        elif isPhiAnalysis and region_name == "CR":
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Phi/Data/CR/"
        elif not isPhiAnalysis and region_name == "SR" :
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Rho/Data/SR/"
        elif not isPhiAnalysis and region_name == "CR" :
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Rho/Data/CR/"
            
    elif not plotOnlyData:
        if isWideRange:
            if isPhiAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HPhi/Signal/WideRange/"
            if isPhiAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZPhi_MadGraph/Signal/WideRange/"

            if isRhoAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HRho/Signal/WideRange/"
            if isRhoAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZRho_MadGraph/Signal/WideRange/"

            if isKAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HK/Signal/WideRange/"
            if isKAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZK/Signal/WideRange/"
        else:
            if isPhiAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HPhi/Signal/NarrowRange/"
            if isPhiAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZPhi_MadGraph/Signal/NarrowRange/"

            if isRhoAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HRho/Signal/NarrowRange/"
            if isRhoAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZRho_MadGraph/Signal/NarrowRange/"

            if isKAnalysis and isHAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/HK/Signal/NarrowRange/"
            if isKAnalysis and isZAnalysis: output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/ZK/Signal/NarrowRange/"

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
    canvas[histo_name].SaveAs(output_dir + histo_name + ".png")