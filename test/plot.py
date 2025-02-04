import ROOT
import math
import copy
import sys
import tdrstyle, CMS_lumi
from ROOT import gROOT

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

signal_magnify = int(sys.argv[1])
CR_magnify = 1. 

plotOnlyData = int(sys.argv[2])

isPhi = int(sys.argv[3]) #note that in python true = 1 and false = 0
print( "#############################")
print( "is Phi = ",isPhi)
print( "#############################")

inputname = sys.argv[4] #CR or SR


list_inputfiles = []
for filename in sys.argv[5:]:
    list_inputfiles.append(filename)


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
signalfile = ROOT.TFile(sys.argv[5])
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



for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    #sample_name = (filename.split("_")[2])[:-5] 
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)

        print("histo_name = ",histo_name)
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

        if isPhi:
            histo_container[-1].SetLineColor(9)   #blue, 2 for red
        else:
            histo_container[-1].SetLineColor(46)   #blue, 2 for red

        histo_container[-1].SetLineWidth(4)   #kind of thick
        histo_container[-1].Scale(1./histo_container[-1].GetEntries()) #normalize to 1
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
    leg1 = ROOT.TLegend(0.65,0.74,0.87,0.97) #right positioning
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
        if isPhi:
            hstack[histo_name].GetXaxis().SetLimits(1.00,1.042)
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

        else:
            hstack[histo_name].GetXaxis().SetLimits(0.5,1.)
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


    hstack[histo_name].Draw("SAME,histo")


    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)
  
    hsignal[histo_name].Draw("SAME,hist")
    if plotOnlyData:
        if isPhi:
            leg1.AddEntry(hsignal[histo_name],"#phi#gamma data","l")
        elif not isPhi:
            leg1.AddEntry(hsignal[histo_name],"#rho#gamma data","l")
    elif not plotOnlyData:
        if isPhi:
            leg1.AddEntry(hsignal[histo_name],"#phi#gamma signal","l")
        elif not isPhi:
            leg1.AddEntry(hsignal[histo_name],"#rho#gamma signal","l")

    #CMS_lumi.CMS_lumi(canvas[histo_name], iPeriod, iPos) #Print integrated lumi and energy information
    leg1.Draw()


    ################################################

    if plotOnlyData:
        if isPhi and inputname == "SR": 
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Phi/Data/SR/"
        elif isPhi and inputname == "CR":
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Phi/Data/CR/"
        elif not isPhi and inputname == "SR" :
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Rho/Data/SR/"
        elif not isPhi and inputname == "CR" :
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Rho/Data/CR/"
            
    elif not plotOnlyData:
        if isPhi: 
            output_dir = "/eos/user/e/eferrand/www/eferrand/ZMesonGamma/Phi/Signal/"
        elif not isPhi:
            output_dir = "/eos/user/g/gumoret/www/HZMesonGamma/Rho/Signal/"

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
    canvas[histo_name].SaveAs(output_dir + histo_name + ".png")
