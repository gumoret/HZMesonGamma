#######################################################################
#                                                                     #
# <<  Hi! I smuggle bad coded functions. But, ehy, they're free! >>   #
#                                                                     #
#######################################################################                                                                                                                                                       

import ROOT
import math
import os
from array import array
import time


#------- Arrays and reader for the BDT -------#

mesonIsoCh_array            = array('f', [0.])
mesonIsoNeu_array           = array('f', [0.])
cosTheta_array              = array('f', [0.])
mesonPt_array               = array('f', [0.])
mesonEta_array              = array('f', [0.])
photonEt_array              = array('f', [0.])
mesonGammaDeltaPhi_array    = array('f', [0.])
lxy_array                   = array('f', [0.])

reader = ROOT.TMVA.Reader("!Color")


#MCfromDATA correction for photons ###############################################################################
ph_ID_scale_name_2024  = "/eos/cms/store/group/phys_egamma/ScaleFactors/Data2024/PhoID/passingMVA122XV1wp80/egammaEffi.txt_EGM2D.root"
ph_ID_scale_file_2024  = ROOT.TFile(ph_ID_scale_name_2024)
ph_ID_scale_histo_2024 = ROOT.TH2F()
ph_ID_scale_histo_2024 = ph_ID_scale_file_2024.Get("EGamma_SF2D")
'''
ph_pixVeto_scale_name_2018  = "scale_factors/HasPix_SummaryPlot_UL18.root"
ph_pixVeto_scale_file_2018  = ROOT.TFile(ph_pixVeto_scale_name_2018)
ph_pixVeto_scale_histo_2018 = ROOT.TH2F()
ph_pixVeto_scale_histo_2018 = ph_pixVeto_scale_file_2018.Get("MVAID/SF_HasPix_MVAID")
'''
class Simplified_Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT,isD0=True):

                
        ###################################################################################
        #                                                                                 #
        #------------------------ Add BDT variables to the reader ------------------------#
        #                                                                                 #
        ###################################################################################
        if not isD0:
            #reader.AddVariable("mesonIsoCh",mesonIsoCh_array)
            #reader.AddVariable("mesonIsoNeu",mesonIsoNeu_array)
            reader.AddVariable("cosTheta",cosTheta_array)
            reader.AddVariable("mesonPt/bosonMass",mesonPt_array)#/bosonMass
            reader.AddVariable("mesonEta",mesonEta_array)
            reader.AddVariable("photonEt/bosonMass",photonEt_array)#/bosonMass
            #reader.AddVariable("mesonGammaDeltaPhi",mesonGammaDeltaPhi_array)
        elif isD0:
            reader.AddVariable("mesonPt/bosonMass",mesonPt_array)#/bosonMass
            reader.AddVariable("mesonEta",mesonEta_array)
            reader.AddVariable("photonEt/bosonMass",photonEt_array)#/bosonMass
            reader.AddVariable("lxy", lxy_array)



        if isBDT: reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    #Get BDT output function ###########################################################################################################

    def get_BDT_output(self,mesonPt,bosonMass,mesonEta,photonEt,lxy):#mesonIsoCh,cosTheta,mesonIsoNeu,bosonMass,photonEt,mesonGammaDeltaPhi
        
        if lxy is None:
            #mesonIsoCh_array[0]           = mesonIsoCh
            #mesonIsoNeu_array[0]          = mesonIsoNeu
            cosTheta_array[0]             = cosTheta
            mesonPt_array[0]              = mesonPt/bosonMass
            mesonEta_array[0]             = mesonEta
            photonEt_array[0]             = photonEt/bosonMass 
            #mesonGammaDeltaPhi_array[0]   = mesonGammaDeltaPhi
        if lxy is not None:
            mesonPt_array[0]              = mesonPt/bosonMass
            mesonEta_array[0]             = mesonEta
            photonEt_array[0]             = photonEt/bosonMass
            lxy_array[0]              = lxy

        return reader.EvaluateMVA("BDT")

    #Photon scaling function ###########################################################################################################
    def get_photon_scale(self,ph_pt, ph_eta):

        local_ph_pt = ph_pt
        if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
            local_ph_pt = 499.
        
        local_ph_eta = ph_eta
        if local_ph_eta >= 2.5:
            local_ph_eta = 2.49
        if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
            local_ph_eta = -2.5

        scale_factor_ID = ph_ID_scale_histo_2024.GetBinContent( ph_ID_scale_histo_2024.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2024.GetYaxis().FindBin(local_ph_pt) )
        ph_ID_err       = ph_ID_scale_histo_2024.GetBinError( ph_ID_scale_histo_2024.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2024.GetYaxis().FindBin(local_ph_pt) )

        #etaBin = 1 if abs(local_ph_eta) < 1.48 else 4
        #scale_factor_pixVeto = ph_pixVeto_scale_histo_2018.GetBinContent(etaBin)
        #ph_pixVeto_err       = ph_pixVeto_scale_histo_2018.GetBinError(etaBin)
        
        scale_factor = scale_factor_ID #* scale_factor_pixVeto
        tot_err      = ph_ID_err #math.sqrt( ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )

        return scale_factor, tot_err
        