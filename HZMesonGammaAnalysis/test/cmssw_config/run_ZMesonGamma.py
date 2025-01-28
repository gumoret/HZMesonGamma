import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the correct one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

################################################################################################################
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
#                       era='2018-Prompt')  
################################################################################################################



#Input source
if options.runningOnData:
   #process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v33') # OLD ONE : 102X_dataRun2_Sep2018ABC_v2##########
   input_path = '/store/data/Run2018B/Tau/MINIAOD/UL2018_MiniAODv2-v2/'#############


else:
   #process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_forMB_v1' )  # OLD ONE : '106X_upgrade2018_realistic_v15_L1v1'
   input_path = '/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/Zrhogamma/'


   #INPUT FILE LIST
   '''                                                                                                                                                                                                   
   For the given path, get the List of all files in the directory tree                                                                                                                               
   '''                                                                                                                                                                                                   
   def getListOfFiles(dirName):                                                                                                                                                                          
       # create a list of file and sub directories                                                                                                                                                       
       # names in the given directory                                                                                                                                                                    
       listOfFile = os.listdir(dirName)                                                                                                                                                                  
       allFiles = list()                                                                                                                                                                                 
       # Iterate over all the entries                                                                                                                                                                    
       for entry in listOfFile:                                                                                                                                                                          
           # Create full path                                                                                                                                                                            
           fullPath = os.path.join(dirName, entry)                                                                                                                                                       
           # If entry is a directory then get the list of files in this directory                                                                                                                        
           if os.path.isdir(fullPath):                                                                                                                                                                   
              allFiles = allFiles + getListOfFiles(fullPath)                                                                                                                                          
           else:                                                                                                                                                                                         
              allFiles.append("file:" + fullPath)                                                                                                                                                                 
                                                                                                                                                                                                          
       return allFiles     

   # Get the list of all files in directory tree at given path                                                                                                                                           
   listOfFiles = getListOfFiles(input_path)                                                                                                                                                          
   #print(listOfFiles)  

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring (listOfFiles), duplicateCheckMode = cms.untracked.string ('noDuplicateCheck'))


# Output file
if options.runningOnData:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("ZMesonGamma_Data.root"))
else:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("ZMesonGamma_Signal.root"))



process.load("HZMesonGammaAnalysis.HZMesonGamma.ZMesonGamma_cfi")
process.ZMesonGamma.runningOnData = options.runningOnData


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################


process.seq = cms.Path(process.ZMesonGamma)
process.schedule = cms.Schedule(process.seq)