import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') #Could be not the correct one, but should contain the one without "condDBv2"

from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#inputFiles={'file:/eos/user/p/pellicci/MesonGamma_root/2023/Zrhogamma_miniAOD/Zrhogamma_2018UL_11.root',
            #'file:/eos/user/p/pellicci/MesonGamma_root/2023/Zrhogamma_miniAOD/Zrhogamma_2018UL_0.root'}

input_path = '/eos/user/p/pellicci/MesonGamma_root/2024/Hphigamma_EvtGen_miniAOD/'


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

listOfFiles = getListOfFiles(input_path)                                                                                                                                                              
##print(listOfFiles) 

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring (listOfFiles), duplicateCheckMode = cms.untracked.string ('noDuplicateCheck'))   


#Output file
process.TFileService = cms.Service("TFileService", fileName = cms.string("ZMesonGammaGen_output.root"))
process.ZMesonGammaGen = cms.EDAnalyzer('ZMesonGammaGen') ##nome
process.p = cms.Path(process.ZMesonGammaGen) ##nome
