import FWCore.ParameterSet.Config as cms 

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('DIGI',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    #  input = cms.untracked.int32(-1)
    input = cms.untracked.int32(10)
)

import glob
input_files = glob.glob('/home/llr/cms/sauvan/DATA_UPG/HGCAL/2110-2112_occupancies/ttbar_D86_12_1_0_pre5_PU200_20211101/*.root')
inputs = ['file:'+f for f in input_files]

# Input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(inputs),
        #  fileNames = cms.untracked.vstring('file:../../../../../data/ttbar_D86_12_1_X_PU200_20210920/Events_6491114_0.root'),
       inputCommands=cms.untracked.vstring(
           'keep *',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
           'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLBarrel_RECO',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLEndcap_RECO',
           'drop MTDTrackingRecHitedmNewDetSetVector_mtdTrackingRecHits__RECO',
           'drop BTLDetIdBTLSampleFTLDataFrameTsSorted_mix_FTLBarrel_HLT',
           'drop ETLDetIdETLSampleFTLDataFrameTsSorted_mix_FTLEndcap_HLT',
           )
       )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    #  fileName = cms.string("../../../../../data/ttbar_D86_12_1_0_pre5_PU200_20211101_ntuple_211104.root")
    fileName = cms.string("ntuple.root")
    )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

#  from SLHCUpgradeSimulations.Configuration.aging import customise_aging_3000
#  process = customise_aging_3000(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import agedHGCal
process = agedHGCal(process)
process.HGCAL_noise_fC.doseMap=cms.string("SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt")
process.HGCAL_noise_heback.doseMap=cms.string("SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt")

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
# Use latest V3 trigger geometry implementation
#  from L1Trigger.L1THGCal.customTriggerGeometry import custom_geometry_decentralized_V11
#  process = custom_geometry_decentralized_V11(process, implementation=2)

#  process.hgcalVFEProducer.ProcessorParameters.connectAllModules = True

process.hgcalConcentratorProducer.ProcessorParameters.threshold_silicon = 0.5
process.hgcalConcentratorProducer.ProcessorParameters.threshold_scintillator = 0.5

# Setup new digi for startup scenario
process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_ee.oldDigi = False
process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_hesi.oldDigi = False
process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_hesc.oldDigi = False
process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_ee.oldDigi = False
process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_hesi.oldDigi = False
process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_hesc.oldDigi = False
#
#  process.hgcalVFEProducer.ProcessorParameters.summationCfg.noiseThreshold = 0.
#
#  scaleByDoseFactor = 1. # 3000 fb-1
#  process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_ee.scaleByDoseFactor = scaleByDoseFactor
#  process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_hesi.scaleByDoseFactor = scaleByDoseFactor
#  process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_ee.scaleByDoseFactor = scaleByDoseFactor
#  process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_hesi.scaleByDoseFactor = scaleByDoseFactor
   
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)
   
# load ntuplizer
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
process.hgcalTriggerNtuplizer.Ntuples = [
        process.ntuple_event,
        #  process.ntuple_digis,
        process.ntuple_triggercells
        ]
#  process.ntuple_triggercells.FillSimEnergy = True
process.ntuple_step = cms.Path(process.hgcalTriggerNtuples)


# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.ntuple_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

