import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('DIGI',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
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
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('/store/relval/CMSSW_13_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/125X_mcRun4_realistic_v5_2026D88noPU-v1/2590000/02fb5828-a118-49d3-8fcb-142debe91c20.root'),
       inputCommands=cms.untracked.vstring(
           'keep *',
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
    fileName = cms.string("ntuple.root")
    )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# load HGCAL TPG simulation
#process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitivesNew_cff')

# Use new Stage 1 processor
from L1Trigger.L1THGCal.customNewProcessors import custom_stage1_truncation
process = custom_stage1_truncation(process)

# Switch to latest trigger geometry containing information on links mapping
from L1Trigger.L1THGCal.customTriggerGeometry import custom_geometry_V11_Imp3
process = custom_geometry_V11_Imp3(process)

process.hgcl1tpg_step = cms.Path(process.L1THGCalTriggerPrimitives)


# load ntuplizer and custom to use collections from Stag1 truncation processor
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.customNtuples import custom_ntuples_stage1_truncation
from L1Trigger.L1THGCalUtilities.customNtuples import custom_ntuples_layer1_truncation
process = custom_ntuples_stage1_truncation(process)
#process = custom_ntuples_layer1_truncation(process)
process.ntuple_step = cms.Path(process.L1THGCalTriggerNtuples)

# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.ntuple_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
