import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('inputFile', '', VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, 'D127 GEN-SIM-DIGI-RAW input file, prefixed with file:')
options.outputFile = 'stage2-emulator-v19.root'
options.maxEvents = 10
options.parseArguments()

if not options.inputFile:
    raise RuntimeError('Pass inputFile=file:/path/to/D127_GEN-SIM-DIGI-RAW.root')

from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9
process = cms.Process('STAGE2EMU', Phase2C26I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D127Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D127_cff')
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

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(options.inputFile))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
from L1Trigger.L1THGCal.customNewProcessors import custom_stage2_comparison
process = custom_stage2_comparison(process)

process.TFileService = cms.Service('TFileService', fileName=cms.string(options.outputFile))
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cfi import (
    ntuple_clusters,
    ntuple_event,
    ntuple_gen,
    ntuple_multiclusters,
)
input_clusters = ntuple_clusters.clone(
    Prefix=cms.untracked.string('cl2d'),
    Multiclusters=cms.InputTag('l1tHGCalBackEndLayer2Producer',
                               'HGCalBackendLayer2ProcessorStage2Emulator'),
)
reference_multiclusters = ntuple_multiclusters.clone(
    Prefix=cms.untracked.string('refcl3d'),
    Multiclusters=cms.InputTag('l1tHGCalBackEndLayer2ProducerReference',
                               'HGCalBackendLayer2Processor3DClustering'),
)
semi_emulator_multiclusters = ntuple_multiclusters.clone(
    Prefix=cms.untracked.string('semicl3d'),
    Multiclusters=cms.InputTag('l1tHGCalBackEndLayer2ProducerSemiEmulator',
                               'HGCalBackendLayer2ProcessorSemiEmulator'),
)
stage2_emulator_multiclusters = ntuple_multiclusters.clone(
    Prefix=cms.untracked.string('s2cl3d'),
    Multiclusters=cms.InputTag('l1tHGCalBackEndLayer2Producer',
                               'HGCalBackendLayer2ProcessorStage2Emulator')
)
process.l1tHGCalTriggerNtuplizer.Ntuples = cms.VPSet(
    ntuple_event,
    ntuple_gen.clone(MCEvent=cms.InputTag('generatorSmeared', '', 'GEN')),
    input_clusters,
    reference_multiclusters,
    semi_emulator_multiclusters,
    stage2_emulator_multiclusters,
)

process.tpg = cms.Path(process.L1THGCalTriggerPrimitives)
process.ntuple = cms.EndPath(process.l1tHGCalTriggerNtuplizer)
