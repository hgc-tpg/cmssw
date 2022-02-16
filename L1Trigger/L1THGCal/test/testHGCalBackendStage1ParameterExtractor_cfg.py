import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('SIM',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
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

'''
process.MessageLogger = cms.Service("MessageLogger",
                                    destinations   = cms.untracked.vstring('myDebugOutputFile'),
                                    myDebugOutputFile = cms.untracked.PSet(
                                        threshold = cms.untracked.string('DEBUG'),      #4
                                        default = cms.untracked.PSet(                   
                                            limit = cms.untracked.int32(-1) #5
                                        ),
                                        debugModules = cms.untracked.vstring(                            #6
                                            'hgcalbackendstage1parameterextractor')           #7
                                    )                                                                       #8
                                )
'''

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:junk.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("test_triggergeom.root")
    )



# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(10.01),
        MinPt = cms.double(9.99),
        PartID = cms.vint32(13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single electron pt 10'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)

process.mix.digitizers = cms.PSet(process.theDigitizersValid)

# load HGCAL TPG simulation
#process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitivesNew_cff')

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Eventually modify default geometry parameters
from L1Trigger.L1THGCal.customTriggerGeometry import custom_geometry_V11_Imp3
process = custom_geometry_V11_Imp3(process)

# Fetch stage 1 truncation parameters from correct config
from L1Trigger.L1THGCal.hgcalBackEndLayer1Producer_cfi import stage1truncation_proc

# ordered u/v coordinated of TCs in a module, for consistent TC index definition
ordered_tcu = [4, 5, 4, 3, 6, 7, 6, 5, 4, 5, 4, 3, 2, 3, 2, 1, 7, 6, 6, 7, 5, 4, 4, 5, 5, 4, 4, 5, 7, 6, 6, 7, 0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3]
ordered_tcv = [0, 1, 1, 0, 2, 3, 3, 2, 2, 3, 3, 2, 0, 1, 1, 0, 7, 7, 6, 6, 7, 7, 6, 6, 5, 5, 4, 4, 5, 5, 4, 4, 3, 2, 3, 4, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6]


process.hgcalbackendstage1parameterextractor = cms.EDAnalyzer(
    "HGCalBackendStage1ParameterExtractor",
    outJSONname = cms.string("test.json"),
    testedFpga = cms.int32(1),
    BackendStage1Params = stage1truncation_proc,
    TriggerGeometryParam = process.hgcalTriggerGeometryESProducer.TriggerGeometry,
    TCcoord_UV = cms.PSet(
        TCu = cms.vuint32(ordered_tcu),
        TCv = cms.vuint32(ordered_tcv)
    ),
    MetaData = cms.PSet(
        CMSSW_version = cms.string('CMSSW_12_3_0_pre4')
    ),
)

process.test_step = cms.Path(process.hgcalbackendstage1parameterextractor)

# Schedule definition
process.schedule = cms.Schedule(process.test_step,process.endjob_step)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
