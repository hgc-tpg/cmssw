import FWCore.ParameterSet.Config as cms

Stage2FileWriter = cms.EDAnalyzer('Stage2FileWriter',
  clusters = cms.untracked.InputTag("hgcalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"),
  outputFilename = cms.untracked.string("HGCS2OutputToL1TFile"),
  format = cms.untracked.string("EMP")
)