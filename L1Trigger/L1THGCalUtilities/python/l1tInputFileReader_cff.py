import FWCore.ParameterSet.Config as cms

L1TInputFileReader = cms.EDProducer('L1TInputFileReader',
  files = cms.vstring(("input-emp-vu9p.001.txt")),
  format = cms.untracked.string("EMP")
)