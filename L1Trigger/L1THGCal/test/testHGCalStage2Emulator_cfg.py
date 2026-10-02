import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCal.l1tHGCalBackEndLayer2Producer_cfi import l1tHGCalBackEndLayer2Producer
from L1Trigger.L1THGCal.customNewProcessors import custom_stage2_emulator

process = cms.Process("STAGE2EMU")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))
process.source = cms.Source("EmptySource")

process.l1tHGCalBackEndLayer2Producer = l1tHGCalBackEndLayer2Producer.clone()
process = custom_stage2_emulator(process)
process.p = cms.Path(process.l1tHGCalBackEndLayer2Producer)
