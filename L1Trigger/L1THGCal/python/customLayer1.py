import FWCore.ParameterSet.Config as cms
from L1Trigger.L1THGCal.l1tHGCalBackEndLayer1Producer_cfi import layer1truncation_proc, layer1truncationfw_proc, truncation_params, truncationfw_params

def custom_layer1_truncation(process):
    parameters = layer1truncation_proc.clone()
    process.l1tHGCalBackEndLayer1Producer.ProcessorParameters = parameters
    process.l1tHGCalBackEndLayer1Producer.BypassBackendMapping = cms.bool(False)
    process.l1tHGCalBackEndLayer2Producer.InputCluster = cms.InputTag('l1tHGCalBackEndLayer1Producer:HGCalBackendLayer1ProcessorTruncation')
    process.l1tHGCalTowerProducer.InputTriggerCells = cms.InputTag('l1tHGCalBackEndLayer1Producer:HGCalBackendLayer1ProcessorTruncation')
    return process

def custom_layer1_truncationfw(process):
    parameters = layer1truncationfw_proc.clone()
    process.l1tHGCalBackEndLayer1Producer.ProcessorParameters = parameters
    process.l1tHGCalBackEndLayer1Producer.BypassBackendMapping = cms.bool(False)
    process.l1tHGCalBackEndLayer2Producer.InputCluster = cms.InputTag('l1tHGCalBackEndLayer1Producer:HGCalBackendLayer1ProcessorTruncationFw')
    process.l1tHGCalTowerProducer.InputTriggerCells = cms.InputTag('l1tHGCalBackEndLayer1Producer:HGCalBackendLayer1ProcessorTruncationFw')
    return process

