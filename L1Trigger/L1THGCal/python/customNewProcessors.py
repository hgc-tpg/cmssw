import FWCore.ParameterSet.Config as cms
from L1Trigger.L1THGCal.l1tHGCalBackEndLayer1Producer_cfi import truncation_params
from L1Trigger.L1THGCal.hgcalBackendLayer2_fwClustering_cfi import layer2ClusteringFw_Params
from L1Trigger.L1THGCal.l1tHGCalBackEndLayer2Producer_cfi import stage2_emulator_proc

def custom_clustering_standalone(process):
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.ProcessorName = cms.string('HGCalBackendLayer2Processor3DClusteringSA')
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.DistributionParameters = truncation_params
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.layer2FwClusteringParameters = layer2ClusteringFw_Params
    return process

def custom_tower_standalone(process):
    process.l1tHGCalTowerProducer.ProcessorParameters.ProcessorName = cms.string('HGCalTowerProcessorSA')
    return process

def custom_stage2_emulator(process):
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters = stage2_emulator_proc.clone()
    return process
