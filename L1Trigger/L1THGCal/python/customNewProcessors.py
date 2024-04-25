import FWCore.ParameterSet.Config as cms
from L1Trigger.L1THGCal.hgcalBackendLayer2_fwClustering_cfi import layer2ClusteringFw_Params

def custom_clustering_standalone(process):
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.ProcessorName = cms.string('HGCalBackendLayer2Processor3DClusteringSA')
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.DistributionParameters = truncation_params
    process.l1tHGCalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.layer2FwClusteringParameters = layer2ClusteringFw_Params
    return process

def custom_tower_standalone(process):
    process.l1tHGCalTowerProducer.ProcessorParameters.ProcessorName = cms.string('HGCalTowerProcessorSA')
    return process
