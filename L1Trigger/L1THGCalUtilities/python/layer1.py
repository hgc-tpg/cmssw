import FWCore.ParameterSet.Config as cms
from L1Trigger.L1THGCal.l1tHGCalBackEndLayer1Producer_cfi import layer1truncationfw_proc, \
    truncationfw_params, \
    layer1truncation_proc, \
    truncation_params

class RozBinTruncation(object):
    def __init__(self,
            maxTcsPerBin=truncation_params.maxTcsPerBin):
        self.processor = layer1truncation_proc.clone(
                truncation_parameters=truncation_params.clone(
                    maxTcsPerBin=maxTcsPerBin
                    )
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalBackEndLayer1Producer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                ProcessorParameters = self.processor,
                BypassBackendMapping = cms.bool(False)
                )
        return producer

class RozBinTruncationFw(object):
    def __init__(self,
            maxTcsPerBin=truncationfw_params.maxTcsPerBin,
            doTruncation=truncationfw_params.doTruncation):
        self.processor = layer1truncationfw_proc.clone(
                truncation_parameters=truncationfw_params.clone(
                maxTcsPerBin=maxTcsPerBin,
                doTruncation=doTruncation
                )
        )

    def __call__(self,process,inputs):
        producer = process.l1tHGCalBackEndLayer1Producer.clone(
            InputTriggerCells = cms.InputTag(inputs),
            ProcessorParameters = self.processor,
            BypassBackendMapping = cms.bool(False)
        )
        return producer
